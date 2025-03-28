library(Biostrings)
library(tidyverse)
library(ggplot2)
library(readxl)
library(phyloseq)
library(vegan)
library(reshape2)
library(ggpubr)
library(patchwork)
library(pheatmap)
library(ggplotify)
# Step 1: Create Phyloseq Objects
#read otu table
otu_table_raw <- read_tsv("./original_data/Table_otu_raw.tsv")
#extract taxonomy column
taxonomy_raw <- otu_table_raw[, "taxonomy", drop = FALSE]
#convert otu_table to matrix
otu_table_processed <- otu_table_raw
otu_table_processed <- otu_table_processed %>% remove_rownames %>%
  column_to_rownames(var = "#OTU ID")
otu_table_matrix <- as.matrix(otu_table_processed)
otu_table_matrix <- apply(otu_table_matrix, 2, as.numeric)
rownames(otu_table_matrix) <- rownames(otu_table_processed)
#parse taxonomy
taxonomy_matrix <- taxonomy_raw %>%
  mutate(taxonomy = gsub("k:|p:|c:|o:|f:|g:|s:", "", taxonomy)) %>%  # 移除前缀
  separate(taxonomy, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ",", fill = "right") %>%
  as.matrix()
rownames(taxonomy_matrix) <- rownames(otu_table_processed)
common_otus <- intersect(rownames(otu_table_matrix), rownames(taxonomy_matrix))
#read metadata
metadata <- read_csv("./original_data/metadata.csv")
head(metadata)
sample_data_phy <- sample_data(metadata)
rownames(sample_data_phy) <- metadata$SampleID
#read otu sequences
otu_sequences <- readDNAStringSet("./original_data/otus.fa")
names(otu_sequences) <- names(otu_sequences)

# create phyloseq object
otu_table_phy <- otu_table(otu_table_matrix, taxa_are_rows = TRUE)
tax_table_phy <- tax_table(taxonomy_matrix)

physeq <- phyloseq(otu_table_phy, tax_table_phy, 
                   sample_data_phy, otu_sequences)

# create a directory to store subplots & table
save_plot_dir <- "r-subplots"
if (!dir.exists(save_plot_dir)) {
  dir.create(save_Plot_dir)
}
save_data_dir <- "r-data"
if (!dir.exists(save_data_dir)) {
  dir.create(save_data_dir)
}
# Step 2: Rarefaction
set.seed(1234)
rarefied_ps <- rarefy_even_depth(physeq, replace = TRUE)

# Step 3: Compute Alpha Diversity
alpha_div <- estimate_richness(rarefied_ps, measures = NULL)
# Merge alpha diversity with metadata
metadata <- sample_data(rarefied_ps) %>% as_tibble(rownames = "SampleID")
alpha_div <- alpha_div %>% as_tibble(rownames = "SampleID")
alpha_div_meta <- left_join(alpha_div, metadata, by = "SampleID")

# Alpha diversity metrics needed for analysis
alpha_metrics <- c("Observed", "Chao1", "ACE", "Shannon", 
                   "Simpson", "InvSimpson", "Fisher")

# variables
grouping_vars <- c("Panda", "Group1", "Group2")

# create empty lists to store plots
significant_plots <- list()
nonsignificant_plots <- list()

for (metric in alpha_metrics) {
  for (group in grouping_vars) {
    # kruskal-wallis tests
    kruskal_res <- kruskal.test(as.formula(paste(metric, "~", group)), 
                                data = alpha_div_meta)
    
    # extract p-value
    p_value <- kruskal_res$p.value
    
    # draw boxplot
    p <- ggplot(alpha_div_meta, aes_string(x = group, y = metric, fill = group)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.5) +
      labs(title = paste(metric, "by", group, "\nP-value", signif(p_value, 3)),
           x = group, y = metric) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # catagorize them according to signif
    if (p_value < 0.05) {
      significant_plots[[paste(metric, group, sep="_")]] <- p
    } else {
      nonsignificant_plots[[paste(metric, group, sep="_")]] <- p
    }
  }
}

# merge boxplots which are significant
if (length(significant_plots) > 0) {
  sig_plot <- wrap_plots(significant_plots) + plot_annotation(title = "Significant Alpha Diversity Boxplots (p < 0.05)")
  print(sig_plot)
  sig_plot_path <- "./r-subplots/sig_plots.pdf"
  ggsave(filename = sig_plot_path, plot = sig_plot, width = 12, height = 10, dpi = 600)
}

# merge boxplots which are not significant
if (length(nonsignificant_plots) > 0) {
  nonsig_plot <- wrap_plots(nonsignificant_plots) + plot_annotation(title = "Non-Significant Alpha Diversity Boxplots (p >= 0.05)")
  print(nonsig_plot)
  nonsig_plot_path <- "./r-subplots/nonsig_plots.pdf"
  ggsave(filename = nonsig_plot_path, plot = nonsig_plot, width = 14, height = 12, dpi = 600)
}



# Step 4: visualize the shannon diversity
# unify the theme
base_theme <- theme_minimal(base_family = "Arial", base_size = 8) +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        strip.text = element_text(size = 8, face = "bold")
  )
color_plate <- c("#9BBBE1", "#BABABC", "#C8C2E4", "#E0B77F", "#F09BA0")
group2_colors <- c("#9BBBE1", "#C8C2E4", "#F09BA0")
shannon_group2 <- ggplot(data = alpha_div_meta, aes(x = Group2, y = Shannon, fill = Group2)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  stat_compare_means(method="kruskal.test", label = "p.format", label.x.npc = 0.1, label.y.npc = 0.95) +
  scale_fill_manual(values = group2_colors) +
  labs(x = "Group", y = "Shannon Index") +
  theme_classic2() +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
print(shannon_group2) #significant, fig1
ggsave(filename = "./r-subplots/shannon_group2.svg", plot = shannon_group2, dpi = 1200)
ggsave(filename = "./r-subplots/shannon_group2.pdf", plot = shannon_group2, dpi = 1200)
# Step5. Compute Beta Diversity
bray_dist <- distance(rarefied_ps, method = "bray")
#PCoA analysis
pcoa_bray <- ordinate(rarefied_ps, method = "PCoA", distance = bray_dist)
##bray curtis ~ Group2
adonis_bray_Group2 <- adonis2(bray_dist ~ Group2, data = metadata)
adonis_bray_Group2_r2_value <- round(adonis_bray_Group2$R2[1], 3)
adonis_bray_Group2_p_value  <- round(adonis_bray_Group2$`Pr(>F)`[1], 3)
bray_Group2_pcoa <- plot_ordination(rarefied_ps, pcoa_bray, 
                                  type = "samples", color = "Group2") +
  scale_color_manual(values = group2_colors) +
  stat_ellipse(type = "t", level = 0.95) +
  annotate("text", x = 0.3, y = 1,
           label = paste("R2 =", adonis_bray_Group2_r2_value, "p =", adonis_bray_Group2_p_value)) +
  theme_classic2()
print(bray_Group2_pcoa)
ggsave(filename = "./r-subplots/bray_Group2_pcoa.svg", plot = bray_Group2_pcoa, dpi = 1200)
ggsave(filename = "./r-subplots/bray_Group2_pcoa.pdf", plot = bray_Group2_pcoa, dpi = 1200)
## heatmap
bray_matrix <- as.matrix(bray_dist)
heatmap_metadata <- data.frame(sample_data(rarefied_ps))
heatmap_metadata$SampleID <- rownames(heatmap_metadata)  # 确保行名为样本名
group_map <- heatmap_metadata %>% select(SampleID, Group1)
#compute average distance of each group1 variable
bray_group_matrix <- bray_matrix %>%
  as.data.frame() %>%
  rownames_to_column("Sample1") %>%
  gather(key = "Sample2", value = "Distance", -Sample1) %>%
  left_join(group_map, by = c("Sample1" = "SampleID")) %>%
  rename(Group1_Sample1 = Group1) %>%
  left_join(group_map, by = c("Sample2" = "SampleID")) %>%
  rename(Group1_Sample2 = Group1) %>%
  group_by(Group1_Sample1, Group1_Sample2) %>%
  summarise(Average_Distance = mean(Distance)) %>%
  spread(key = Group1_Sample2, value = Average_Distance)
#convert it to matrix
bray_group_matrix <- column_to_rownames(bray_group_matrix, var = "Group1_Sample1")
bray_group_matrix <- as.matrix(bray_group_matrix)
#plot group1 level heatmap
bray_group_heatmap <- pheatmap(bray_group_matrix, 
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = colorRampPalette(c("#9BBBE1", "white", "#F09BA0"))(50),
         display_numbers = FALSE)
gg_bray_group_heatmap <- as.ggplot(bray_group_heatmap)
ggsave(filename = "./r-subplots/bray_group_heatmap.svg", plot = gg_bray_group_heatmap, dpi = 1200)
ggsave(filename = "./r-subplots/bray_group_heatmap.pdf", plot = gg_bray_group_heatmap, dpi = 1200)
# Step 6. Species Stacked Bar Plot
tax_level <- "Genus"
ps_rel <- transform_sample_counts(rarefied_ps, function(x) x / sum(x))
# ASV to dataframe
asv_df <- psmelt(ps_rel)
# extract metadata and rearrange as Group1
species_metadata <- data.frame(sample_data(rarefied_ps))
species_metadata$SampleID <- rownames(species_metadata)
# make sure the order of Group1
species_metadata$Group1 <- factor(species_metadata$Group1, levels = c("S1", "S2", "S3", "E1", "E2", "E3", "E4", "E5"))

# compute relative abundance
genus_abundance <- asv_df %>%
  group_by(Group1, Genus) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# top 10 genus
top_genera <- genus_abundance %>%
  group_by(Genus) %>%
  summarise(TotalAbundance = sum(Abundance)) %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 10) %>%
  pull(Genus)

# "Others" categorize other genus as "Others"
genus_abundance <- genus_abundance %>%
  mutate(Genus = ifelse(Genus %in% top_genera, Genus, "Others"))

# compute ratio again
genus_abundance <- genus_abundance %>%
  group_by(Group1, Genus) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  group_by(Group1) %>%
  mutate(Percentage = Abundance / sum(Abundance) * 100)

# set color plate
species_color_plate <- c("#9BBBE1", "#C8C2E4", "#E0B77F", "#F09BA0",
                         "#A4D8A4", "#F4B400", "#FF6F61", "#B565A7", "#2E86C1", "#58D68D", "#BABACB")
names(species_color_plate) <- c(top_genera, "Others")

# plot stacked bar plot
genus_plot <- ggplot(genus_abundance, aes(x = Group1, y = Percentage, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +
  scale_y_continuous(expand = c(0.01, 0)) +
  scale_fill_manual(values = species_color_plate) +
  labs(x = "Group1", y = "Relative Abundance (%)", fill = "Genus") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        axis.title.x = element_blank())
print(genus_plot)
ggsave(filename = "./r-subplots/genus_plot.svg", plot = genus_plot, dpi = 1200)
ggsave(filename = "./r-subplots/genus_plot.pdf", plot = genus_plot, dpi = 1200)
# convert genus_abundance to long type
genus_abundance_wide <- genus_abundance %>%
  select(Group1, Genus, Percentage) %>%
  pivot_wider(names_from = Genus, values_from = Percentage, values_fill = list(Percentage = 0))
library(writexl)
write_xlsx(genus_abundance_wide, "./r-data/top10_other_abundance.xlsx")
## plot the stacked bar plot for the 24 samples
asv_df$Group1 <- factor(asv_df$Group1, 
                        levels = c("S1", "S2", "S3", "E1", "E2", "E3", "E4", "E5"))
all_genus_abundance <- asv_df %>%
  group_by(Sample, Group1, Genus) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(Percentage = Abundance / sum(Abundance) * 100)
top_10_genera <- all_genus_abundance %>%
  group_by(Genus) %>%
  summarise(TotalAbundance = sum(Percentage)) %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 10) %>%
  pull(Genus)
all_genus_abundance <- all_genus_abundance %>%
  mutate(Genus = ifelse(Genus %in% top_10_genera, Genus, "Others"))
all_genus_abundance <- all_genus_abundance %>%
  group_by(Sample, Group1, Genus) %>%
  summarise(Percentage = sum(Percentage), .groups = "drop")
all_genus_abundance <- all_genus_abundance %>%
  arrange(Group1, Sample)
order = c("S1.1st", "S1.2nd", "S1.3rd",
          "S2.1st", "S2.2nd", "S2.3rd",
          "S3.1st", "S3.2nd", "S3.3rd",
          "E1.1st", "E1.2nd", "E1.3rd",
          "E2.1st", "E2.2nd", "E2.3rd",
          "E3.1st", "E3.2nd", "E3.3rd",
          "E4.1st", "E4.2nd", "E4.3rd",
          "E5.1st", "E5.2nd", "E5.3rd")
new_all_genus_abundance <- all_genus_abundance %>%
  mutate(newSample = factor(Sample, levels = order))
all_sample_genus_plot <- ggplot(new_all_genus_abundance, aes(x = newSample, y = Percentage, fill = Genus)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = species_color_plate) +
    labs(x = "Sample", y = "Relative Abundance (%)", fill = "Genus") +
    scale_y_continuous(expand = c(0.01, 0)) +
    theme_classic() +
    theme(
      strip.background =  element_blank(),
      strip.placement = "outside",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(size = 12, face = "bold"),
      panel.spacing = unit(0.01, "lines"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank()) +
      facet_grid(~ Group1, scales = "free", switch = "both")
print(all_sample_genus_plot)
ggsave(filename = "./r-subplots/all_sample_genus_plot.svg", plot = all_sample_genus_plot, dpi = 1200, width = 10)
ggsave(filename = "./r-subplots/all_sample_genus_plot.pdf", plot = all_sample_genus_plot, dpi = 1200, width = 10)
# Step.7 ARG heatmap
arg_types = read.csv(file = "./original_data/rpkm.type.txt", sep="\t")
#rename macrolide-lincosamide-streptogramin as mls
#rename other_peptide_antibiotics as opa
arg_types[10, "type"] <- "mls"
arg_types[14, "type"] <- "opa"
new_col_names <- c("type", "C", "S1", "S3", "S2")
colnames(arg_types) <- new_col_names
rownames(arg_types) <- arg_types$type
arg_type_matrix <- as.matrix(arg_types[, -1])
arg_type_heatmap_column <- pheatmap(arg_type_matrix, scale = "none",
         cluster_rows = FALSE, cluster_cols = FALSE,
         cellwidth = 10, cellheight = 10,
         color = colorRampPalette(c( "#2E86C1", "white", "#F4B400"))(100),
         display_numbers = FALSE
         )
ggsave(filename = "./r-subplots/arg_type_heatmap_column.svg", plot = as.ggplot(arg_type_heatmap_column), 
       dpi = 1200)
ggsave(filename = "./r-subplots/arg_type_heatmap_column.pdf", plot = as.ggplot(arg_type_heatmap_column), 
       dpi = 1200)
arg_type_heatmap_row <- pheatmap(arg_type_matrix, scale = "none",
         cluster_rows = FALSE, cluster_cols = FALSE,
         cellwidth = 10, cellheight = 10,
         color = colorRampPalette(c( "#2E86C1", "white", "#F4B400"))(100),
         display_numbers = FALSE
)
ggsave(filename = "./r-subplots/arg_type_heatmap_row.svg", plot = as.ggplot(arg_type_heatmap_row), dpi = 1200)
ggsave(filename = "./r-subplots/arg_type_heatmap_row.pdf", plot = as.ggplot(arg_type_heatmap_row), dpi = 1200)

