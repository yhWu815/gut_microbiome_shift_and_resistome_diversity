library(tidyverse)
#import 16S sequencing report
# Four-12.12-4.R1.fq.gz
amplicon_report <- read_tsv("./original_data/16S_pacbio_report.tsv")
#average num_seqs
avg_num_seqs <- amplicon_report %>% 
  summarise(mean_num_seqs = mean(num_seqs, na.rm = T))
# 37015.92
#average length
avg_length <- amplicon_report %>% 
  summarise(
    total_len = sum(sum_len),
    total_num = sum(num_seqs),
    total_avg_len = total_len / total_num
  )
# 1504.175

#import shotgun metagenomic raw data sequencing report
shotgun_report <- read_tsv("./original_data/meta_illumina_report_raw.tsv")
num_Gb_per_sample <- shotgun_report %>% 
  mutate(sum_len_Gb = sum_len / 1000000000) %>% 
  summarise(mean_sum_len_Gb = sum(sum_len_Gb) / 4)
# 12.46312
avg_num_seqs_per_sample <- shotgun_report %>% 
  summarise(mean_num_seqs = 2 * mean(num_seqs, na.rm = T))
