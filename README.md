# Data, Code and Figures

## Operating System Requirement

- Ubuntu 20.04 LTS (Prefer)

- Windows Subsystem Linux: Ubuntu 20.04 LTS

## Mafft and Fasttree version

- `mafft`: MAFFT v7.453 (2019/Nov/8)

- `fasttree`: FastTree 2.1.11 Double precision (No SSE3)

How to install?

```
sudo apt install mafft
sudo apt install fasttree
```

## Python virtual environment: Conda environment

How to install?

`conda env create -n analysis_16s --file analysis_16s.yml`

The details of conda config files are listed below.

```
name: analysis_16s
channels:
  - bioconda
  - conda-forge
  - defaults
dependencies:
  - biopython=1.79
  - fastp=0.23.4
  - fastqc=0.12.1
  - jupyter_contrib_nbextensions=0.7.0
  - lefse=1.1.2
  - matplotlib=3.2.2
  - nose=1.3.7
  - notebook=6.1.1
  - numpy=1.21.6
  - openpyxl=3.0.10
  - pandas=1.3.5
  - pyfastx=1.1.0
  - python=3.7.12
  - requests=2.31.0
  - scikit-bio=0.5.7
  - scikit-learn=1.0.2
  - scipy=1.7.3
  - seaborn=0.11.2
  - statannotations=0.6.0
  - unifrac=1.1.1
```

## Data files

- **Table_otu.raw.txt.xls**: OTU table processed by Uparse(97% similarity) on the [Magigene Cloud Platform](http://cloud.magigene.com/). Chloroplasts, mitochondria, Unclassified, archaea otus were removed.

- **otus.fa**: Representative otu sequences processed by Uparse(97% similarity) on the [Magigene Cloud Platform](http://cloud.magigene.com/). Chloroplasts, mitochondria, Unclassified, archaea otus were removed.

- **otus_aligned.fasta**: Aligned otu sequences by mafft. 

  `Bash: mafft otus.fa > otus_aligned.fa
`
- **otus_fasttree.tre**: Construct phylogenic tree from `otus_aligned.fa`. 

  `Bash: fasttree -nt otus_aligned.fa > otus_fasttree.tre`

- **otu_table_samples_rearranged.tsv**: rarefacted otu table of 24 samples as the input of biom package.

- **otu_table_samples.biom**: BIOM-Format otu table for Unifrac Distance Matrix. 

- **otu_table_individuals_rearranged_merged.tsv**: rarefacted otu table of 8 subgroups as the input of biom package.

- **otu_table_individuals.biom**: BIOM-Format otu table for Unifrac Distance Matrix.

- **rpkm.type.txt**: Abundance of Antibiotic-resistant genes (ARGs) grouped by types processed by [ARGs-OAP v3.0](https://github.com/xinehc/args_oap) 

- **rpkm.subtype.txt**: ARGs grouped by subtypes processed by [ARGs-OAP v3.0](https://github.com/xinehc/args_oap)

- **rpkm.genes.txt**: ARGs grouped by genes processed by [ARGs-OAP v3.0](https://github.com/xinehc/args_oap)

- **lefse_results.txt.xls**: lefse analysis result produced by [Magigene Cloud Platform](http://cloud.magigene.com/).

## Scripts and jupyter notebook files

- **data_analysis_visualization.ipynb**: Script for Relative Abundance, Microbial Diversity analysis and Data visualization.

- **SpeciesAbundance.py**: Script for calculating the species abundance for given otu table at given genus

- **BetaDiversity.py**: Script for calculating the distance matrix and PCoA analysis.

- **args_oap.sh**: Shell script for ARGs analysis.

## Conda Environment Configuration file

- **analysis_16s.yml**: conda environment config file with all necessary packages.

## Figures

- **LDA.plot**: LDA bar plot produced by [Magigene Cloud Platform](http://cloud.magigene.com/).

- **lefse_cladogram_plot.pdf**: LEfSe cladogram plot produced by [Magigene Cloud Platform](http://cloud.magigene.com/).
