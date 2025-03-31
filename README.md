# Data, Code and Figures

## Operating System Requirement

- Ubuntu 20.04 LTS (Prefer)

- Windows Subsystem Linux: Ubuntu 20.04 LTS


## Files in the `original_data` folder

- **Table_otu.raw.tsv**: OTU table processed by Uparse(97% similarity) on the [Magigene Cloud Platform](http://cloud.magigene.com/). Chloroplasts, mitochondria, Unclassified, archaea otus were removed.

- **otus.fa**: Representative otu sequences processed by Uparse(97% similarity) on the [Magigene Cloud Platform](http://cloud.magigene.com/). Chloroplasts, mitochondria, Unclassified, archaea otus were removed.

- **rpkm.type.txt**: Abundance of Antibiotic-resistant genes (ARGs) grouped by types processed by [ARGs-OAP v3.0](https://github.com/xinehc/args_oap) 

- **rpkm.subtype.txt**: ARGs grouped by subtypes processed by [ARGs-OAP v3.0](https://github.com/xinehc/args_oap)

- **rpkm.genes.txt**: ARGs grouped by genes processed by [ARGs-OAP v3.0](https://github.com/xinehc/args_oap)

- **lefse_results.txt.xls**: lefse analysis result produced by [Magigene Cloud Platform](http://cloud.magigene.com/).

- **args_oap.sh**: Shell script for ARGs analysis.

- **LDA.plot**: LDA bar plot produced by [Magigene Cloud Platform](http://cloud.magigene.com/).

- **lefse_cladogram_plot.pdf**: LEfSe cladogram plot produced by [Magigene Cloud Platform](http://cloud.magigene.com/).

- **metadata.csv**: metadata of five pandas

## Code files

- **Figure1.R**: R script for Figure 1

- **multidrug_polymyxin.ipynb**: Python notebook for multidrug and polymyxin ARG subtype heatmap plot.

## Data and Figure files

- **Figure1.pdf & Figure1.svg**: Figure 1

- **FigureS1.pdf & FigureS1.svg**: Figure S1

- **Table1.xlsx**: Table 1

- **Table2.xlsx**: Table 2
