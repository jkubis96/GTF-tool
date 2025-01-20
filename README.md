### GTF-tool - R library for advanced genome annotation file adjustment

<br />



<p align="right">
<img  src="https://github.com/jkubis96/Logos/blob/main/logos/jbs_current.png?raw=true" alt="drawing" width="180" />
</p>


### Author: Jakub Kubiś

<div align="left">
 Institute of Bioorganic Chemistry<br />
 Polish Academy of Sciences<br />
 Department of Molecular Neurobiology<br />



<br />


## Description

GTF-tool provides advanced functionalities for processing genomic annotation files in GTF (Gene Transfer Format), GFF3, and GFF2 formats. The package is designed to handle large datasets from sources like NCBI, Ensembl, and GENCODE. 

It includes tools for loading, filtering, sorting, optimizing, and extending annotation data. Key features include:

- Efficient loading and filtering of annotation files based on genetic elements (e.g., genes, exons, CDS).
- Sorting annotations by chromosome, strand, and genomic position.
- Normalizing and optimizing GTF data by resolving naming conflicts and filling missing annotations.
- Predicting and adding 5' and 3' untranslated regions (UTRs) dynamically based on genomic proximity and neighboring elements.

The package leverages parallel processing for high performance, making it ideal 
for large-scale genomic analyses.


#### Installation

```
install.packages("https://github.com/jkubis96/GTF-tool/raw/refs/heads/main/packages/GTF.tool_0.1.0.tar.gz", repos = NULL, type = "source")
```


#### Loading

```
library(GTF.tool)
```


#### Documentation

* [GTF-tool](https://jkubis96.github.io/GTF-tool/)

<br />

<br />


#### Have fun JBS©