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

GTF.tool provides advanced functionalities for processing genomic annotation files in GTF (Gene Transfer Format) and GFF formats. The package is designed to handle large datasets from sources like NCBI, Ensembl, and GENCODE. 

It includes tools for loading, filtering, sorting, optimizing, and extending annotation data. Key features include:

- Efficient loading and filtering of annotation files based on genetic elements (e.g., genes, exons, CDS).
- Sorting annotations by chromosome, strand, and genomic position.
- Normalizing and optimizing GTF data by resolving naming conflicts and filling missing annotations.
- Creating CDS (coding sequence) regions based on exon coordinates grouped by transcript.
- Inferring intron regions by identifying gaps between exons within the same transcript.

The package leverages parallel processing for high performance, making it ideal for large-scale genomic analyses.


#### Installation

```
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install the package from URL with dependencies
devtools::install_url("https://github.com/jkubis96/GTF-tool/raw/refs/heads/extended_version/packages/GTF.tool_0.1.0.tar.gz", dependencies = TRUE)
```


#### Loading

```
library(GTF.tool)
```


#### Documentation

* [GTF-tool](https://jkubis96.github.io/GTF-tool/)

<br />


#### Example in R

```
GTF <- load_annotation('annotation.gtf')

GTF2 <- create_GTF_df(GTF, optimize = TRUE, shift = 100000)
 
GTF3 <- add_CDS(GTF2)
 
GTF4 <- add_introns(GTF3)
 
GTF5 <- add_UTR(GTF4, five_prime_utr = 400 , three_prime_utr = 1000)

GTF6 <- create_full_GTF(GTF5)

write.table(GTF6, 'new_annotation.gtf', quote = F, sep = '\t', col.names = F, row.names = F)

GTF7 <- create_reduced_GTF(GTF5)

write.table(GTF7, 'reduced_new_annotation.gtf', quote = F, sep = '\t', col.names = T, row.names = F)

GTF8 <- refflat_create(GTF5)

write.table(GTF8, 'annotation.refflat', quote = F, sep = '\t', col.names = F, row.names = F)
```

<br />

<br />



#### Have fun JBS©