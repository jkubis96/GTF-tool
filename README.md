### GTF-tool - R library for advanced genome annotation file adjustment

<br />



<p align="right">
<img  src="https://github.com/jkubis96/Logos/blob/main/logos/jbs_current.png?raw=true" alt="drawing" width="180" />
</p>


### Author: Jakub Kubiś

<div align="left">
 Institute of Bioorganic Chemistry<br />
 Polish Academy of Sciences<br />



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
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Install the package from URL with dependencies
remotes::install_url("https://github.com/jkubis96/GTF-tool/raw/refs/heads/main/packages/GTF.tool_0.1.3.tar.gz", dependencies=TRUE)
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
GTF <- load_annotation("GTF.tool/tests/test_anno.gtf")
 
GTF2 <- create_GTF_df(GTF, optimize = TRUE, shift = 100000)
 
GTF3 <- add_CDS(GTF2)
  
GTF4 <- add_introns(GTF3)
 
GTF5 <- add_UTR(GTF4, five_prime_utr = 300, three_prime_utr = 800)

GTF6 <- create_full_GTF(GTF5)
 
write.table(GTF6, 'new_annotation.gtf', quote = F, sep = '\t', col.names = F, row.names = F)

GTF7 <- create_reduced_GTF(GTF5)

write.table(GTF7, 'reduced_new_annotation.gtf', quote = F, sep = '\t', col.names = T, row.names = F)

ref_flat <- refflat_create(GTF5)
 
write.table(ref_flat, 'annotation.refflat', quote = F, sep = '\t', col.names = F, row.names = F)
```

<br />

<br />



#### Have fun JBS©