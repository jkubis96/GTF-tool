packages <- c(
  "stringr",
  "dplyr",
  "doSNOW",
  "foreach",
  "parallel",
  "doParallel",
  "data.table",
  "readr"
)

installed <- packages %in% installed.packages()[, "Package"]
if (any(!installed)) {
  install.packages(packages[!installed], dependencies = TRUE)
}

invisible(lapply(packages, library, character.only = TRUE))


source("../R/GTF.tool.R")

GTF <- load_annotation("test_anno.gtf")


GTF2 <- create_GTF_df(GTF, optimize = TRUE, shift = 100000)


GTF5 <- add_UTR(GTF2, five_prime_utr = 300, three_prime_utr = 800)

clean_df <- GTF5 %>%
  group_by(transcript_id, annotationType) %>%
  mutate(has_jbio = any(source == "JBIO-predicted")) %>%
  filter(!(has_jbio & source != "JBIO-predicted")) %>%
  ungroup() %>%
  select(-has_jbio)



clean_df_tmp <- clean_df[clean_df$gene_id == "ENSG00000187634.13", ]

GTF_tmp <- GTF5[GTF5$gene_id == "ENSG00000187634.13", ]

ref_flat <- refflat_create(GTF2)
