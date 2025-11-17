#' Add CDS (Coding Sequence) Annotations to Genomic Data
#'
#' @description
#' This function generates coding sequence (CDS) regions for each transcript based on its associated exons.
#' It assumes that input data includes properly assigned exons and their transcript IDs. The output consists
#' of CDS regions suitable for further analysis or export in GTF/GFF format.
#'
#' @param input A data frame containing genomic data. The data frame should have the following columns:
#'   - `chr`: Chromosome identifier
#'   - `start`: Start position of the annotation
#'   - `end`: End position of the annotationdevtools::build()
#'   - `strand`: Strand information ('+' or '-')
#'   - `annotationType`: Type of annotation (e.g., 'EXON', 'CDS')
#'   - `gene_name`: Name of the associated gene
#'
#' @param genetic_elements A character vector of annotation types to consider when checking for existing elements
#'   (default: c("CDS", "TRANSCRIPT", "MRNA")). This is used to determine whether a given element type
#'   already exists for a transcript; if it does, the CDS will not be created for that transcript.
#'
#' @return A data frame with the original input data and additional rows for the CDS.
#'   Each added row includes the following fields:
#'   - `source`: "JBIO-predicted" for newly added annotations
#'   - `annotationType`: Indicates 'CDS'
#'   - `start` and `end`: Updated start and end positions for the CDS
#'   - `strand`: Strand information copied from the input data
#'   - Other fields as present in the input data
#'
#' @details
#' The function iterates over unique chromosomes and strand orientations, calculating CDS
#' for each gene.
#'
#' @examples
#'
#' # Run the function
#' output_data <- add_CDS(input, genetic_elements = c("CDS", "TRANSCRIPT", "MRNA"))
#'
#' @import stringr dplyr doSNOW foreach parallel doParallel data.table
#' @export
add_CDS <- function(input, genetic_elements = c("TRANSCRIPT", "MRNA", "CDS")) {
  set.seed(123)

  input <- sort_alias(input)

  options(scipen = 999)

  iterations <- length(unique(input$chr))
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)


  CPU <- max(1, detectCores() - 2)

  cl <- makeCluster(CPU)


  registerDoParallel(cl)
  registerDoSNOW(cl)

  chromosomes <- unique(input$chr)

  cat("\n\n CDS adding... \n\n")
  tmp_final <- foreach(chr = chromosomes, .packages = c("dplyr"), .options.snow = opts, .combine = rbind) %dopar% {
    df <- data.frame(matrix(ncol = ncol(input), nrow = 0))
    colnames(df) <- colnames(input)

    for (sen in c("+", "-")) {
      tmp <- input[input$chr %in% chr & input$strand %in% sen, ]

      gen_list <- unique(tmp$gene_name)

      if (length(gen_list) == 0) {
        next
      }


      for (i in 1:length(gen_list)) {
        tmp_transcripts <- tmp[tmp$gene_name %in% gen_list[i], ]

        tmp_transcripts <- tmp_transcripts[tmp_transcripts$transcript_id != "", ]

        tr_vec <- unique(tmp_transcripts$transcript_id)

        if (length(tr_vec) == 0) {
          next
        }

        for (tri in tr_vec) {
          tmp_transcripts_vec <- tmp_transcripts[tmp_transcripts$transcript_id %in% tri, ]

          locus <- NULL
          for (g in genetic_elements) {
            if (g %in% toupper(tmp_transcripts_vec$annotationType)) {
              locus <- g
              break
            }
          }

          if (is.null(locus)) {
            tmp_ex <- tmp_transcripts_vec[toupper(tmp_transcripts_vec$annotationType) %in% "EXON", ]

            tmp_tr <- tmp_ex[1, ]

            tmp_tr$source <- "JBIO-predicted"
            tmp_tr$annotationType <- "CDS"
            tmp_tr$score <- "."
            tmp_tr$phase <- "."
            tmp_tr$start <- as.numeric(min(tmp_ex$start))
            tmp_tr$end <- as.numeric(max(tmp_ex$end))
            tmp_tr$metadata <- ""

            df[nrow(df) + 1, ] <- tmp_tr[1, ]
          }
        }
      }
    }

    return(df)
  }


  close(pb)
  stopCluster(cl)

  tmp_final <- rbind(input, tmp_final)
  tmp_final <- sort_alias(tmp_final)

  return(tmp_final)
}



#' Add Intronic Regions Annotations to Genomic Data
#'
#' @description
#' This function generates intronic regions for each transcript based on its associated exons.
#' It assumes that input data includes properly assigned exons and their transcript IDs. The output consists
#' of intronic regions suitable for further analysis or export in GTF/GFF format.
#'
#' @param input A data frame containing genomic data. The data frame should have the following columns:
#'   - `chr`: Chromosome identifier
#'   - `start`: Start position of the annotation
#'   - `end`: End position of the annotation
#'   - `strand`: Strand information ('+' or '-')
#'   - `annotationType`: Type of annotation (e.g., 'EXON', 'CDS')
#'   - `gene_name`: Name of the associated gene
#'
#'
#' @return A data frame with the original input data and additional rows for the introns.
#'   Each added row includes the following fields:
#'   - `source`: "JBIO-predicted" for newly added annotations
#'   - `annotationType`: Indicates 'intron'
#'   - `start` and `end`: Updated start and end positions for the intron
#'   - `strand`: Strand information copied from the input data
#'   - Other fields as present in the input data
#'
#' @details
#' The function iterates over unique chromosomes and strand orientations, calculating intron positions
#' for each gene.
#'
#' @examples
#'
#' # Run the function
#' output_data <- add_introns(input)
#'
#' @import stringr dplyr doSNOW foreach parallel doParallel data.table
#' @export
add_introns <- function(input) {
  set.seed(123)

  input <- sort_alias(input)

  options(scipen = 999)

  iterations <- length(unique(input$chr))
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)


  CPU <- max(1, detectCores() - 2)

  cl <- makeCluster(CPU)


  registerDoParallel(cl)
  registerDoSNOW(cl)

  chromosomes <- unique(input$chr)

  cat("\n\n introns adding... \n\n")
  tmp_final <- foreach(chr = chromosomes, .packages = c("dplyr"), .options.snow = opts, .combine = rbind) %dopar% {
    df <- data.frame(matrix(ncol = ncol(input), nrow = 0))
    colnames(df) <- colnames(input)

    for (sen in c("+", "-")) {
      tmp <- input[input$chr %in% chr & input$strand %in% sen, ]

      gen_list <- unique(tmp$gene_name)

      if (length(gen_list) == 0) {
        next
      }


      for (i in 1:length(gen_list)) {
        tmp_transcripts <- tmp[tmp$gene_name %in% gen_list[i], ]

        tmp_transcripts <- tmp_transcripts[tmp_transcripts$transcript_id != "", ]

        tr_vec <- unique(tmp_transcripts$transcript_id)

        if (length(tr_vec) == 0) {
          next
        }


        for (tri in tr_vec) {
          tmp_transcripts_vec <- tmp_transcripts[tmp_transcripts$transcript_id %in% tri, ]


          if ("EXON" %in% toupper(tmp_transcripts_vec$annotationType)) {
            tmp_ex <- tmp_transcripts_vec[toupper(tmp_transcripts_vec$annotationType) %in% "EXON", ]

            if (nrow(tmp_ex) > 1) {
              tmp_ex <- tmp_ex[order(tmp_ex$start), ]

              for (exi in 1:nrow(tmp_ex)) {
                if (exi != 1) {
                  tmp_tr <- tmp_ex[1, ]

                  tmp_tr$source <- "JBIO-predicted"
                  tmp_tr$annotationType <- "intron"
                  tmp_tr$score <- "."
                  tmp_tr$phase <- "."
                  tmp_tr$start <- as.numeric(tmp_ex$end[exi - 1]) + 1
                  tmp_tr$end <- as.numeric(tmp_ex$start[exi]) - 1
                  tmp_tr$metadata <- ""

                  df[nrow(df) + 1, ] <- tmp_tr[1, ]
                }
              }
            }
          }
        }
      }
    }

    return(df)
  }


  close(pb)
  stopCluster(cl)

  tmp_final <- rbind(input, tmp_final)
  tmp_final <- sort_alias(tmp_final)

  return(tmp_final)
}




#' Load and Filter Annotation Data (GTF/GFF)
#'
#' @description
#' This function reads an annotation file in GTF (Gene Transfer Format) and GFF, or similar formats. It supports files from
#' multiple sources such as NCBI, Ensembl, and GENCODE. The function parses the file contents and optionally filters the data
#' based on specified genetic elements. It uses various libraries for efficient data manipulation and parallel computing.
#'
#' @param path Character. The file path to the annotation file (GTF, GFF, etc.) to be loaded. The file should be tab-delimited
#' and may include comments prefixed with `#`.
#' @param genetic_elements Character vector (optional). A vector of genetic elements (e.g., "gene", "exon", "CDS") to filter the data.
#' If `NaN` (default), no filtering is performed.
#'
#' @return
#' A data frame (tibble) containing the parsed and optionally filtered annotation data with the following columns:
#' \itemize{
#'   \item `X1` - Character. Chromosome or scaffold name.
#'   \item `X2` - Character. Source or annotation tool used.
#'   \item `X3` - Character. Type of genetic element (e.g., gene, exon).
#'   \item `X4` - Integer. Start position of the feature.
#'   \item `X5` - Integer. End position of the feature.
#'   \item `X6` - Character. Score or confidence level (if available).
#'   \item `X7` - Character. Strand information (`+` or `-`).
#'   \item `X8` - Character. Phase (e.g., 0, 1, 2 for coding sequences).
#'   \item `X9` - Character. Additional attributes in the format of key-value pairs.
#' }
#'
#' @details
#' The function supports annotation files in GTF/GFF formats from widely used sources such as NCBI, Ensembl, and GENCODE.
#' It uses `readr` for efficient file reading and supports filtering based on case-insensitive matching of genetic elements.
#'
#' @examples
#' # Load an annotation file without filtering:
#' annotation_data <- load_annotation("path/to/annotation_file.gtf")
#'
#' # Load an annotation file and filter for genes and exons:
#' annotation_data <- load_annotation("path/to/annotation_file.gtf", genetic_elements = c("gene", "exon"))
#'
#' @import readr stringr dplyr doSNOW foreach parallel doParallel
#' @export
load_annotation <- function(path, genetic_elements = NaN) {
  set.seed(123)

  options(scipen = 999)


  cat("\n\n Data loading... \n\n")


  GTF <- read_delim(path,
    delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE, comment = "#",
    col_types = cols(
      X1 = col_character(),
      X2 = col_character(),
      X3 = col_character(),
      X4 = col_integer(),
      X5 = col_integer(),
      X6 = col_character(),
      X7 = col_character(),
      X8 = col_character(),
      X9 = col_character()
    )
  )

  if (!is.na(genetic_elements)) {
    GTF <- GTF[toupper(GTF$X3) %in% toupper(genetic_elements), ]
  }



  return(GTF)
}


#' Sort Annotation Data by Chromosome and Strand
#'
#' @description
#' This function sorts annotation data by chromosome and strand. It organizes data based on user-defined chromosome
#' ordering and sorts each strand (`+` or `-`) by start and end positions. The function is designed to work with
#' data frames that follow the format of GTF or GFF files.
#' Data must be pre-loaded using load_annotation() or after calling create_GTF_df(), followed by add_UTR(), create_full_GTF(), or create_reduced_GTF().
#'
#' @param input Data frame. The input annotation data. The first column should represent chromosomes, the fourth
#' column should represent start positions, the fifth column end positions, and the seventh column strands (`+` or `-`).
#' @param chromosomes_queue Character vector (optional). A predefined order of chromosomes for sorting.
#' If `NaN` (default), chromosomes will be ordered based on their unique appearance in the input data.
#'
#' @return
#' A data frame containing the sorted annotation data, where:
#' \itemize{
#'   \item Data is grouped and ordered by chromosome as defined in `chromosomes_queue`.
#'   \item Each chromosome group is further sorted by strand (`+` and `-`).
#'   \item Within each strand, data is sorted by start and end positions in ascending order.
#' }
#'
#' @details
#' This function is particularly useful for processing annotation data (e.g., GTF, GFF files) where sorting by
#' genomic position and strand is required. The default chromosome order is determined by the unique appearance
#' of chromosomes in the input data, but users can provide a specific order using `chromosomes_queue`.
#'
#' The function processes each chromosome and strand separately, ensuring data is correctly grouped and sorted.
#'
#' @examples
#' # Example with default chromosome order:
#' sorted_data <- sort_alias(input_data, chromosomes_queue = NaN)
#'
#' # Example with a specific chromosome order:
#' sorted_data <- sort_alias(input_data, chromosomes_queue = c("chr1", "chr2", "chrX", "chrY"))
#'
#' @export
sort_alias <- function(input, chromosomes_queue = NaN) {
  set.seed(123)

  options(scipen = 999)

  chro_col <- colnames(input)[1]
  start_col <- colnames(input)[4]
  end_col <- colnames(input)[5]
  starnd <- colnames(input)[7]
  ti <- colnames(input)[11]



  if (is.na(chromosomes_queue)) {
    chromosomes_queue <- unique(input[[chro_col]])
  }

  output <- data.frame()

  for (ch in chromosomes_queue) {
    for (s in c("+", "-")) {
      cat("\n", paste("CHROMOSOME:", ch, "| STRAND:", s, " sorting...                              "))

      tmp <- input[(input[[chro_col]] %in% ch & input[[starnd]] %in% s), ]

      tmp <- tmp[order(tmp[[start_col]], -tmp[[end_col]]), ]

      tmp$gene_name <- factor(tmp$gene_name, levels = unique(tmp$gene_name))


      tmp <- tmp %>%
        group_by(gene_name) %>%
        arrange(.data[[ti]], .by_group = TRUE) %>%
        ungroup()

      tmp$gene_name <- as.character(tmp$gene_name)

      output <- rbind(output, tmp)
    }
  }

  return(output)
}




#' @title Helper for create_GTF_df
#' @description Internal function used inside create_GTF_df.
#' @keywords internal
#' @import stringr dplyr doSNOW foreach parallel doParallel data.table
optimize_gtf <- function(df, shift = 100000) {
  options(scipen = 999)

  df$lp <- 1:length(df$X1)

  # repair duplicated genes names from different loci

  duplicated_names <- c()


  # find annotation shifted on shift parameter (same chromosome different location)

  for (pm in c("+", "-")) {
    for (ch in unique(df$X1)) {
      tmp <- df[(df$X1 %in% ch & df$X7 %in% pm), ]
      tmp <- tmp[tmp$gene_name != "", ]

      if (nrow(tmp) > 0) {
        tmp <- tmp %>%
          group_by(gene_name) %>%
          summarise(
            start_min = min(X4),
            start_max = max(X4),
            diff = start_max - start_min
          )

        to_add <- as.vector(tmp$gene_name[tmp$diff >= shift])
        if (length(to_add) > 0) {
          duplicated_names <- c(duplicated_names, to_add)
        }
      }
    }
  }


  # find annotation form different strands

  tmp <- df %>%
    group_by(gene_name) %>%
    summarise(
      has_plus = any(X7 == "+"),
      has_minus = any(X7 == "-"),
      both_strands = has_plus & has_minus
    )


  if (TRUE %in% tmp$both_strands) {
    duplicated_names <- c(duplicated_names, as.vector(tmp$gene_name[tmp$both_strands == TRUE]))
  }


  # find annotations from different chromosomes

  tmp <- df %>%
    group_by(gene_name) %>%
    summarise(
      unique_X1 = n_distinct(X1),
      multiple_X1 = unique_X1 > 1
    )


  if (TRUE %in% tmp$multiple_X1) {
    duplicated_names <- c(duplicated_names, as.vector(tmp$gene_name[tmp$multiple_X1 == TRUE]))
  }



  duplicated_names <- unique(duplicated_names)

  duplicated_names <- duplicated_names[duplicated_names != ""]

  if (length(duplicated_names) > 0) {
    # data preparing

    dedup_df <- df[, colnames(df) %in% c("X1", "X7", "X3", "X4", "X5", "gene_name", "transcript_id", "lp")]
    dedup_df <- dedup_df[dedup_df$gene_name %in% duplicated_names, ]
    # dedup_df <- sort_alias(dedup_df)


    # chcek if in on different chromosome the same annotation

    dedup_df <- dedup_df %>%
      group_by(gene_name) %>%
      mutate(check_chr = n_distinct(X1) > 1) %>%
      ungroup()

    # chcek if on both strand the same annotation

    dedup_df <- dedup_df %>%
      group_by(gene_name) %>%
      mutate(check_str = n_distinct(X7) > 1) %>%
      ungroup()

    dedup_df$renamed <- gsub(" ", "", gsub('"', "", dedup_df$gene_name))


    dedup_df$renamed[dedup_df$check_chr == TRUE] <- paste0(dedup_df$renamed[dedup_df$check_chr == TRUE], ".", dedup_df$X1[dedup_df$check_chr == TRUE])


    dedup_df$renamed[dedup_df$check_str == TRUE] <- paste0(dedup_df$renamed[dedup_df$check_str == TRUE], ".str", dedup_df$X7[dedup_df$check_str == TRUE])

    dedup_df <- dedup_df[, !colnames(dedup_df) %in% c("check_str", "check_chr")]


    duplicated_names <- unique(dedup_df$renamed)


    cat("\n\n Duplicated genes repairing...             \n\n ")

    results_list <- list()
    result_index <- 1

    for (gen in seq_along(duplicated_names)) {
      if (length(duplicated_names[gen]) == 0) {
        next
      }

      cat("\r", paste("GENE:", duplicated_names[gen], "| PROGRESS:", round((gen / length(duplicated_names) * 100), 2), "%                           "))

      tmp <- dedup_df[dedup_df$renamed %in% duplicated_names[gen], ]

      group <- 1

      for (i in seq_len(nrow(tmp))) {
        if (nrow(tmp) == 1) {
          tmp$group <- group
          results_list[[result_index]] <- tmp
          result_index <- result_index + 1
        } else if (i == 1) {
          tmp1 <- tmp[i, ]
          tmp1$group <- group
          results_list[[result_index]] <- tmp1
          result_index <- result_index + 1
        } else {
          tmp1 <- tmp[i - 1, ]
          tmp2 <- tmp[i, ]

          if (tmp2$X1[1] == tmp1$X1[1] &&
            tmp2$X4[1] - tmp1$X5[length(tmp1$X5)] > shift &&
            tmp2$transcript_id[1] != tmp1$transcript_id[1]) {
            group <- group + 1
          }

          tmp2$group <- group
          results_list[[result_index]] <- tmp2
          result_index <- result_index + 1
        }
      }
    }


    # global_df <- do.call(rbind, results_list)
    global_df <- rbindlist(results_list, use.names = TRUE, fill = TRUE)

    global_df <- global_df %>%
      group_by(renamed) %>%
      mutate(check = n_distinct(group) > 1) %>%
      ungroup()


    global_df$renamed[global_df$check == TRUE] <- paste0(global_df$renamed[global_df$check == TRUE], ".var", global_df$group[global_df$check == TRUE])

    global_df <- global_df[, !colnames(global_df) %in% c("check", "group")]


    rm(dedup_df)


    # rename genes

    idx <- match(global_df$lp, df$lp)
    valid <- !is.na(idx)
    df$gene_name[idx[valid]] <- global_df$renamed[valid]

    rm(global_df)
  }


  return(df[, !colnames(df) %in% c("lp")])
}


#' @title Convert and Process GTF DataFrame
#'
#' @description
#' The `create_GTF_df` function processes a GTF (Gene Transfer Format) input dataset to extract and optimize
#' gene, transcript, and associated annotations. It handles various GTF formats such as those from GENCODE,
#' Ensembl, NCBI, or custom sources. The function also resolves duplicated gene names, unifies annotations,
#' and provides flexibility for optimization.
#' Data must be pre-loaded using load_annotation() and followed by sort_alias().
#'
#'
#' @param input A data frame containing GTF data. The GTF file must be pre-loaded as a data frame and
#' should have at least 9 columns with annotation data in the 9th column.
#' @param optimize Logical (default: TRUE). If `TRUE`, the function performs optimization steps
#' including filling missing annotations, removing redundant rows, and unifying gene names.
#' @param shift Numeric (default: 100000). Determines the threshold for resolving duplicated gene names
#' based on genomic locus proximity.
#'
#' @details
#' This function is designed to handle different formats of GTF annotations:
#' - **GENCODE/Ensembl**: Extracts `gene_id`, `gene_name`, `transcript_id`, and `transcript_name`.
#' - **NCBI**: Extracts `GeneID`, `gene_name`, and `GenBank` transcript identifiers.
#' - **Custom Format**: Parses annotations containing `gene:` and `transcript:` prefixes.
#'
#' The function includes:
#' - Optimization to merge missing or inconsistent annotations across the dataset.
#' - Detection and repair of duplicated gene names appearing at different loci or strands.
#' - Parallelized processing for improved performance.
#'
#' @return A processed data frame with standardized GTF fields, including:
#' - `gene_name`
#' - `gene_id`
#' - `transcript_id`
#' - `transcript_name`
#'
#' Additional columns depend on the input data frame structure.
#'
#' @examples
#'
#' # Process the GTF data
#' processed_gtf <- create_GTF_df(gtf_data, optimize = TRUE, shift = 50000)
#'
#' @import stringr dplyr doSNOW foreach parallel doParallel data.table
#' @export
create_GTF_df <- function(input, optimize = TRUE, shift = 100000) {
  set.seed(123)

  options(scipen = 999)

  cat("\n\n GTF converting... \n\n")

  df <- input[, 1:8]

  df$bind <- 1:nrow(df)

  input$bind <- 1:nrow(df)


  if (TRUE %in% grepl("gene_name", input$X9) | TRUE %in% grepl("gene_id", input$X9)) {
    # GENECODE & https://www.ensembl.org/index.html
    #############################################################


    # gene id
    df$gene_id_check <- grepl("gene_id", input$X9)

    df$gene_id <- ifelse(
      grepl("gene_id", input$X9),
      gsub(" ", "", gsub('"', "", gsub(".*=", "", gsub(";.*", "", gsub(".*gene_id", "", input$X9))))),
      ""
    )


    df$gene_id_check <- df$gene_id != ""


    # gene name
    df$gene_name_check <- grepl("gene_name", input$X9)

    df$gene_name <- ifelse(
      grepl("gene_name", input$X9),
      gsub(" ", "", gsub('"', "", gsub(".*=", "", gsub(";.*", "", gsub(".*gene_name", "", input$X9))))),
      ""
    )


    df$gene_name_check <- df$gene_name != ""



    # transcript name
    df$transcript_name_check <- grepl("transcript_name", input$X9)

    df$transcript_name <- ifelse(
      grepl("transcript_name", input$X9),
      gsub(" ", "", gsub('"', "", gsub(".*=", "", gsub(".*?\\s", "", gsub(";.*", "", gsub(".*transcript_name", "", input$X9)))))),
      ""
    )


    df$transcript_name_check <- df$transcript_name != ""



    # transcript id
    df$transcript_id_check <- grepl("transcript_id", input$X9)

    df$transcript_id <- ifelse(
      grepl("transcript_id", input$X9),
      gsub(" ", "", gsub('"', "", gsub(".*=", "", gsub(".*?\\s", "", gsub(";.*", "", gsub(".*transcript_id", "", input$X9)))))),
      ""
    )

    df$transcript_id_check <- df$transcript_id != ""

    # gene_type

    type <- grepl("gene_type", input$X9)
    biotype <- grepl("gene_biotype", input$X9)


    if (TRUE %in% type) {
      df$gene_type <- ifelse(
        type,
        gsub(" ", "", gsub('"', "", gsub(".*=", "", gsub(";.*", "", gsub(".*gene_type", "", input$X9))))),
        ""
      )
    } else if (TRUE %in% biotype) {
      df$gene_type <- ifelse(
        biotype,
        gsub(" ", "", gsub('"', "", gsub(".*=", "", gsub(";.*", "", gsub(".*gene_biotype", "", input$X9))))),
        ""
      )
    } else {
      df$gene_type <- ""
    }



    # reapiring lack of gene_names


    df <- df[!(df$gene_name_check == FALSE & df$gene_id_check == FALSE & df$transcript_id_check == FALSE & df$transcript_name_check == FALSE), ]

    df$gene_id[df$gene_id_check == FALSE] <- df$transcript_id[df$gene_id_check == FALSE]


    tmp_gn <- distinct(df[, c("gene_id", "gene_name")])

    df$gene_name <- tmp_gn$gene_name[match(df$gene_id, tmp_gn$gene_id)]

    df$gene_name_check <- df$gene_name != ""

    df$gene_name[df$gene_name_check == FALSE] <- df$gene_id[df$gene_name_check == FALSE]

    df$transcript_name[df$transcript_name_check == FALSE] <- df$gene_name[df$transcript_name_check == FALSE]


    tmp_type <- distinct(df[, c("gene_type", "gene_name")])


    df$gene_type <- tmp_type$gene_type[match(df$gene_name, tmp_type$gene_name)]

    rm(tmp_type)


    ############################################################################


    df <- df[, !colnames(df) %in% c("gene_id_check", "gene_name_check", "transcript_name_check", "transcript_id_check")]


    if (optimize) {
      df <- optimize_gtf(df = df, shift = shift)
    }
  } else if (TRUE %in% grepl("=gene-", input$X9) | TRUE %in% grepl("gene=", input$X9)) {
    # NCBI
    ################################################################


    # gene id
    df$gene_id_check <- grepl("GeneID", input$X9)

    df$gene_id <- ifelse(
      df$gene_id_check,
      gsub(".*GeneID:([0-9]+).*", "\\1", input$X9),
      ""
    )

    df$gene_id_check <- df$gene_id != ""


    # gene name 1
    df$gene_name_check1 <- grepl("=gene-", input$X9)

    df$gene_name1 <- ifelse(
      df$gene_name_check1,
      gsub(" ", "", gsub('"', "", gsub(";.*", "", gsub(".*=gene-", "", input$X9)))),
      ""
    )


    df$gene_name_check1 <- df$gene_name1 != ""


    # gene name 1
    df$gene_name_check2 <- grepl(";gene=", input$X9)

    df$gene_name2 <- ifelse(
      df$gene_name_check2,
      gsub(" ", "", gsub('"', "", gsub(";.*", "", gsub(".*;gene=", "", input$X9)))),
      ""
    )


    df$gene_name_check2 <- df$gene_name2 != ""


    df$gene_name <- ""
    df$gene_name[df$gene_name_check2 == TRUE] <- df$gene_name2[df$gene_name_check2 == TRUE]
    df$gene_name[df$gene_name_check1 == TRUE] <- df$gene_name1[df$gene_name_check1 == TRUE]



    # transcript name

    df$transcript_name <- ""


    # transcript id
    df$transcript_id_check1 <- grepl("Genbank:", input$X9)

    df$transcript_id1 <- ifelse(
      df$transcript_id_check1,
      gsub(" ", "", gsub('"', "", gsub(";.*", "", gsub(",.*", "", gsub(".*Genbank:", "", input$X9))))),
      ""
    )

    df$transcript_id_check1 <- df$transcript_id1 != ""



    # transcript id
    df$transcript_id_check2 <- grepl("transcript_id=", input$X9)

    df$transcript_id2 <- ifelse(
      df$transcript_id_check2,
      gsub(" ", "", gsub('"', "", gsub(";.*", "", gsub(",.*", "", gsub(".*transcript_id=", "", input$X9))))),
      ""
    )

    df$transcript_id_check2 <- df$transcript_id2 != ""



    df$transcript_id <- ""
    df$transcript_id[df$transcript_id_check2 == TRUE] <- df$transcript_id2[df$transcript_id_check2 == TRUE]
    df$transcript_id[df$transcript_id_check1 == TRUE] <- df$transcript_id1[df$transcript_id_check1 == TRUE]




    # gene_type

    type <- grepl("gene_type", input$X9)
    biotype <- grepl("gene_biotype", input$X9)


    if (TRUE %in% type) {
      df$gene_type <- ifelse(
        type,
        gsub(" ", "", gsub('"', "", gsub(".*=", "", gsub(";.*", "", gsub(".*gene_type", "", input$X9))))),
        ""
      )
    } else if (TRUE %in% biotype) {
      df$gene_type <- ifelse(
        biotype,
        gsub(" ", "", gsub('"', "", gsub(".*=", "", gsub(";.*", "", gsub(".*gene_biotype", "", input$X9))))),
        ""
      )
    } else {
      df$gene_type <- ""
    }




    # reapiring lack of gene_names


    df <- df[!(df$gene_name_check1 == FALSE & df$gene_name_check2 == FALSE & df$gene_id_check == FALSE & df$transcript_id_check1 == FALSE & df$transcript_id_check2 == FALSE), ]

    df$gene_id[df$gene_id_check == FALSE] <- df$transcript_id[df$gene_id_check == FALSE]


    tmp_gn <- distinct(df[, c("gene_id", "gene_name")])

    df$gene_name <- tmp_gn$gene_name[match(df$gene_id, tmp_gn$gene_id)]


    # df$gene_name[df$gene_name_check1 == FALSE | df$gene_name_check2 == FALSE] <-  df$gene_id[df$gene_name_check1 == FALSE | df$gene_name_check2 == FALSE]

    df$transcript_name[df$transcript_id_check1 == TRUE | df$transcript_id_check2 == TRUE] <- df$gene_name[df$transcript_id_check1 == TRUE | df$transcript_id_check2 == TRUE]

    tmp_type <- distinct(df[, c("gene_type", "gene_name")])

    df$gene_type <- tmp_type$gene_type[match(df$gene_name, tmp_type$gene_name)]

    rm(tmp_type)


    ############################################################################

    df <- df[, !colnames(df) %in% c("gene_id_check", "gene_name_check1", "gene_name_check2", "transcript_id_check1", "transcript_id_check2", "gene_name1", "gene_name2", "transcript_id1", "transcript_id2")]


    if (optimize) {
      df <- optimize_gtf(df = df, shift = shift)
    }
  }


  df$metadata <- input$X9[match(df$bind, input$bind)]

  df <- df[, !colnames(df) %in% c("bind")]


  colnames(df) <- c("chr", "source", "annotationType", "start", "end", "score", "strand", "phase", "gene_id", "gene_name", "transcript_name", "transcript_id", "gene_type", "metadata")

  return(df)
}

#' Add UTR (Untranslated Region) Annotations to Genomic Data
#'
#' @description
#' This function extends genomic annotations by adding predicted 5' and 3' untranslated regions (UTRs)
#' based on input genomic data. It identifies genes within specified genetic elements and adjusts UTR lengths
#' dynamically based on proximity to neighboring genes.
#' Data must be pre-loaded using load_annotation() followed by create_GTF_df().
#'
#' @param input A data frame containing genomic data. The data frame should have the following columns:
#'   - `chr`: Chromosome identifier
#'   - `start`: Start position of the annotation
#'   - `end`: End position of the annotation
#'   - `strand`: Strand information ('+' or '-')
#'   - `annotationType`: Type of annotation (e.g., 'EXON', 'CDS')
#'   - `gene_name`: Name of the associated gene
#'
#' @param five_prime_utr_length Integer, the default length of the 5' UTR to add (default: 400).
#' @param three_prime_utr_length Integer, the default length of the 3' UTR to add (default: 800).
#' @param genetic_elements A character vector of annotation types to consider for UTR extension
#'   (default: c("EXON", "CDS", "TRANSCRIPT", "MRNA")).
#'
#' @return A data frame with the original input data and additional rows for the predicted UTRs and transcripts.
#'   Each added row includes the following fields:
#'   - `source`: "JBIO-predicted" for newly added annotations
#'   - `annotationType`: Indicates 'five_prime_UTR', 'three_prime_UTR', or 'transcript'
#'   - `start` and `end`: Updated start and end positions for the UTRs or transcript
#'   - `strand`: Strand information copied from the input data
#'   - Other fields as present in the input data
#'
#' @details
#' The function iterates over unique chromosomes and strand orientations, calculating appropriate UTR lengths
#' for each gene. It dynamically adjusts the UTR lengths based on available space between genes to avoid overlap.
#'
#' @examples
#'
#' # Run the function
#' output_data <- add_UTR(input, five_prime_utr_length = 400, three_prime_utr_length = 800, genetic_elements = c("EXON", "CDS", "TRANSCRIPT", "MRNA"))
#'
#' @import stringr dplyr doSNOW foreach parallel doParallel data.table
#' @export
add_UTR <- function(input, five_prime_utr_length = 400, three_prime_utr_length = 800, biotype = "protein_coding", transcript_limit = NULL, meta_string = NULL, genetic_elements = c("TRANSCRIPT", "MRNA", "CDS")) {
  genetic_elements <- toupper(genetic_elements)
  set.seed(123)

  options(scipen = 999)

  cat("\n\n UTRs sequence extending...             \n\n")

  tmp_final <- list()
  il <- 0
  final <- input

  # input$sort_val <- 1:length(input$chr)
  chromosomes <- unique(input$chr)

  CDS <- final[toupper(input$annotationType) %in% toupper(genetic_elements), ]


  to_extend <- CDS[CDS$transcript_id != "", ]


  if (!is.null(meta_string)) {
    to_extend <- to_extend[grepl(pattern = toupper(meta_string), toupper(to_extend$metadata)), ]
  }

  if (!is.null(transcript_limit)) {
    to_extend$diff <- to_extend$end - to_extend$start
    to_extend <- to_extend[to_extend$diff > transcript_limit, ]
  }


  if (!is.null(biotype)) {
    if (biotype %in% CDS$gene_type) {
      to_extend <- to_extend[to_extend$gene_type %in% biotype, ]
    } else {
      cat(paste0("Provided biotype: ", as.character(biotype), " does not exist in data!"))
    }
  }

  to_extend_genes <- unique(to_extend$gene_name)

  rm(to_extend)

  for (chr in chromosomes) {
    for (sen in c("+", "-")) {
      tmp_genes <- CDS[CDS$chr %in% chr & CDS$strand %in% sen, ]

      gen_list <- unique(tmp_genes$gene_name)

      if (length(gen_list) > 0) {
        for (i in 1:length(gen_list)) {
          if (!gen_list[i] %in% to_extend_genes) {
            next
          }

          cat("\r", paste("LOC:", chr, "| STRAND:", sen, "| PROGRESS:", round((i / length(gen_list) * 100), 2), "%           "))


          tmp_transcripts <- tmp_genes[tmp_genes$gene_name %in% gen_list[i], ]

          locus <- NULL
          for (g in genetic_elements) {
            if (g %in% toupper(tmp_transcripts$annotationType)) {
              locus <- g
              break
            }
          }

          if (is.null(locus)) {
            next
          }

          tmp_transcripts <- tmp_transcripts[toupper(tmp_transcripts$annotationType) %in% locus, ]


          # debug

          # if (i == 12) {break}

          # tmp variables

          for (tix in rownames(tmp_transcripts)) {
            tmp <- tmp_transcripts[tix, ]


            five_prime_utr_length_cor_minus <- NaN
            five_prime_utr_length_cor_plus <- NaN
            three_prime_utr_length_cor_minus <- NaN
            three_prime_utr_length_cor_plus <- NaN


            if (i == 1) {
              # first element


              tmp2 <- tmp_genes[tmp_genes$gene_name %in% gen_list[i + 1], ]


              if (min(tmp2$start) > (max(tmp$end) + ((five_prime_utr_length + three_prime_utr_length) + 2))) {
                if (sen == "+") {
                  five_prime_utr_length_cor_plus <- round(five_prime_utr_length, 0)
                  three_prime_utr_length_cor_plus <- round(three_prime_utr_length, 0)
                } else {
                  five_prime_utr_length_cor_minus <- round(five_prime_utr_length, 0)
                  three_prime_utr_length_cor_minus <- round(three_prime_utr_length, 0)
                }
              } else if (min(tmp2$start) > max(tmp$end) + ((five_prime_utr_length + three_prime_utr_length) / 2)) {
                if (sen == "+") {
                  five_prime_utr_length_cor_plus <- round(five_prime_utr_length, 0)
                  three_prime_utr_length_cor_plus <- round(three_prime_utr_length * 0.5, 0)
                } else {
                  five_prime_utr_length_cor_minus <- round(five_prime_utr_length * 0.5, 0)
                  three_prime_utr_length_cor_minus <- round(three_prime_utr_length, 0)
                }
              } else if (min(tmp2$start) > max(tmp$end) + ((five_prime_utr_length + three_prime_utr_length) / 3)) {
                if (sen == "+") {
                  five_prime_utr_length_cor_plus <- round(five_prime_utr_length, 0)
                  three_prime_utr_length_cor_plus <- round(three_prime_utr_length * 0.3, 0)
                } else {
                  five_prime_utr_length_cor_minus <- round(five_prime_utr_length * 0.3, 0)
                  three_prime_utr_length_cor_minus <- round(three_prime_utr_length, 0)
                }
              } else if (min(tmp2$start) > max(tmp$end)) {
                if (sen == "+") {
                  five_prime_utr_length_cor_plus <- five_prime_utr_length
                  three_prime_utr_length_cor_plus <- round((min(tmp2$start) - max(tmp$end) - 2) / 2, 0)
                  if (three_prime_utr_length_cor_plus < 0) {
                    three_prime_utr_length_cor_plus <- 0
                  }
                } else {
                  five_prime_utr_length_cor_minus <- round((min(tmp2$start) - max(tmp$end) - 2) / 2, 0)
                  if (five_prime_utr_length_cor_minus < 0) {
                    five_prime_utr_length_cor_minus <- 0
                  }
                  three_prime_utr_length_cor_minus <- three_prime_utr_length
                }
              } else {
                if (sen == "+") {
                  five_prime_utr_length_cor_plus <- five_prime_utr_length
                } else {
                  three_prime_utr_length_cor_minus <- three_prime_utr_length
                }
              }




              tmp0_UTR5 <- tmp
              tmp0_UTR5$source <- "JBIO-predicted"
              tmp0_UTR5$annotationType <- "five_prime_UTR"
              tmp0_UTR5$score <- "."
              tmp0_UTR5$phase <- "."


              tmp0_UTR3 <- tmp
              tmp0_UTR3$source <- "JBIO-predicted"
              tmp0_UTR3$annotationType <- "three_prime_UTR"
              tmp0_UTR3$score <- "."
              tmp0_UTR3$phase <- "."

              tmp0_transcript <- tmp
              tmp0_transcript$source <- "JBIO-predicted"
              tmp0_transcript$annotationType <- "transcript"
              tmp0_transcript$score <- "."
              tmp0_transcript$phase <- "."

              if (sen == "+") {
                # UTR5
                if (!is.na(five_prime_utr_length_cor_plus)) {
                  utr5_start <- as.numeric(min(tmp$start) - five_prime_utr_length_cor_plus)
                  if (utr5_start < 0) {
                    utr5_start <- 0
                  }

                  tmp0_UTR5$start <- utr5_start
                  tmp0_UTR5$end <- as.numeric(min(tmp$start) - 1)

                  # tmp_final <- rbind(tmp_final, tmp0_UTR5)
                  il <- il + 1
                  tmp_final[[il]] <- tmp0_UTR5
                }


                # UTR3
                if (!is.na(three_prime_utr_length_cor_plus)) {
                  tmp0_UTR3$start <- as.numeric(max(tmp$end) + 1)
                  tmp0_UTR3$end <- as.numeric(max(tmp$end) + three_prime_utr_length_cor_plus)

                  # tmp_final <- rbind(tmp_final, tmp0_UTR3)
                  il <- il + 1
                  tmp_final[[il]] <- tmp0_UTR3
                }


                # TRANSCRIPT

                if (!is.na(three_prime_utr_length_cor_plus) & !is.na(five_prime_utr_length_cor_plus)) {
                  tmp0_transcript$start <- as.numeric(tmp0_UTR5$start)
                  tmp0_transcript$end <- as.numeric(tmp0_UTR3$end)

                  # tmp_final <- rbind(tmp_final, tmp0_transcript)
                  il <- il + 1
                  tmp_final[[il]] <- tmp0_transcript
                } else if (is.na(three_prime_utr_length_cor_plus) & !is.na(five_prime_utr_length_cor_plus)) {
                  tmp0_transcript$start <- as.numeric(tmp0_UTR5$start)
                  tmp0_transcript$end <- as.numeric(max(tmp$end))

                  # tmp_final <- rbind(tmp_final, tmp0_transcript)
                  il <- il + 1
                  tmp_final[[il]] <- tmp0_transcript
                } else if (!is.na(three_prime_utr_length_cor_plus) & is.na(five_prime_utr_length_cor_plus)) {
                  tmp0_transcript$start <- as.numeric(min(tmp$start))
                  tmp0_transcript$end <- as.numeric(tmp0_UTR3$end)

                  # tmp_final <- rbind(tmp_final, tmp0_transcript)
                  il <- il + 1
                  tmp_final[[il]] <- tmp0_transcript
                } else {
                  tmp0_transcript$start <- as.numeric(min(tmp$start))
                  tmp0_transcript$end <- as.numeric(max(tmp$end))

                  # tmp_final <- rbind(tmp_final, tmp0_transcript)
                  il <- il + 1
                  tmp_final[[il]] <- tmp0_transcript
                }
              } else if (sen == "-") {
                # UTR3
                if (!is.na(three_prime_utr_length_cor_minus)) {
                  utr3_start <- as.numeric(min(tmp$start) - three_prime_utr_length_cor_minus)
                  if (utr3_start < 0) {
                    utr3_start <- 0
                  }

                  tmp0_UTR3$start <- utr3_start
                  tmp0_UTR3$end <- as.numeric(min(tmp$start) - 1)

                  # tmp_final <- rbind(tmp_final, tmp0_UTR3)
                  il <- il + 1
                  tmp_final[[il]] <- tmp0_UTR3
                }


                # UTR5
                if (!is.na(five_prime_utr_length_cor_minus)) {
                  tmp0_UTR5$start <- as.numeric(max(tmp$end) + 1)
                  tmp0_UTR5$end <- as.numeric(max(tmp$end) + five_prime_utr_length_cor_minus)

                  # tmp_final <- rbind(tmp_final, tmp0_UTR5)
                  il <- il + 1
                  tmp_final[[il]] <- tmp0_UTR5
                }

                # TRANSCRIPT
                if (!is.na(five_prime_utr_length_cor_minus) & !is.na(three_prime_utr_length_cor_minus)) {
                  tmp0_transcript$start <- as.numeric(tmp0_UTR3$start)
                  tmp0_transcript$end <- as.numeric(tmp0_UTR5$end)

                  # tmp_final <- rbind(tmp_final, tmp0_transcript)
                  il <- il + 1
                  tmp_final[[il]] <- tmp0_transcript
                } else if (is.na(five_prime_utr_length_cor_minus) & !is.na(three_prime_utr_length_cor_minus)) {
                  tmp0_transcript$start <- as.numeric(tmp0_UTR3$start)
                  tmp0_transcript$end <- as.numeric(max(tmp$end))

                  # tmp_final <- rbind(tmp_final, tmp0_transcript)
                  il <- il + 1
                  tmp_final[[il]] <- tmp0_transcript
                } else if (!is.na(five_prime_utr_length_cor_minus) & is.na(three_prime_utr_length_cor_minus)) {
                  tmp0_transcript$start <- as.numeric(min(tmp$start))
                  tmp0_transcript$end <- as.numeric(tmp0_UTR5$end)

                  # tmp_final <- rbind(tmp_final, tmp0_transcript)
                  il <- il + 1
                  tmp_final[[il]] <- tmp0_transcript
                } else {
                  tmp0_transcript$start <- as.numeric(min(tmp$start))
                  tmp0_transcript$end <- as.numeric(max(tmp$end))

                  # tmp_final <- rbind(tmp_final, tmp0_transcript)
                  il <- il + 1
                  tmp_final[[il]] <- tmp0_transcript
                }
              }
            } else if (i < length(gen_list)) {
              # middle element

              tmp_p <- tmp_genes[tmp_genes$gene_name %in% gen_list[i - 1], ]

              tmp2 <- tmp_genes[tmp_genes$gene_name %in% gen_list[i + 1], ]


              if (min(tmp2$start) > (max(tmp$end) + (five_prime_utr_length + three_prime_utr_length) + 2)) {
                three_prime_utr_length_cor_plus <- three_prime_utr_length

                five_prime_utr_length_cor_minus <- five_prime_utr_length

                ################################################################################################
                if (max(tmp_p$end) < (min(tmp$start) - ((three_prime_utr_length + five_prime_utr_length) + 2)) & sen == "+") {
                  five_prime_utr_length_cor_plus <- round(five_prime_utr_length, 0)
                } else if (max(tmp_p$end) < (min(tmp$start) - (five_prime_utr_length + three_prime_utr_length) + 2) & sen == "-") {
                  three_prime_utr_length_cor_minus <- round(three_prime_utr_length, 0)
                  #
                } else if (max(tmp_p$end) < (min(tmp$start) - ((three_prime_utr_length + five_prime_utr_length) / 2)) & sen == "+") {
                  five_prime_utr_length_cor_plus <- round(five_prime_utr_length * 0.5, 0)
                } else if (max(tmp_p$end) < (min(tmp$start) - ((five_prime_utr_length + three_prime_utr_length) / 2)) & sen == "-") {
                  three_prime_utr_length_cor_minus <- round(three_prime_utr_length * 0.5, 0)
                  #
                } else if (max(tmp_p$end) < (min(tmp$start) - ((three_prime_utr_length + five_prime_utr_length) / 3)) & sen == "+") {
                  five_prime_utr_length_cor_plus <- round(five_prime_utr_length * 0.3, 0)
                } else if (max(tmp_p$end) < (min(tmp$start) - ((five_prime_utr_length + three_prime_utr_length) / 3)) & sen == "-") {
                  three_prime_utr_length_cor_minus <- round(three_prime_utr_length * 0.3, 0)
                } else if (max(tmp_p$end) < min(tmp$start)) {
                  if (sen == "+") {
                    five_prime_utr_length_cor_plus <- round((min(tmp$start) - max(tmp_p$end) - 2) / 2, 0)
                    if (five_prime_utr_length_cor_plus < 0) {
                      five_prime_utr_length_cor_plus <- 0
                    }
                  } else {
                    three_prime_utr_length_cor_minus <- round((min(tmp$start) - max(tmp_p$end) - 2) / 2, 0)
                    if (three_prime_utr_length_cor_minus < 0) {
                      three_prime_utr_length_cor_minus <- 0
                    }
                  }
                }
                #################################################################################################
              } else if (min(tmp2$start) > max(tmp$end) + ((five_prime_utr_length + three_prime_utr_length) / 2)) {
                three_prime_utr_length_cor_plus <- round(three_prime_utr_length * 0.5, 0)

                five_prime_utr_length_cor_minus <- round(five_prime_utr_length * 0.5, 0)

                ################################################################################################
                if (max(tmp_p$end) > (min(tmp$start) - ((three_prime_utr_length + five_prime_utr_length) + 2)) & sen == "+") {
                  five_prime_utr_length_cor_plus <- round(five_prime_utr_length, 0)
                } else if (max(tmp_p$end) > (min(tmp$start) - (five_prime_utr_length + three_prime_utr_length) + 2) & sen == "-") {
                  three_prime_utr_length_cor_minus <- round(three_prime_utr_length, 0)
                  #
                } else if (max(tmp_p$end) < (min(tmp$start) - ((three_prime_utr_length + five_prime_utr_length) / 2)) & sen == "+") {
                  five_prime_utr_length_cor_plus <- round(five_prime_utr_length * 0.5, 0)
                } else if (max(tmp_p$end) < (min(tmp$start) - ((five_prime_utr_length + three_prime_utr_length) / 2)) & sen == "-") {
                  three_prime_utr_length_cor_minus <- round(three_prime_utr_length * 0.5, 0)
                  #
                } else if (max(tmp_p$end) < (min(tmp$start) - ((three_prime_utr_length + five_prime_utr_length) / 3)) & sen == "+") {
                  five_prime_utr_length_cor_plus <- round(five_prime_utr_length * 0.3, 0)
                } else if (max(tmp_p$end) < (min(tmp$start) - ((five_prime_utr_length + three_prime_utr_length) / 3)) & sen == "-") {
                  three_prime_utr_length_cor_minus <- round(three_prime_utr_length * 0.3, 0)
                } else if (max(tmp_p$end) < min(tmp$start)) {
                  if (sen == "+") {
                    five_prime_utr_length_cor_plus <- round((min(tmp$start) - max(tmp_p$end) - 2) / 2, 0)
                    if (five_prime_utr_length_cor_plus < 0) {
                      five_prime_utr_length_cor_plus <- 0
                    }
                  } else {
                    three_prime_utr_length_cor_minus <- round((min(tmp$start) - max(tmp_p$end) - 2) / 2, 0)
                    if (three_prime_utr_length_cor_minus < 0) {
                      three_prime_utr_length_cor_minus <- 0
                    }
                  }
                }
                #################################################################################################
              } else if (min(tmp2$start) > max(tmp$end) + ((five_prime_utr_length + three_prime_utr_length) / 3)) {
                three_prime_utr_length_cor_plus <- round(three_prime_utr_length * 0.3, 0)

                five_prime_utr_length_cor_minus <- round(five_prime_utr_length * 0.3, 0)

                ################################################################################################
                if (max(tmp_p$end) > (min(tmp$start) - ((three_prime_utr_length + five_prime_utr_length) + 2)) & sen == "+") {
                  five_prime_utr_length_cor_plus <- round(five_prime_utr_length, 0)
                } else if (max(tmp_p$end) < (min(tmp$start) - ((five_prime_utr_length + three_prime_utr_length) + 2)) & sen == "-") {
                  three_prime_utr_length_cor_minus <- round(three_prime_utr_length, 0)
                  #
                } else if (max(tmp_p$end) < (min(tmp$start) - ((three_prime_utr_length + five_prime_utr_length) / 2)) & sen == "+") {
                  five_prime_utr_length_cor_plus <- round(five_prime_utr_length * 0.5, 0)
                } else if (max(tmp_p$end) < (min(tmp$start) + ((five_prime_utr_length + three_prime_utr_length) / 2)) & sen == "-") {
                  three_prime_utr_length_cor_minus <- round(three_prime_utr_length * 0.5, 0)
                  #
                } else if (max(tmp_p$end) < (min(tmp$start) - ((three_prime_utr_length + five_prime_utr_length) / 3)) & sen == "+") {
                  five_prime_utr_length_cor_plus <- round(five_prime_utr_length * 0.3, 0)
                } else if (max(tmp_p$end) < (min(tmp$start) - ((five_prime_utr_length + three_prime_utr_length) / 3)) & sen == "-") {
                  three_prime_utr_length_cor_minus <- round(three_prime_utr_length * 0.3, 0)
                } else if (max(tmp_p$end) < min(tmp$start)) {
                  if (sen == "+") {
                    five_prime_utr_length_cor_plus <- round((min(tmp$start) - max(tmp_p$end) - 2) / 2, 0)
                    if (five_prime_utr_length_cor_plus < 0) {
                      five_prime_utr_length_cor_plus <- 0
                    }
                  } else {
                    three_prime_utr_length_cor_minus <- round((min(tmp$start) - max(tmp_p$end) - 2) / 2, 0)
                    if (three_prime_utr_length_cor_minus < 0) {
                      three_prime_utr_length_cor_minus <- 0
                    }
                  }
                }
                #################################################################################################
              } else if (min(tmp2$start) > max(tmp$end)) {
                three_prime_utr_length_cor_plus <- round((min(tmp2$start) - max(tmp$end) - 2) / 2, 0)
                if (three_prime_utr_length_cor_plus < 0) {
                  three_prime_utr_length_cor_plus <- 0
                }

                five_prime_utr_length_cor_minus <- round((min(tmp2$start) - max(tmp$end) - 2) / 2, 0)
                if (five_prime_utr_length_cor_minus < 0) {
                  five_prime_utr_length_cor_minus <- 0
                }

                ################################################################################################
                if (max(tmp_p$end) < (min(tmp$start) - ((three_prime_utr_length + five_prime_utr_length) + 2)) & sen == "+") {
                  five_prime_utr_length_cor_plus <- round(five_prime_utr_length, 0)
                } else if (max(tmp_p$end) < (min(tmp$start) - ((five_prime_utr_length + three_prime_utr_length) + 2)) & sen == "-") {
                  three_prime_utr_length_cor_minus <- round(three_prime_utr_length, 0)
                  #
                } else if (max(tmp_p$end) < (min(tmp$start) - ((three_prime_utr_length + five_prime_utr_length) / 2)) & sen == "+") {
                  five_prime_utr_length_cor_plus <- round(five_prime_utr_length * 0.5, 0)
                } else if (max(tmp_p$end) < (min(tmp$start) - ((five_prime_utr_length + three_prime_utr_length) / 2)) & sen == "-") {
                  three_prime_utr_length_cor_minus <- round(three_prime_utr_length * 0.5, 0)
                  #
                } else if (max(tmp_p$end) < (min(tmp$start) - ((three_prime_utr_length + five_prime_utr_length) / 3)) & sen == "+") {
                  five_prime_utr_length_cor_plus <- round(five_prime_utr_length * 0.3, 0)
                } else if (max(tmp_p$end) < (min(tmp$start) - ((five_prime_utr_length + three_prime_utr_length) / 3)) & sen == "-") {
                  three_prime_utr_length_cor_minus <- round(three_prime_utr_length * 0.3, 0)
                } else if (max(tmp_p$end) < min(tmp$start)) {
                  if (sen == "+") {
                    five_prime_utr_length_cor_plus <- round((min(tmp$start) - max(tmp_p$end) - 2) / 2, 0)
                    if (five_prime_utr_length_cor_plus < 0) {
                      five_prime_utr_length_cor_plus <- 0
                    }
                  } else {
                    three_prime_utr_length_cor_minus <- round((min(tmp$start) - max(tmp_p$end) - 2) / 2, 0)
                    if (three_prime_utr_length_cor_minus < 0) {
                      three_prime_utr_length_cor_minus <- 0
                    }
                  }
                }
                #################################################################################################
              } else if (max(tmp_p$end) < min(tmp$start)) {
                if (max(tmp_p$end) < (min(tmp$start) - ((three_prime_utr_length + five_prime_utr_length) + 2)) & sen == "+") {
                  five_prime_utr_length_cor_plus <- round(five_prime_utr_length, 0)
                } else if (max(tmp_p$end) < (min(tmp$start) - ((five_prime_utr_length + three_prime_utr_length) + 2)) & sen == "-") {
                  three_prime_utr_length_cor_minus <- round(three_prime_utr_length, 0)
                  #
                } else if (max(tmp_p$end) < (min(tmp$start) - ((three_prime_utr_length + five_prime_utr_length) / 2)) & sen == "+") {
                  five_prime_utr_length_cor_plus <- round(five_prime_utr_length * 0.5, 0)
                } else if (max(tmp_p$end) < (min(tmp$start) - ((five_prime_utr_length + three_prime_utr_length) / 2)) & sen == "-") {
                  three_prime_utr_length_cor_minus <- round(three_prime_utr_length * 0.5, 0)
                  #
                } else if (max(tmp_p$end) < (min(tmp$start) - ((three_prime_utr_length + five_prime_utr_length) / 3)) & sen == "+") {
                  five_prime_utr_length_cor_plus <- round(five_prime_utr_length * 0.3, 0)
                } else if (max(tmp_p$end) < (min(tmp$start) - ((five_prime_utr_length + three_prime_utr_length) / 3)) & sen == "-") {
                  three_prime_utr_length_cor_minus <- round(three_prime_utr_length * 0.3, 0)
                } else if (max(tmp_p$end) < min(tmp$start)) {
                  if (sen == "+") {
                    five_prime_utr_length_cor_plus <- round((min(tmp$start) - max(tmp_p$end) - 2) / 2, 0)
                    if (five_prime_utr_length_cor_plus < 0) {
                      five_prime_utr_length_cor_plus <- 0
                    }
                  } else {
                    three_prime_utr_length_cor_minus <- round((min(tmp$start) - max(tmp_p$end) - 2) / 2, 0)
                    if (three_prime_utr_length_cor_minus < 0) {
                      three_prime_utr_length_cor_minus <- 0
                    }
                  }
                }
              }



              tmp0_UTR5 <- tmp
              tmp0_UTR5$source <- "JBIO-predicted"
              tmp0_UTR5$annotationType <- "five_prime_UTR"
              tmp0_UTR5$score <- "."
              tmp0_UTR5$phase <- "."


              tmp0_UTR3 <- tmp
              tmp0_UTR3$source <- "JBIO-predicted"
              tmp0_UTR3$annotationType <- "three_prime_UTR"
              tmp0_UTR3$score <- "."
              tmp0_UTR3$phase <- "."

              tmp0_transcript <- tmp
              tmp0_transcript$source <- "JBIO-predicted"
              tmp0_transcript$annotationType <- "transcript"
              tmp0_transcript$score <- "."
              tmp0_transcript$phase <- "."

              if (sen == "+") {
                # UTR5
                if (!is.na(five_prime_utr_length_cor_plus)) {
                  tmp0_UTR5$start <- as.numeric(min(tmp$start) - five_prime_utr_length_cor_plus)
                  tmp0_UTR5$end <- as.numeric(min(tmp$start) - 1)

                  # tmp_final <- rbind(tmp_final, tmp0_UTR5)
                  il <- il + 1
                  tmp_final[[il]] <- tmp0_UTR5
                }


                # UTR3
                if (!is.na(three_prime_utr_length_cor_plus)) {
                  tmp0_UTR3$start <- as.numeric(max(tmp$end) + 1)
                  tmp0_UTR3$end <- as.numeric(max(tmp$end) + three_prime_utr_length_cor_plus)

                  # tmp_final <- rbind(tmp_final, tmp0_UTR3)
                  il <- il + 1
                  tmp_final[[il]] <- tmp0_UTR3
                }


                # TRANSCRIPT
                if (!is.na(three_prime_utr_length_cor_plus) & !is.na(five_prime_utr_length_cor_plus)) {
                  tmp0_transcript$start <- as.numeric(tmp0_UTR5$start)
                  tmp0_transcript$end <- as.numeric(tmp0_UTR3$end)

                  # tmp_final <- rbind(tmp_final, tmp0_transcript)
                  il <- il + 1
                  tmp_final[[il]] <- tmp0_transcript
                } else if (is.na(three_prime_utr_length_cor_plus) & !is.na(five_prime_utr_length_cor_plus)) {
                  tmp0_transcript$start <- as.numeric(tmp0_UTR5$start)
                  tmp0_transcript$end <- as.numeric(max(tmp$end))

                  # tmp_final <- rbind(tmp_final, tmp0_transcript)
                  il <- il + 1
                  tmp_final[[il]] <- tmp0_transcript
                } else if (!is.na(three_prime_utr_length_cor_plus) & is.na(five_prime_utr_length_cor_plus)) {
                  tmp0_transcript$start <- as.numeric(min(tmp$start))
                  tmp0_transcript$end <- as.numeric(tmp0_UTR3$end)

                  # tmp_final <- rbind(tmp_final, tmp0_transcript)
                  il <- il + 1
                  tmp_final[[il]] <- tmp0_transcript
                } else {
                  tmp0_transcript$start <- as.numeric(min(tmp$start))
                  tmp0_transcript$end <- as.numeric(max(tmp$end))

                  # tmp_final <- rbind(tmp_final, tmp0_transcript)
                  il <- il + 1
                  tmp_final[[il]] <- tmp0_transcript
                }
              } else if (sen == "-") {
                # UTR3
                if (!is.na(three_prime_utr_length_cor_minus)) {
                  tmp0_UTR3$start <- as.numeric(min(tmp$start) - three_prime_utr_length_cor_minus)
                  tmp0_UTR3$end <- as.numeric(min(tmp$start) - 1)

                  # tmp_final <- rbind(tmp_final, tmp0_UTR3)
                  il <- il + 1
                  tmp_final[[il]] <- tmp0_UTR3
                }


                # UTR5
                if (!is.na(five_prime_utr_length_cor_minus)) {
                  tmp0_UTR5$start <- as.numeric(max(tmp$end) + 1)
                  tmp0_UTR5$end <- as.numeric(max(tmp$end) + five_prime_utr_length_cor_minus)

                  # tmp_final <- rbind(tmp_final, tmp0_UTR5)
                  il <- il + 1
                  tmp_final[[il]] <- tmp0_UTR5
                }


                # TRANSCRIPT
                if (!is.na(five_prime_utr_length_cor_minus) & !is.na(three_prime_utr_length_cor_minus)) {
                  tmp0_transcript$start <- as.numeric(tmp0_UTR3$start)
                  tmp0_transcript$end <- as.numeric(tmp0_UTR5$end)

                  # tmp_final <- rbind(tmp_final, tmp0_transcript)
                  il <- il + 1
                  tmp_final[[il]] <- tmp0_transcript
                } else if (is.na(five_prime_utr_length_cor_minus) & !is.na(three_prime_utr_length_cor_minus)) {
                  tmp0_transcript$start <- as.numeric(tmp0_UTR3$start)
                  tmp0_transcript$end <- as.numeric(max(tmp$end))

                  # tmp_final <- rbind(tmp_final, tmp0_transcript)
                  il <- il + 1
                  tmp_final[[il]] <- tmp0_transcript
                } else if (!is.na(five_prime_utr_length_cor_minus) & is.na(three_prime_utr_length_cor_minus)) {
                  tmp0_transcript$start <- as.numeric(min(tmp$start))
                  tmp0_transcript$end <- as.numeric(tmp0_UTR5$end)

                  # tmp_final <- rbind(tmp_final, tmp0_transcript)
                  il <- il + 1
                  tmp_final[[il]] <- tmp0_transcript
                } else {
                  tmp0_transcript$start <- as.numeric(min(tmp$start))
                  tmp0_transcript$end <- as.numeric(max(tmp$end))

                  # tmp_final <- rbind(tmp_final, tmp0_transcript)
                  il <- il + 1
                  tmp_final[[il]] <- tmp0_transcript
                }
              }
            } else if (i == length(gen_list)) {
              # final element


              tmp2 <- tmp_genes[tmp_genes$gene_name %in% gen_list[i - 1], ]


              # End UTR

              if (max(tmp2$end) < max(tmp$end)) {
                three_prime_utr_length_cor_plus <- round(three_prime_utr_length, 0)
                five_prime_utr_length_cor_minus <- round(five_prime_utr_length, 0)
              }


              # tmp2 -> tmp1 +
              if (max(tmp2$end) < (min(tmp$start) - ((three_prime_utr_length + five_prime_utr_length) + 2)) & sen == "+") {
                five_prime_utr_length_cor_plus <- round(five_prime_utr_length, 0)
              } else if (max(tmp2$end) < (min(tmp$start) - ((three_prime_utr_length + five_prime_utr_length) / 2)) & sen == "+") {
                five_prime_utr_length_cor_plus <- round(five_prime_utr_length * 0.5, 0)
              } else if (max(tmp2$end) < (min(tmp$start) - ((three_prime_utr_length + five_prime_utr_length) / 3)) & sen == "+") {
                five_prime_utr_length_cor_plus <- round(five_prime_utr_length * 0.3, 0)
              } else if (max(tmp2$end) < min(tmp$start) & sen == "+") {
                five_prime_utr_length_cor_plus <- round((min(tmp$start) - max(tmp2$end) - 2) / 2, 0)
                if (five_prime_utr_length_cor_plus < 0) {
                  five_prime_utr_length_cor_plus <- 0
                }
              }


              # tmp2 -> tmp1 -
              if (max(tmp2$end) < (min(tmp$start) - ((five_prime_utr_length + three_prime_utr_length) + 2)) & sen == "-") {
                three_prime_utr_length_cor_minus <- round(three_prime_utr_length, 0)
              } else if (max(tmp2$end) < (min(tmp$start) - ((five_prime_utr_length + three_prime_utr_length) / 2)) & sen == "-") {
                three_prime_utr_length_cor_minus <- round(three_prime_utr_length * 0.5, 0)
              } else if (max(tmp2$end) < (min(tmp$start) - ((five_prime_utr_length + three_prime_utr_length) / 3)) & sen == "-") {
                three_prime_utr_length_cor_minus <- round(three_prime_utr_length * 0.3, 0)
              } else if (max(tmp2$end) < min(tmp$start) & sen == "-") {
                three_prime_utr_length_cor_minus <- round((min(tmp$start) - max(tmp2$end) - 2) / 2, 0)
                if (three_prime_utr_length_cor_minus < 0) {
                  three_prime_utr_length_cor_minus <- 0
                }
              }



              tmp0_UTR5 <- tmp
              tmp0_UTR5$source <- "JBIO-predicted"
              tmp0_UTR5$annotationType <- "five_prime_UTR"
              tmp0_UTR5$score <- "."
              tmp0_UTR5$phase <- "."


              tmp0_UTR3 <- tmp
              tmp0_UTR3$source <- "JBIO-predicted"
              tmp0_UTR3$annotationType <- "three_prime_UTR"
              tmp0_UTR3$score <- "."
              tmp0_UTR3$phase <- "."

              tmp0_transcript <- tmp
              tmp0_transcript$source <- "JBIO-predicted"
              tmp0_transcript$annotationType <- "transcript"
              tmp0_transcript$score <- "."
              tmp0_transcript$phase <- "."

              if (sen == "+") {
                # UTR5
                if (!is.na(five_prime_utr_length_cor_plus)) {
                  utr5_start <- as.numeric(min(tmp$start) - five_prime_utr_length_cor_plus)
                  if (utr5_start < 0) {
                    utr5_start <- 0
                  }

                  tmp0_UTR5$start <- utr5_start
                  tmp0_UTR5$end <- as.numeric(min(tmp$start) - 1)

                  # tmp_final <- rbind(tmp_final, tmp0_UTR5)
                  il <- il + 1
                  tmp_final[[il]] <- tmp0_UTR5
                }



                # UTR3
                if (!is.na(three_prime_utr_length_cor_plus)) {
                  tmp0_UTR3$start <- as.numeric(max(tmp$end) + 1)
                  tmp0_UTR3$end <- as.numeric(max(tmp$end) + three_prime_utr_length_cor_plus)

                  # tmp_final <- rbind(tmp_final, tmp0_UTR3)
                  il <- il + 1
                  tmp_final[[il]] <- tmp0_UTR3
                }


                # TRANSCRIPT
                if (!is.na(five_prime_utr_length_cor_plus) & !is.na(three_prime_utr_length_cor_plus)) {
                  tmp0_transcript$start <- as.numeric(tmp0_UTR5$start)
                  tmp0_transcript$end <- as.numeric(tmp0_UTR3$end)

                  # tmp_final <- rbind(tmp_final, tmp0_transcript)
                  il <- il + 1
                  tmp_final[[il]] <- tmp0_transcript
                } else if (is.na(five_prime_utr_length_cor_plus) & !is.na(three_prime_utr_length_cor_plus)) {
                  tmp0_transcript$start <- as.numeric(min(tmp$start))
                  tmp0_transcript$end <- as.numeric(tmp0_UTR3$end)

                  # tmp_final <- rbind(tmp_final, tmp0_transcript)
                  il <- il + 1
                  tmp_final[[il]] <- tmp0_transcript
                } else if (!is.na(five_prime_utr_length_cor_plus) & is.na(three_prime_utr_length_cor_plus)) {
                  tmp0_transcript$start <- as.numeric(tmp0_UTR5$start)
                  tmp0_transcript$end <- as.numeric(max(tmp$end))

                  # tmp_final <- rbind(tmp_final, tmp0_transcript)
                  il <- il + 1
                  tmp_final[[il]] <- tmp0_transcript
                } else {
                  tmp0_transcript$start <- as.numeric(min(tmp$start))
                  tmp0_transcript$end <- as.numeric(max(tmp$end))

                  # tmp_final <- rbind(tmp_final, tmp0_transcript)
                  il <- il + 1
                  tmp_final[[il]] <- tmp0_transcript
                }
              } else if (sen == "-") {
                # UTR3
                if (!is.na(three_prime_utr_length_cor_minus)) {
                  utr3_start <- as.numeric(min(tmp$start) - three_prime_utr_length_cor_minus)
                  if (utr3_start < 0) {
                    utr3_start <- 0
                  }

                  tmp0_UTR3$start <- utr3_start
                  tmp0_UTR3$end <- as.numeric(min(tmp$start - 1))

                  # tmp_final <- rbind(tmp_final, tmp0_UTR3)
                  il <- il + 1
                  tmp_final[[il]] <- tmp0_UTR3
                }


                # UTR5
                if (!is.na(five_prime_utr_length_cor_minus)) {
                  tmp0_UTR5$start <- as.numeric(max(tmp$end + 1))
                  tmp0_UTR5$end <- as.numeric(max(tmp$end) + five_prime_utr_length_cor_minus)

                  # tmp_final <- rbind(tmp_final, tmp0_UTR5)
                  il <- il + 1
                  tmp_final[[il]] <- tmp0_UTR5
                }


                # TRANSCRIPT
                if (!is.na(five_prime_utr_length_cor_minus) & !is.na(three_prime_utr_length_cor_minus)) {
                  tmp0_transcript$start <- as.numeric(tmp0_UTR3$start)
                  tmp0_transcript$end <- as.numeric(tmp0_UTR5$end)

                  # tmp_final <- rbind(tmp_final, tmp0_transcript)
                  il <- il + 1
                  tmp_final[[il]] <- tmp0_transcript
                } else if (is.na(five_prime_utr_length_cor_minus) & !is.na(three_prime_utr_length_cor_minus)) {
                  tmp0_transcript$start <- as.numeric(tmp0_UTR3$start)
                  tmp0_transcript$end <- as.numeric(max(tmp$end))

                  # tmp_final <- rbind(tmp_final, tmp0_transcript)
                  il <- il + 1
                  tmp_final[[il]] <- tmp0_transcript
                } else if (!is.na(five_prime_utr_length_cor_minus) & is.na(three_prime_utr_length_cor_minus)) {
                  tmp0_transcript$start <- as.numeric(min(tmp$start))
                  tmp0_transcript$end <- as.numeric(tmp0_UTR5$end)

                  # tmp_final <- rbind(tmp_final, tmp0_transcript)
                  il <- il + 1
                  tmp_final[[il]] <- tmp0_transcript
                } else {
                  tmp0_transcript$start <- as.numeric(min(tmp$start))
                  tmp0_transcript$end <- as.numeric(max(tmp$end))

                  # tmp_final <- rbind(tmp_final, tmp0_transcript)
                  il <- il + 1
                  tmp_final[[il]] <- tmp0_transcript
                }
              }
            }
          }
        }
      }
    }
  }



  tmp_final_2 <- rbindlist(tmp_final, use.names = TRUE, fill = TRUE)


  tmp_final_2$diff <- tmp_final_2$end - tmp_final_2$start
  tmp_final_2 <- tmp_final_2[tmp_final_2$diff > 10, ]

  tmp_final_2[, diff := NULL]


  final <- rbind(final, tmp_final_2)

  rm(tmp_final, tmp_final_2)

  output <- sort_alias(final)


  return(output)
}




#' Create Refflat-Style Annotations from Genomic Data
#'
#' @description
#' This function generates a refflat-style data frame from genomic annotations. The refflat format
#' represents gene structures, including transcripts and exons, for use in genomic analyses.
#' Data must be pre-loaded using load_annotation() followed by create_GTF_df().
#'
#' @param input A data frame containing genomic data with the following columns:
#'   - `chr`: Chromosome identifier
#'   - `start`: Start position of the annotation
#'   - `end`: End position of the annotation
#'   - `strand`: Strand information ('+' or '-')
#'   - `gene_name`: Name of the associated gene
#'   - `gene_id`: Unique identifier for the gene
#'   - `transcript_name`: Name of the associated transcript
#'   - `transcript_id`: Unique identifier for the transcript
#'   - `annotationType`: Type of annotation (e.g., 'EXON', 'TRANSCRIPT')
#'
#' @param geneName A string specifying the column name representing gene names (default: 'gene_name').
#' @param name A string specifying the column name representing gene IDs (default: 'gene_id').
#' @param genetic_elements Character vector (optional). A vector of genetic element types (e.g., 'CDS', 'GENE', 'MRNA') to include when creating transcripts for a RefFlat file.
#'
#' @return A data frame in refflat format with the following columns:
#'   - `geneName`: Gene name (based on the `geneName` parameter)
#'   - `name`: Gene ID (based on the `name` parameter)
#'   - `chrom`: Chromosome identifier
#'   - `strand`: Strand information
#'   - `txStart`: Start position of the transcript
#'   - `txEnd`: End position of the transcript
#'   - `cdsStart`: Start position of the coding sequence (CDS)
#'   - `cdsEnd`: End position of the CDS
#'   - `exonCount`: Number of exons in the transcript
#'   - `exonStarts`: Comma-separated list of exon start positions
#'   - `exonEnds`: Comma-separated list of exon end positions
#'
#' @details
#' The function processes input genomic data in parallel using multiple CPU cores. It filters and processes
#' annotations for transcripts and exons, ensuring that overlapping regions are resolved. Transcripts and
#' exons are matched, and the resulting structure is formatted in the refflat style.
#'
#' @examples
#' #
#' # Run the function
#' refflat_data <- refflat_create(input, geneName = "gene_name", name = "gene_id", genetic_elements = c("CDS", "GENE", "MRNA"))
#'
#' @import stringr dplyr doSNOW foreach parallel doParallel
#' @export
refflat_create <- function(input, geneName = "gene_name", name = "transcript_id",
                           genetic_elements = c("TRANSCRIPT", "MRNA", "CDS", "GENE")) {
  set.seed(123)
  input <- sort_alias(input)
  options(scipen = 999)
  iterations <- length(unique(input$chr))
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  CPU <- max(1, detectCores() - 2)
  cl <- makeCluster(CPU)
  registerDoParallel(cl)
  registerDoSNOW(cl)
  input <- input[, c(
    "chr", "start", "end", "gene_name", "gene_id",
    "transcript_name", "transcript_id", "strand", "annotationType"
  )]
  input <- distinct(input)
  input$sort_val <- 1:length(input$chr)
  input$diff <- input$end - input$start
  chromosomes <- unique(input$chr)
  cat("\n\n Refflat creating ... \n\n")
  results <- foreach(
    chr = chromosomes, .packages = c("dplyr"),
    .options.snow = opts, .combine = rbind
  ) %dopar% {
    df <- data.frame(matrix(ncol = 11, nrow = 0))
    colnames(df) <- c(
      "geneName", "name", "chrom", "strand",
      "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount",
      "exonStarts", "exonEnds"
    )
    for (sen in c("+", "-")) {
      tmp <- input[input$chr %in% chr & input$strand %in%
        sen, ]
      gen_list <- unique(tmp$gene_name)
      if (length(gen_list) == 0) {
        break
      }
      for (i in 1:length(gen_list)) {
        tmp_transcripts <- tmp[tmp$gene_name %in% gen_list[i], ]
        locus <- NULL
        for (g in genetic_elements) {
          if (g %in% toupper(tmp_transcripts$annotationType)) {
            locus <- g
            break
          }
        }
        if (is.null(locus)) {
          next
        }
        tmp_exons <- tmp_transcripts[toupper(tmp_transcripts$annotationType) %in%
          c("EXON"), ]
        tmp_transcripts <- tmp_transcripts[toupper(tmp_transcripts$annotationType) %in%
          c(locus), ]
        if (length(tmp_transcripts$chr) == 0) {
          next
        }
        for (tix in rownames(tmp_transcripts)) {
          corr_transcript <- tmp_transcripts[tix, ]
          curr_exons <- tmp_exons[tmp_exons$transcript_id %in%
            corr_transcript$transcript_id, ]
          if (nrow(curr_exons) > 0) {
            geneName <- gsub("\"", "", as.character(corr_transcript$gene_name[1]))
            name <- gsub("\"", "", as.character(corr_transcript$transcript_id[1]))
            chrom <- as.character(corr_transcript$chr[1])
            strand <- as.character(corr_transcript$strand[1])
            txStart <- format(as.numeric(corr_transcript$start[1]), scientific = FALSE)
            txEnd <- format(as.numeric(corr_transcript$end[1]), scientific = FALSE)
            cdsStart <- format(as.numeric(min(curr_exons$start)), scientific = FALSE)
            cdsEnd <- format(as.numeric(max(curr_exons$end)), scientific = FALSE)
            exonCount <- as.character(length(curr_exons$start))

            exonStarts <- c()

            for (es in curr_exons$start) {
              exonStarts <- c(exonStarts, format(as.numeric(es), scientific = FALSE))
            }

            exonStarts <- paste0(exonStarts, collapse = ",")

            rm(es)

            exonEnds <- c()

            for (ee in curr_exons$end) {
              exonEnds <- c(exonEnds, format(as.numeric(ee), scientific = FALSE))
            }

            exonEnds <- paste0(exonEnds, collapse = ",")

            rm(ee)

            df[nrow(df) + 1, ] <- c(
              geneName, name, chrom, strand,
              txStart, txEnd, cdsStart, cdsEnd,
              exonCount, exonStarts, exonEnds
            )
          } else {
            geneName <- as.character(corr_transcript$gene_name[1])
            name <- as.character(corr_transcript$transcript_id[1])
            if (length(name) == 0) {
              name <- as.character(corr_transcript$gene_id[1])
            }
            chrom <- as.character(corr_transcript$chr[1])
            strand <- as.character(corr_transcript$strand[1])
            txStart <- format(as.numeric(corr_transcript$start[1]), scientific = FALSE)
            txEnd <- format(as.numeric(corr_transcript$end[1]), scientific = FALSE)
            cdsStart <- format(as.numeric(corr_transcript$start[1]), scientific = FALSE)
            cdsEnd <- format(as.numeric(corr_transcript$end[1]), scientific = FALSE)
            exonCount <- 1

            exonStarts <- format(as.numeric(curr_exons$start[1]), scientific = FALSE)
            exonEnds <- format(as.numeric(curr_exons$end[1]), scientific = FALSE)

            df[nrow(df) + 1, ] <- c(
              as.character(gsub(
                "\"",
                "", geneName
              )), as.character(gsub(
                "\"",
                "", name
              )), chrom, strand, txStart,
              txEnd, cdsStart,
              cdsEnd, exonCount,
              exonStarts, exonEnds
            )
          }
        }
      }
    }
    return(df)
  }
  close(pb)
  stopCluster(cl)

  return(results)
}



#' Generate Full GTF Annotations
#'
#' @description
#' This function creates a GTF-style data frame by combining and formatting essential genomic annotation fields.
#' Data must be pre-loaded using load_annotation() followed by create_GTF_df() and / or add_UTR(), create_full_GTF(), or create_reduced_GTF().
#'
#' @param input A data frame containing genomic data with the following columns:
#'   - `chr`: Chromosome identifier
#'   - `start`: Start position of the annotation
#'   - `end`: End position of the annotation
#'   - `strand`: Strand information ('+' or '-')
#'   - `transcript_name`: Name of the associated transcript
#'   - `transcript_id`: Unique identifier for the transcript
#'   - `gene_name`: Name of the associated gene
#'   - `gene_id`: Unique identifier for the gene
#'   - `annotationType`: Type of annotation (e.g., 'EXON', 'TRANSCRIPT')
#'   - `source`: Source of the transcript
#'   - `gene_id`: Unique identifier for the gene
#'   - `gene_name`: Name of the associated gene
#'   - `transcript_name`: Name of the associated transcript
#'   - `transcript_id`: Unique identifier for the transcript
#'
#' @return A GTF-style data frame including the first eight columns of the input and a new `combine` column.
#' The `combine` column contains concatenated annotation fields in the format:
#' `gene_id`, `gene_name`, `transcript_name`, and `transcript_id`.
#'
#' @examples
#'
#' # Run the function
#' gtf_data <- create_full_GTF(input_data)
#'
#' @export
create_full_GTF <- function(input) {
  set.seed(123)

  options(scipen = 999)

  output <- input[, 1:8]

  gene_id <- gsub(" ", "", input$gene_id)
  gene_name <- gsub(" ", "", input$gene_name)
  transcript_name <- gsub(" ", "", input$transcript_name)
  transcript_id <- gsub(" ", "", input$transcript_id)
  gene_type <- gsub(" ", "", input$gene_type)


  combine_row <- function(id, name, tr_name, tr_id, gt) {
    parts <- c(
      if (nchar(id) > 0) paste0("gene_id ", id, ";"),
      if (nchar(name) > 0) paste0("gene_name ", name, ";"),
      if (nchar(tr_name) > 0) paste0("transcript_name ", tr_name, ";"),
      if (nchar(tr_id) > 0) paste0("transcript_id ", tr_id, ";"),
      if (nchar(gt) > 0) paste0("gene_type ", gt, ";")
    )
    paste(parts, collapse = " ")
  }

  output$combine <- mapply(combine_row, gene_id, gene_name, transcript_name, transcript_id, gene_type)
  return(output)
}




#' Generate Reduced GTF Annotations
#'
#' @description
#' This function creates a simplified GTF-style data frame by selecting key columns and adding annotation type and transcript source.
#' Data must be pre-loaded using load_annotation() followed by create_GTF_df() and / or add_UTR(), create_full_GTF(), or create_reduced_GTF().
#'
#' @param input A data frame containing genomic data with the following columns:
#'   - `chr`: Chromosome identifier
#'   - `start`: Start position of the annotation
#'   - `end`: End position of the annotation
#'   - `strand`: Strand information ('+' or '-')
#'   - `transcript_name`: Name of the associated transcript
#'   - `transcript_id`: Unique identifier for the transcript
#'   - `gene_name`: Name of the associated gene
#'   - `gene_id`: Unique identifier for the gene
#'   - `annotationType`: Type of annotation (e.g., 'EXON', 'TRANSCRIPT')
#'   - `source`: Source of the transcript
#'   - `gene_id`: Unique identifier for the gene
#'   - `gene_name`: Name of the associated gene
#'   - `transcript_name`: Name of the associated transcript
#'   - `transcript_id`: Unique identifier for the transcript
#'
#' @return A data frame containing the following columns:
#'   - `chr`: Chromosome identifier
#'   - `start`: Start position of the annotation
#'   - `end`: End position of the annotation
#'   - `strand`: Strand information
#'   - `transcript_name`: Transcript name
#'   - `transcript_id`: Transcript ID
#'   - `gene_name`: Gene name
#'   - `gene_id`: Gene ID
#'   - `annotationType`: Annotation type
#'   - `transcriptType`: Transcript source
#'
#' @examples
#'
#' # Run the function
#' reduced_gtf_data <- create_reduced_GTF(input_data)
#'
#' @export
create_reduced_GTF <- function(input) {
  set.seed(123)

  options(scipen = 999)

  output <- input %>%
    dplyr::select(chr, start, end, strand, transcript_name, transcript_id, gene_name, gene_id, gene_type)
  output$annotationType <- input$annotationType
  output$transcriptType <- input$source


  return(output)
}
