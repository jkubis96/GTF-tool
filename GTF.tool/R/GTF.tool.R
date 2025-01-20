#' Load and Filter Annotation Data (GTF/GFF)
#'
#' @description
#' This function reads an annotation file in GTF (Gene Transfer Format), GFF3, GFF2, or similar formats. It supports files from 
#' multiple sources such as NCBI, Ensembl, and GENCODE. The function parses the file contents and optionally filters the data 
#' based on specified genetic elements. It uses various libraries for efficient data manipulation and parallel computing.
#'
#' @param path Character. The file path to the annotation file (GTF, GFF3, GFF2, etc.) to be loaded. The file should be tab-delimited 
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
#' The function supports annotation files in GTF, GFF3, and GFF2 formats from widely used sources such as NCBI, Ensembl, and GENCODE. 
#' It uses `readr` for efficient file reading and supports filtering based on case-insensitive matching of genetic elements.
#' 
#' @examples
#' # Load an annotation file without filtering:
#' annotation_data <- load_annotation("path/to/annotation_file.gtf")
#'
#' # Load an annotation file and filter for genes and exons:
#' annotation_data <- load_annotation("path/to/annotation_file.gtf", genetic_elements = c("gene", "exon"))
#'
#' @import readr stringr dplyr doSNOW foreach doParallel
#' @export
load_annotation <- function(path, genetic_elements = NaN) {
  
  set.seed(123)
  
  options(scipen = 999)


  cat('\n\n Data loading... \n\n')


  GTF <- read_delim(path, delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE, comment = '#', 
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
                    ))
  
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
  
  
  
  if (is.na(chromosomes_queue)) {
    chromosomes_queue = unique(input[[chro_col]])
  }
  
  output <- data.frame()
  
  for (ch in chromosomes_queue) {
    for (s in c('+', '-')) {
      
      cat('\n',paste('CHROMOSOME:', ch , '| STRAND:', s,' sorting...                              '))
      
    
      tmp <- input[(input[[chro_col]] %in% ch & input[[starnd]] %in% s),]
      tmp <- tmp[order(tmp[[start_col]], tmp[[end_col]]),]
      output <- rbind(output, tmp)
      
    }
    
  }
  
  return(output)
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
#' @import stringr dplyr doSNOW foreach doParallel
#' @export
create_GTF_df <- function(input, optimize = TRUE, shift = 100000) {

    set.seed(123)
    
    options(scipen = 999)
    
    cat('\n\n GTF converting... \n\n')

    df <- input[,1:8]
    
   
    
    if (TRUE %in% unique(grepl('gene_name', input$X9))) {
      
    #GENECODE & https://www.ensembl.org/index.html
    #############################################################
    
      # gene id
      df$gene_id_check <- grepl('gene_id', input$X9)
      
      df$gene_id <- ifelse(
        grepl('gene_id', input$X9),
        gsub(' ', '', gsub('"', '',gsub(".*=", "", gsub(";.*", "", gsub(".*gene_id", "", input$X9))))),
        ""
      )
      
      
      df$gene_id_check <- df$gene_id != ""
      
      
      # gene name
      df$gene_name_check <- grepl('gene_name', input$X9)
      
      df$gene_name <- ifelse(
        grepl('gene_name', input$X9),
        gsub(' ', '', gsub('"', '', gsub(".*=", "", gsub(";.*", "", gsub(".*gene_name", "", input$X9))))),
        ""
      )
      
      
      df$gene_name_check <- df$gene_name != ""
      
      
      
      # transcript name
      df$transcript_name_check <- grepl('transcript_name', input$X9)
      
      df$transcript_name <- ifelse(
        grepl('transcript_name', input$X9),
        gsub(' ', '', gsub('"', '', gsub(".*=", "", gsub(".*?\\s", "", gsub(";.*", "", gsub(".*transcript_name", "", input$X9)))))),
        ""
      )
      
      
      df$transcript_name_check <- df$transcript_name != ""
      
      
      
      # transcript id
      df$transcript_id_check <-grepl('transcript_id', input$X9)
      
      df$transcript_id <- ifelse(
        grepl('transcript_id', input$X9),
        gsub(' ', '', gsub('"', '', gsub(".*=", "", gsub(".*?\\s", "", gsub(";.*", "", gsub(".*transcript_id", "", input$X9)))))),
        ""
      )
      
      df$transcript_id_check <- df$transcript_id != ""
      
      
    
    if (optimize) {
    
    df$gene_id[df$gene_id == ""] <- df$transcript_id[df$gene_id == ""]
    
    df <- df[!(df$gene_name_check == FALSE & df$gene_id_check == FALSE & df$transcript_id_check == FALSE),  ]
    
    df$gene_name[df$gene_name_check == FALSE] <-  df$gene_id[df$gene_name_check == FALSE ]
    

    df$transcript_id[df$transcript_id == ""] <- df$gene_id[df$transcript_id == ""]
    
    df$transcript_name[df$transcript_name_check == FALSE] <- df$gene_name[df$transcript_name_check == FALSE]
    
    #repaire gene names
    gene_names <- df[,colnames(df) %in% c('gene_name', 'gene_id', 'transcript_id', 'transcript_name')]
    gene_names <- distinct(gene_names)
   
    
    cat('\n\n Unifying the names of genes... \n\n')
    
    iterations <- length(gene_names$gene_name)
    pb <- txtProgressBar(max = iterations, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    df <- df %>% distinct()
    CPU <- detectCores() - 2
    
    cl <- makeCluster(CPU)
    
    
    registerDoParallel(cl)
    registerDoSNOW(cl)
    
    
   
  
    df2 <- foreach(n = 1:length(gene_names$gene_name), .options.snow = opts, .combine=rbind) %dopar% {
      tmp <- df[df$gene_id %in% gene_names$gene_id[n],] 
      tmp$gene_name <- gene_names$gene_name[n]
      tmp
      
    }
    
    close(pb)
    stopCluster(cl)  
    
    
    df <- distinct(df2)
    
    rm(df2)
    
    
    
    
    df$lp <- 1:length(df$X1)
    
    #repair duplicated genes names from different loci
    
    duplicated_names <- unique(df$gene_name[duplicated(df$gene_name)])
    
    
    if (length(duplicated_names) > 0) {
      dedup_df <- df[, colnames(df) %in% c('X1','X7', 'X3', 'X4', 'X5', 'gene_name', 'lp')]
      dedup_df <- dedup_df[dedup_df$gene_name %in% duplicated_names,]
      dedup_df$renamed <- gsub(' ', '', gsub('"', '', dedup_df$gene_name))
      dedup_df$renamed <- paste0(dedup_df$renamed, '.',  dedup_df$X1)
      
      genes <- unique(duplicated_names)
      global_df <- data.frame(matrix(ncol = length(colnames(dedup_df)), nrow = 0))
      cat('\n\n Duplicated genes repairing...             \n\n ')
      
      for (strand in c('+','-')) {
        tmp_strand <- dedup_df[dedup_df$X7 %in% strand,]
        for (gen in genes){
          tmp_ch <- tmp_strand[tmp_strand$gene_name %in% gen, ]
          for (c in unique(tmp_ch$X1)) {
            tmp <- tmp_ch[tmp_ch$X1 %in% c,]
            tmp <- tmp[order(tmp$X4), ]
            group <- 1
            for (i in 1:length(tmp$gene_name)) {
              cat('\r',paste('GENE:',gen ,'| STRAND:', strand , '| PROGRESS:', round((i/length(tmp$gene_name)*100),2), '%           '))
              
              if (length(tmp$gene_name) == 1) {
                
                tmp$group <- group
                global_df <- rbind(global_df, tmp)

              } else if (i  == 1) {
                
                tmp1 <- tmp[i,]
                tmp1$group <- group
                global_df <- rbind(global_df, tmp1)
                
              } else if (i <= length(tmp$gene_name)) {
                tmp1 <- tmp[i-1,]
             
                tmp2 <- tmp[i,]
                if (tmp2$X1[1] == tmp1$X1[1] & tmp2$X4[1] - shift >  tmp1$X5[1]) {
                  group <- group + 1
                  tmp2$group <- group
                  global_df <- rbind(global_df, tmp2)
                } else {
                  tmp2$group <- group
                  global_df <- rbind(global_df, tmp2)
                }
                
              } 
            }
          } 
        }
      }
      
      
      
      
      
      
      global_df <- global_df %>%
        group_by(renamed) %>%
        filter(mean(group) != 1) %>%
        ungroup()
      
      
      
      
      global_df$renamed <- paste0(global_df$renamed, '.', global_df$group)
      
      

      global_df <- global_df %>%
        group_by(gene_name, X1) %>%
        mutate(
          renamed = if (n_distinct(renamed) == 1) {
            substr(renamed, 1, nchar(renamed) - 2)
          } else {
            renamed
          }
        ) %>%
        ungroup()
      
      
      
      
      repaired_df <- data.frame()
      for (g in unique(global_df$renamed)) {
        tmp <- global_df[global_df$renamed %in% g,]
        if (length(unique(tmp$X7)) > 1) {
          tmp <- global_df[global_df$renamed %in% g,]
          tmp$renamed <- paste0(tmp$renamed, '-str.', tmp$X7)
          repaired_df <- rbind(repaired_df, tmp)
        } else {
          repaired_df <- rbind(repaired_df, tmp)
        }
      }
      
      
      
      global_df <- repaired_df
      
      rm(repaired_df)
      
      #rename genes
      
      for (new_name in 1:length(global_df$lp)) {
        cat('\r',paste('Gene renaming -> PROGRESS:', round((new_name/length(global_df$lp)*100),2), '%                            '))
        df$gene_name[df$lp == global_df$lp[new_name]] <- global_df$renamed[global_df$lp == global_df$lp[new_name]]
      }
      
      
      
      rm(global_df)
      
    }
    
    }
    
    df <- df[,!colnames(df) %in% c('lp', 'gene_id_check', 'gene_name_check', 'transcript_name_check', 'transcript_id_check')]
    
    } else if (TRUE %in% unique(grepl('=gene-', input$X9))) {
      
    #NCBI
    ################################################################
    
     
      
      # gene id
      df$gene_id_check <- grepl('GeneID', input$X9)
      
      df$gene_id <- ifelse(
        grepl('GeneID', input$X9),
        gsub(' ', '', gsub('"', '', gsub(".*?\\s", "", gsub(",.*", "", gsub(".*GeneID:", "", input$X9))))),
        ""
      )
      
      df$gene_id_check <- df$gene_id != ""
      
      
      # gene name 1
      df$gene_name_check <- grepl('=gene-', input$X9)
      
      df$gene_name <- ifelse(
        grepl('gene_name', input$X9),
        gsub(' ', '', gsub('"', '', gsub(";.*", "", gsub(".*=gene-", "", input$X9)))),
        ""
      )
      
      
      df$gene_name_check <- df$gene_name != ""
      
      
      # gene name 1
      df$gene_name_check2 <- grepl(';gene=', input$X9)
      
      df$gene_name2 <- ifelse(
        grepl('gene_name', input$X9),
        gsub(' ', '', gsub('"', '', gsub(";.*", "", gsub(".*;gene=", "", input$X9)))),
        ""
      )
      
      
      df$gene_name_check2 <- df$gene_name2 != ""
      
      # transcript name
      
      df$transcript_name <- NA
      
      
      # transcript id
      df$transcript_id_check <-grepl('Genbank:', input$X9)
      
      df$transcript_id <- ifelse(
        grepl('transcript_id', input$X9),
        gsub(' ', '', gsub('"', '', gsub(";.*", "", gsub(",.*", "", gsub(".*Genbank:", "", input$X9))))),
        ""
      )
      
      df$transcript_id_check <- df$transcript_id != ""
      
      
      if (optimize) {
      #
      
      
      df$gene_id[df$gene_id == ""] <- df$transcript_id[df$gene_id == ""]
      
      df$gene_name[df$gene_name_check == FALSE & df$gene_name_check2 == TRUE] <- df$gene_name2[df$gene_name_check2 == TRUE & df$gene_name_check == FALSE]
      
      df$gene_name_check[df$gene_name_check2 == TRUE] <- TRUE
      
      df <- df[!(df$gene_name_check == FALSE & df$gene_id_check == FALSE & df$transcript_id_check == FALSE),  ]
      
      df$gene_name[df$gene_name_check == FALSE] <-  df$gene_id[df$gene_name_check == FALSE ]
      
      df$transcript_id[df$transcript_id == ""] <- df$gene_id[df$transcript_id == ""]
      
      df$transcript_name <- df$transcript_id
      
      
      #repaire gene names
      gene_names <- df[,colnames(df) %in% c('gene_name', 'gene_id', 'transcript_id', 'transcript_name')]
      gene_names <- distinct(gene_names)
      
      
      cat('\n\n Unifying the names of genes... \n\n')
      
      iterations <- length(gene_names$gene_name)
      pb <- txtProgressBar(max = iterations, style = 3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
      
      df <- df %>% distinct()
      CPU <- detectCores() - 2
      
      cl <- makeCluster(CPU)
      
      
      registerDoParallel(cl)
      registerDoSNOW(cl)
      
      
      
      
      df2 <- foreach(n = 1:length(gene_names$gene_name), .options.snow = opts, .combine=rbind) %dopar% {
        tmp <- df[df$gene_id %in% gene_names$gene_id[n],] 
        tmp$gene_name <- gene_names$gene_name[n]
        tmp
        
      }
      
      close(pb)
      stopCluster(cl)  
      
      
      df <- distinct(df2)
      
      rm(df2)
      
      
      
      
      df$lp <- 1:length(df$X1)
      
      #repair duplicated genes names from different loci
      
      duplicated_names <- unique(df$gene_name[duplicated(df$gene_name)])
      
      
      if (length(duplicated_names) > 0) {
        dedup_df <- df[, colnames(df) %in% c('X1','X7', 'X3', 'X4', 'X5', 'gene_name', 'lp')]
        dedup_df <- dedup_df[dedup_df$gene_name %in% duplicated_names,]
        dedup_df$renamed <- gsub(' ', '', gsub('"', '', dedup_df$gene_name))
        dedup_df$renamed <- paste0(dedup_df$renamed, '.',  dedup_df$X1)
        
        genes <- unique(duplicated_names)
        global_df <- data.frame(matrix(ncol = length(colnames(dedup_df)), nrow = 0))
        cat('\n\n Duplicated genes repairing...             \n\n ')
        
        for (strand in c('+','-')) {
          tmp_strand <- dedup_df[dedup_df$X7 %in% strand,]
          for (gen in genes){
            tmp_ch <- tmp_strand[tmp_strand$gene_name %in% gen, ]
            for (c in unique(tmp_ch$X1)) {
              tmp <- tmp_ch[tmp_ch$X1 %in% c,]
              tmp <- tmp[order(tmp$X4), ]
              group <- 1
              for (i in 1:length(tmp$gene_name)) {
                cat('\r',paste('GENE:',gen ,'| STRAND:', strand , '| PROGRESS:', round((i/length(tmp$gene_name)*100),2), '%           '))
                
                if (length(tmp$gene_name) == 1) {
                  
                  tmp$group <- group
                  global_df <- rbind(global_df, tmp)
                  
                } else if (i  == 1) {
                  
                  tmp1 <- tmp[i,]
                  tmp1$group <- group
                  global_df <- rbind(global_df, tmp1)
                  
                } else if (i <= length(tmp$gene_name)) {
                  tmp1 <- tmp[i-1,]
                  
                  tmp2 <- tmp[i,]
                  if (tmp2$X1[1] == tmp1$X1[1] & tmp2$X4[1] - shift >  tmp1$X5[1]) {
                    group <- group + 1
                    tmp2$group <- group
                    global_df <- rbind(global_df, tmp2)
                  } else {
                    tmp2$group <- group
                    global_df <- rbind(global_df, tmp2)
                  }
                  
                } 
              }
            } 
          }
        }
        
        
        
        
        
        
        global_df <- global_df %>%
          group_by(renamed) %>%
          filter(mean(group) != 1) %>%
          ungroup()
        
        
        
        
        global_df$renamed <- paste0(global_df$renamed, '.', global_df$group)
        
        
        
        global_df <- global_df %>%
          group_by(gene_name, X1) %>%
          mutate(
            renamed = if (n_distinct(renamed) == 1) {
              substr(renamed, 1, nchar(renamed) - 2)
            } else {
              renamed
            }
          ) %>%
          ungroup()
        
        
        
        
        repaired_df <- data.frame()
        for (g in unique(global_df$renamed)) {
          tmp <- global_df[global_df$renamed %in% g,]
          if (length(unique(tmp$X7)) > 1) {
            tmp <- global_df[global_df$renamed %in% g,]
            tmp$renamed <- paste0(tmp$renamed, '-str.', tmp$X7)
            repaired_df <- rbind(repaired_df, tmp)
          } else {
            repaired_df <- rbind(repaired_df, tmp)
          }
        }
        
        
        
        global_df <- repaired_df
        
        rm(repaired_df)
        
        #rename genes
        
        for (new_name in 1:length(global_df$lp)) {
          cat('\r',paste('Gene renaming -> PROGRESS:', round((new_name/length(global_df$lp)*100),2), '%                            '))
          df$gene_name[df$lp == global_df$lp[new_name]] <- global_df$renamed[global_df$lp == global_df$lp[new_name]]
        }
        
        
        
        rm(global_df)
        
      }
      
      }
      
      df <- df[,!colnames(df) %in% c('lp', 'gene_id_check', 'gene_name_check', 'transcript_name_check', 'transcript_id_check')]
      
  
    } else if (TRUE %in% unique(grepl('gene:', input$X9)) & TRUE %in% unique(grepl('transcript:', input$X9))) {
      
      #CUSTOM
      ################################################################
      
    

      # gene name
      df$gene_name_check <- grepl('gene:', input$X9)
      
      df$gene_name <- ifelse(
        grepl('gene_name', input$X9),
        gsub(' ', '', gsub('"', '', gsub(";.*", "", gsub(".*gene:", "", input$X9)))),
        ""
      )
      
      
      df$gene_name_check <- df$gene_name != ""
      
      # transcript name
      df$transcript_name_check <- grepl('transcript:', input$X9)
      
      df$transcript_name <- ifelse(
        grepl('transcript_name', input$X9),
        gsub(' ', '', gsub('"', '', gsub(";.*", "", gsub(",.*", "", gsub(".*transcript:", "", input$X9))))),
        ""
      )
      

      
      # transcript id
      # gene id
      
      
      df$transcript_id <- df$transcript_name
      
      df$gene_id <- df$gene_name
      
      
      if (optimize) {
        
      
      df <- df[!(df$gene_name_check == FALSE & df$transcript_name_check == FALSE),  ]
      
      df$gene_name[df$gene_name_check == FALSE] <-  df$transcript_name[df$gene_name_check == FALSE ]
      df$transcript_name[df$transcript_name_check == FALSE] <-  df$gene_name[df$transcript_name_check == FALSE ]
      
      #
      
   
      
      
        
      #repaire gene names
      gene_names <- df[,colnames(df) %in% c('gene_name', 'gene_id', 'transcript_id', 'transcript_name')]
      gene_names <- distinct(gene_names)
      
      
      cat('\n\n Unifying the names of genes... \n\n')
      
      iterations <- length(gene_names$gene_name)
      pb <- txtProgressBar(max = iterations, style = 3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
      
      df <- df %>% distinct()
      CPU <- detectCores() - 2
      
      cl <- makeCluster(CPU)
      
      
      registerDoParallel(cl)
      registerDoSNOW(cl)
      
      
      
      
      df2 <- foreach(n = 1:length(gene_names$gene_name), .options.snow = opts, .combine=rbind) %dopar% {
        tmp <- df[df$gene_id %in% gene_names$gene_id[n],] 
        tmp$gene_name <- gene_names$gene_name[n]
        tmp
        
      }
      
      close(pb)
      stopCluster(cl)  
      
      
      df <- distinct(df2)
      
      rm(df2)
      
      
      
      
      df$lp <- 1:length(df$X1)
      
      #repair duplicated genes names from different loci
      
      duplicated_names <- unique(df$gene_name[duplicated(df$gene_name)])
      
      
      if (length(duplicated_names) > 0) {
        dedup_df <- df[, colnames(df) %in% c('X1','X7', 'X3', 'X4', 'X5', 'gene_name', 'lp')]
        dedup_df <- dedup_df[dedup_df$gene_name %in% duplicated_names,]
        dedup_df$renamed <- gsub(' ', '', gsub('"', '', dedup_df$gene_name))
        dedup_df$renamed <- paste0(dedup_df$renamed, '.',  dedup_df$X1)
        
        genes <- unique(duplicated_names)
        global_df <- data.frame(matrix(ncol = length(colnames(dedup_df)), nrow = 0))
        cat('\n\n Duplicated genes repairing...             \n\n ')
        
        for (strand in c('+','-')) {
          tmp_strand <- dedup_df[dedup_df$X7 %in% strand,]
          for (gen in genes){
            tmp_ch <- tmp_strand[tmp_strand$gene_name %in% gen, ]
            for (c in unique(tmp_ch$X1)) {
              tmp <- tmp_ch[tmp_ch$X1 %in% c,]
              tmp <- tmp[order(tmp$X4), ]
              group <- 1
              for (i in 1:length(tmp$gene_name)) {
                cat('\r',paste('GENE:',gen ,'| STRAND:', strand , '| PROGRESS:', round((i/length(tmp$gene_name)*100),2), '%           '))
                
                if (length(tmp$gene_name) == 1) {
                  
                  tmp$group <- group
                  global_df <- rbind(global_df, tmp)
                  
                } else if (i  == 1) {
                  
                  tmp1 <- tmp[i,]
                  tmp1$group <- group
                  global_df <- rbind(global_df, tmp1)
                  
                } else if (i <= length(tmp$gene_name)) {
                  tmp1 <- tmp[i-1,]
                  
                  tmp2 <- tmp[i,]
                  if (tmp2$X1[1] == tmp1$X1[1] & tmp2$X4[1] - shift >  tmp1$X5[1]) {
                    group <- group + 1
                    tmp2$group <- group
                    global_df <- rbind(global_df, tmp2)
                  } else {
                    tmp2$group <- group
                    global_df <- rbind(global_df, tmp2)
                  }
                  
                } 
              }
            } 
          }
        }
        
        
        
        
        
        
        global_df <- global_df %>%
          group_by(renamed) %>%
          filter(mean(group) != 1) %>%
          ungroup()
        
        
        
        
        global_df$renamed <- paste0(global_df$renamed, '.', global_df$group)
        
        
        
        global_df <- global_df %>%
          group_by(gene_name, X1) %>%
          mutate(
            renamed = if (n_distinct(renamed) == 1) {
              substr(renamed, 1, nchar(renamed) - 2)
            } else {
              renamed
            }
          ) %>%
          ungroup()
        
        
        
        
        repaired_df <- data.frame()
        for (g in unique(global_df$renamed)) {
          tmp <- global_df[global_df$renamed %in% g,]
          if (length(unique(tmp$X7)) > 1) {
            tmp <- global_df[global_df$renamed %in% g,]
            tmp$renamed <- paste0(tmp$renamed, '-str.', tmp$X7)
            repaired_df <- rbind(repaired_df, tmp)
          } else {
            repaired_df <- rbind(repaired_df, tmp)
          }
        }
        
        
        
        global_df <- repaired_df
        
        rm(repaired_df)
        
        #rename genes
        
        for (new_name in 1:length(global_df$lp)) {
          cat('\r',paste('Gene renaming -> PROGRESS:', round((new_name/length(global_df$lp)*100),2), '%                            '))
          df$gene_name[df$lp == global_df$lp[new_name]] <- global_df$renamed[global_df$lp == global_df$lp[new_name]]
        }
        
        
        
        rm(global_df)
        
      }
      
      }
      
      df <- df[,!colnames(df) %in% c('lp', 'gene_id_check', 'gene_name_check', 'transcript_name_check', 'transcript_id_check')]
      
      
    }
    
    
    
    colnames(df) <- c('chr','source','annotationType','start','end','score','strand','phase','gene_id','gene_name','transcript_name','transcript_id')

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
#' output_data <- add_UTR(input, five_prime_utr_length = 400, three_prime_utr_length = 800, genetic_elements = c("EXON", "CDS", 'TRANSCRIPT', 'MRNA'))
#'
#' @import stringr dplyr doSNOW foreach doParallel
#' @export
add_UTR <- function(input, five_prime_utr_length = 400, three_prime_utr_length = 800, genetic_elements = c("EXON", "CDS", 'TRANSCRIPT', 'MRNA')) {
  
  set.seed(123)
  
  options(scipen = 999)

  cat('\n\n UTRs sequence extending...             \n\n')
  
  tmp_final <- data.frame()
  final <- input
  
  # input$sort_val <- 1:length(input$chr)
  chromosomes <- unique(input$chr)
  
  CDS <- final[toupper(input$annotationType) %in% toupper(genetic_elements), ]
  
  
  for (chr in chromosomes) {
  
  tmp_primary <- CDS[CDS$chr %in% chr,]
  # tmp_primary <- tmp_primary[order(tmp_primary$start),]
  # tmp_primary$lp <- 1:length(tmp_primary$chr)
  
  
  for (sen in c('+','-')) {
  tmp <- tmp_primary[tmp_primary$strand %in% sen,]
  gen_list <- unique(tmp$gene_name)

  
  if (length(gen_list) > 0) {
    
  
  for (i in 1:length(gen_list)) {
          # debug
          # if (i == 12) {break}
  
          # tmp variables
    
          five_prime_utr_length_cor_minus = NaN
          five_prime_utr_length_cor_plus = NaN
          three_prime_utr_length_cor_minus = NaN
          three_prime_utr_length_cor_plus = NaN
            
          cat('\r',paste('LOC:',chr, '| STRAND:', sen , '| PROGRESS:', round((i/length(gen_list)*100),2), '%           '))
          
          if (i == 1) {
            
            # first element
            
            tmp1 <- tmp[tmp$gene_name %in% gen_list[i],]
            
            tmp1_min <- tmp1[tmp1$start == min(tmp1$start),]
            tmp1_min <- tmp1_min[1,]
            tmp1_max <- tmp1[tmp1$end == max(tmp1$end),]
            tmp1_max <- tmp1_max[1,]
            
            
            tmp2 <- tmp[tmp$gene_name %in% gen_list[i+1],] 
            
            
            if (min(tmp2$start) > (max(tmp1$end) + (five_prime_utr_length + three_prime_utr_length)+2)) {
              
              if (sen == '+') {
                five_prime_utr_length_cor_plus <- round(five_prime_utr_length,0)
                three_prime_utr_length_cor_plus <- round(three_prime_utr_length,0)
              } else {
                five_prime_utr_length_cor_minus <- round(five_prime_utr_length,0)
                three_prime_utr_length_cor_minus <- round(three_prime_utr_length,0)
              }
              
            } else if (min(tmp2$start) > max(tmp1$end) + (five_prime_utr_length + three_prime_utr_length)/2) {
              
              if (sen == '+') {
                five_prime_utr_length_cor_plus <- round(five_prime_utr_length,0)
                three_prime_utr_length_cor_plus <- round(three_prime_utr_length*0.5,0)
              } else {
                five_prime_utr_length_cor_minus <- round(five_prime_utr_length*0.5,0)
                three_prime_utr_length_cor_minus <- round(three_prime_utr_length,0)
              }
              
            } else if (min(tmp2$start) > max(tmp1$end) + (five_prime_utr_length + three_prime_utr_length)/3) {
              
              if (sen == '+') {
                five_prime_utr_length_cor_plus <- round(five_prime_utr_length,0)
                three_prime_utr_length_cor_plus <- round(three_prime_utr_length*0.3,0)
              } else {
                five_prime_utr_length_cor_minus <- round(five_prime_utr_length*0.3,0)
                three_prime_utr_length_cor_minus <- round(three_prime_utr_length,0)
                
              }
              
            } else if (min(tmp2$start) > max(tmp1$end)){
              
              if (sen == '+') {
                five_prime_utr_length_cor_plus <- five_prime_utr_length
                three_prime_utr_length_cor_plus <- round((min(tmp2$start) - max(tmp1$end) -2) / 2,0)
                if (three_prime_utr_length_cor_plus < 0) {
                  three_prime_utr_length_cor_plus <- 0
                }
              } else {
                five_prime_utr_length_cor_minus <- round((min(tmp2$start) - max(tmp1$end) -2) / 2, 0)
                if (five_prime_utr_length_cor_minus < 0) {
                  five_prime_utr_length_cor_minus <- 0
                }
                three_prime_utr_length_cor_minus <- three_prime_utr_length
                
              }
              
              
            } else {
              
              if (sen == '+') {
                five_prime_utr_length_cor_plus <- five_prime_utr_length
               
              } else {
                
                three_prime_utr_length_cor_minus <- three_prime_utr_length
                
              }
              
              
            }
            
            
            
            
            tmp0 <- tmp1
            tmp0 <- tmp0[tmp0$start == min(tmp0$start),]
            tmp0 <- tmp0[tmp0$end == max(tmp0$end),]
            tmp0 <- tmp0[1,]
            
            tmp0_UTR5 <- tmp0
            tmp0_UTR5$source <- 'JBIO-predicted'
            tmp0_UTR5$annotationType <- 'five_prime_UTR'
            tmp0_UTR5$score <- '.'
            tmp0_UTR5$phase <- '.'
            
            
            tmp0_UTR3 <- tmp0
            tmp0_UTR3$source <- 'JBIO-predicted'
            tmp0_UTR3$annotationType <- 'three_prime_UTR'
            tmp0_UTR3$score <- '.'
            tmp0_UTR3$phase <- '.'
            
            tmp0_transcript <- tmp0
            tmp0_transcript$source <- 'JBIO-predicted'
            tmp0_transcript$annotationType <- 'transcript'
            tmp0_transcript$score <- '.'
            tmp0_transcript$phase <- '.'
            
            if (sen == '+') {
              
              
              #UTR5
              if (!is.na(five_prime_utr_length_cor_plus)) {
                
                utr5_start <- as.integer(tmp1_min$start - five_prime_utr_length_cor_plus)
                if (utr5_start < 0) {
                  utr5_start <- 0
                }
                
                tmp0_UTR5$start <- utr5_start
                tmp0_UTR5$end <- as.integer(tmp1_min$start - 1)
                
                tmp_final <- rbind(tmp_final, tmp0_UTR5)

              }
              
              
              #UTR3
              if (!is.na(three_prime_utr_length_cor_plus)) {
                
                tmp0_UTR3$start <- as.integer(tmp1_max$end + 1)
                tmp0_UTR3$end <- as.integer(tmp1_max$end + three_prime_utr_length_cor_plus)
                
                tmp_final <- rbind(tmp_final, tmp0_UTR3)

              }
              
              
              #TRANSCRIPT
              
              if (!is.na(three_prime_utr_length_cor_plus) & !is.na(five_prime_utr_length_cor_plus)) {
                
              tmp0_transcript$start <- as.integer(tmp0_UTR5$start)
              tmp0_transcript$end <- as.integer(tmp0_UTR3$end)
              
              tmp_final <- rbind(tmp_final, tmp0_transcript)
              tmp <- rbind(tmp, tmp0_transcript)
              
              } else if (is.na(three_prime_utr_length_cor_plus) & !is.na(five_prime_utr_length_cor_plus)) {
                
                tmp0_transcript$start <- as.integer(tmp0_UTR5$start)
                tmp0_transcript$end <- as.integer(tmp1_max$end)
                
                tmp_final <- rbind(tmp_final, tmp0_transcript)
                tmp <- rbind(tmp, tmp0_transcript)
                
              } else if (!is.na(three_prime_utr_length_cor_plus) & is.na(five_prime_utr_length_cor_plus)) {
                
                tmp0_transcript$start <- as.integer(tmp1_min$start)
                tmp0_transcript$end <- as.integer(tmp0_UTR3$end)
                
                tmp_final <- rbind(tmp_final, tmp0_transcript)
                tmp <- rbind(tmp, tmp0_transcript)
                
              } else {
                
                tmp0_transcript$start <- as.integer(tmp1_min$start)
                tmp0_transcript$end <- as.integer(tmp1_max$end)
                
                tmp_final <- rbind(tmp_final, tmp0_transcript)
                tmp <- rbind(tmp, tmp0_transcript)
                
              }
              
              
            } else if (sen == '-') {
              
              #UTR3
              if (!is.na(three_prime_utr_length_cor_minus)) {
                
                utr3_start <- as.integer(tmp1_min$start - three_prime_utr_length_cor_minus)
                if (utr3_start < 0) {
                  utr3_start <- 0
                }
                
                tmp0_UTR3$start <- utr3_start
                tmp0_UTR3$end <- as.integer(tmp1_min$start - 1)
                
                tmp_final <- rbind(tmp_final, tmp0_UTR3)

              }
              
              
              #UTR5
              if (!is.na(five_prime_utr_length_cor_minus)) {
                
                tmp0_UTR5$start <- as.integer(tmp1_max$end + 1)
                tmp0_UTR5$end <- as.integer(tmp1_max$end + five_prime_utr_length_cor_minus)
                
                tmp_final <- rbind(tmp_final, tmp0_UTR5)

              }
              
              #TRANSCRIPT
              if (!is.na(five_prime_utr_length_cor_minus) & !is.na(three_prime_utr_length_cor_minus)) {
                
                tmp0_transcript$start <- as.integer(tmp0_UTR3$start)
                tmp0_transcript$end <- as.integer(tmp0_UTR5$end)
                
                tmp_final <- rbind(tmp_final, tmp0_transcript)
                tmp <- rbind(tmp, tmp0_transcript)
                
              } else if (is.na(five_prime_utr_length_cor_minus) & !is.na(three_prime_utr_length_cor_minus)) {
                
                tmp0_transcript$start <- as.integer(tmp0_UTR3$start)
                tmp0_transcript$end <- as.integer(tmp1_max$end)
                
                tmp_final <- rbind(tmp_final, tmp0_transcript)
                tmp <- rbind(tmp, tmp0_transcript)
                
              } else if (!is.na(five_prime_utr_length_cor_minus) & is.na(three_prime_utr_length_cor_minus)) {
                
                tmp0_transcript$start <- as.integer(tmp1_min$start)
                tmp0_transcript$end <- as.integer(tmp0_UTR5$end)
                
                tmp_final <- rbind(tmp_final, tmp0_transcript)
                tmp <- rbind(tmp, tmp0_transcript)
                
              } else {
                
                tmp0_transcript$start <- as.integer(tmp1_min$start)
                tmp0_transcript$end <- as.integer(tmp1_max$end)
                
                tmp_final <- rbind(tmp_final, tmp0_transcript)
                tmp <- rbind(tmp, tmp0_transcript)
                
              }
              
            }
            
            
          } else if (i < length(gen_list)) {
            
            # middle element
            
            tmp_p <- tmp[tmp$gene_name %in% gen_list[i-1],]
            
            tmp1 <- tmp[tmp$gene_name %in% gen_list[i],]
            
            
            tmp1_min <- tmp1[tmp1$start == min(tmp1$start),]
            tmp1_min <- tmp1_min[1,]
            tmp1_max <- tmp1[tmp1$end == max(tmp1$end),]
            tmp1_max <- tmp1_max[1,]
            
            
            tmp2 <- tmp[tmp$gene_name %in% gen_list[i+1],] 
            
            
            if (min(tmp2$start) > (max(tmp1$end) + (five_prime_utr_length + three_prime_utr_length)+2)) {
              
              three_prime_utr_length_cor_plus <- three_prime_utr_length
      
              five_prime_utr_length_cor_minus <- five_prime_utr_length
              
              ################################################################################################
              if (max(tmp_p$end) < (min(tmp1$start) - (five_prime_utr_length)+2) & sen == '+') {
                five_prime_utr_length_cor_plus <- round(five_prime_utr_length,0)
              } else if (max(tmp_p$end) < (min(tmp1$start) - (three_prime_utr_length)+2) & sen == '-') {
                three_prime_utr_length_cor_minus <- round(three_prime_utr_length,0)
                #
              } else if (max(tmp_p$end) < (min(tmp1$start) - (five_prime_utr_length )/2) & sen == '+') {
                five_prime_utr_length_cor_plus <- round(five_prime_utr_length*0.5,0)
              } else if (max(tmp_p$end) < (min(tmp1$start) - (three_prime_utr_length )/2) & sen == '-') {
                three_prime_utr_length_cor_minus <- round(three_prime_utr_length*0.5,0)
                #
              } else if (max(tmp_p$end) < (min(tmp1$start) - (five_prime_utr_length)/3) & sen == '+') {
                five_prime_utr_length_cor_plus <- round(five_prime_utr_length*0.3,0)
              } else if (max(tmp_p$end) < (min(tmp1$start) - (three_prime_utr_length)/3) & sen == '-') {
                three_prime_utr_length_cor_minus <- round(three_prime_utr_length*0.3,0)
                
              } else if (max(tmp_p$end) < min(tmp1$start)) {
                
                if (sen == '+') {
                  five_prime_utr_length_cor_plus <- min(tmp1$start) - max(tmp_p$end) -2
                  if (five_prime_utr_length_cor_plus < 0) {
                    five_prime_utr_length_cor_plus <- 0
                  }
                  
                } else {
                  three_prime_utr_length_cor_minus <- min(tmp1$start) - max(tmp_p$end) -2
                  if (three_prime_utr_length_cor_minus < 0) {
                    three_prime_utr_length_cor_minus <- 0
                  }
                  
                }
                
              }
              #################################################################################################
              
              
            } else if (min(tmp2$start) > max(tmp1$end) + (five_prime_utr_length + three_prime_utr_length)/2) {
              
              three_prime_utr_length_cor_plus <- round(three_prime_utr_length*0.5,0)
        
              five_prime_utr_length_cor_minus <- round(five_prime_utr_length*0.5,0)
              
              ################################################################################################
              if (max(tmp_p$end) > (min(tmp1$start) - (five_prime_utr_length)+2) & sen == '+') {
                five_prime_utr_length_cor_plus <- round(five_prime_utr_length,0)
              } else if (max(tmp_p$end) > (min(tmp1$start) - (three_prime_utr_length)+2) & sen == '-') {
                three_prime_utr_length_cor_minus <- round(three_prime_utr_length,0)
                #
              } else if (max(tmp_p$end) < (min(tmp1$start) - (five_prime_utr_length )/2) & sen == '+') {
                five_prime_utr_length_cor_plus <- round(five_prime_utr_length*0.5,0)
              } else if (max(tmp_p$end) < (min(tmp1$start) - (three_prime_utr_length )/2) & sen == '-') {
                three_prime_utr_length_cor_minus <- round(three_prime_utr_length*0.5,0)
                #
              } else if (max(tmp_p$end) < (min(tmp1$start) - (five_prime_utr_length)/3) & sen == '+') {
                five_prime_utr_length_cor_plus <- round(five_prime_utr_length*0.3,0)
              } else if (max(tmp_p$end) < (min(tmp1$start) - (three_prime_utr_length)/3) & sen == '-') {
                three_prime_utr_length_cor_minus <- round(three_prime_utr_length*0.3,0)
              } else if (max(tmp_p$end) < min(tmp1$start)) {
                
                if (sen == '+') {
                  five_prime_utr_length_cor_plus <- min(tmp1$start) - max(tmp_p$end) -2
                  if (five_prime_utr_length_cor_plus < 0) {
                    five_prime_utr_length_cor_plus <- 0
                  }
                  
                } else {
                  three_prime_utr_length_cor_minus <- min(tmp1$start) - max(tmp_p$end) -2
                  if (three_prime_utr_length_cor_minus < 0) {
                    three_prime_utr_length_cor_minus <- 0
                  }
                  
                }
                
              }
              #################################################################################################
              
              
              
            } else if (min(tmp2$start) > max(tmp1$end) + (five_prime_utr_length + three_prime_utr_length)/3) {
              
              three_prime_utr_length_cor_plus <- round(three_prime_utr_length*0.3,0)
              
              five_prime_utr_length_cor_minus <- round(five_prime_utr_length*0.3,0)
              
              ################################################################################################
              if (max(tmp_p$end) > (min(tmp1$start) - (five_prime_utr_length)+2) & sen == '+') {
                five_prime_utr_length_cor_plus <- round(five_prime_utr_length,0)
              } else if (max(tmp_p$end) < (min(tmp1$start) - (three_prime_utr_length)+2) & sen == '-') {
                three_prime_utr_length_cor_minus <- round(three_prime_utr_length,0)
                #
              } else if (max(tmp_p$end) < (min(tmp1$start) - (five_prime_utr_length )/2) & sen == '+') {
                five_prime_utr_length_cor_plus <- round(five_prime_utr_length*0.5,0)
              } else if (max(tmp_p$end) < (min(tmp1$start) + (three_prime_utr_length )/2) & sen == '-') {
                three_prime_utr_length_cor_minus <- round(three_prime_utr_length*0.5,0)
                #
              } else if (max(tmp_p$end) < (min(tmp1$start) - (five_prime_utr_length)/3) & sen == '+') {
                five_prime_utr_length_cor_plus <- round(five_prime_utr_length*0.3,0)
              } else if (max(tmp_p$end) < (min(tmp1$start) - (three_prime_utr_length)/3) & sen == '-') {
                three_prime_utr_length_cor_minus <- round(three_prime_utr_length*0.3,0)
              } else if (max(tmp_p$end) < min(tmp1$start)) {
                
                if (sen == '+') {
                  five_prime_utr_length_cor_plus <- min(tmp1$start) - max(tmp_p$end) -2
                  if (five_prime_utr_length_cor_plus < 0) {
                    five_prime_utr_length_cor_plus <- 0
                  }
                } else {
                  three_prime_utr_length_cor_minus <- min(tmp1$start) - max(tmp_p$end) -2
                  if (three_prime_utr_length_cor_minus < 0) {
                    three_prime_utr_length_cor_minus <- 0
                  }
                  
                }
                
              }
              #################################################################################################
              
              
            } else if (min(tmp2$start) > max(tmp1$end)) {
              
            
              
              three_prime_utr_length_cor_plus <- round((min(tmp2$start) - max(tmp1$end) - 2)/2,0)
                if (three_prime_utr_length_cor_plus < 0) {
                  three_prime_utr_length_cor_plus <- 0
                }
              
              five_prime_utr_length_cor_minus <- round((min(tmp2$start) - max(tmp1$end) - 2)/2,0)
                if (five_prime_utr_length_cor_minus < 0) {
                  five_prime_utr_length_cor_minus <- 0
                }

              ################################################################################################
              if (max(tmp_p$end) < (min(tmp1$start) - (five_prime_utr_length)+2) & sen == '+') {
                five_prime_utr_length_cor_plus <- round(five_prime_utr_length,0)
              } else if (max(tmp_p$end) < (min(tmp1$start) - (three_prime_utr_length)+2) & sen == '-') {
                three_prime_utr_length_cor_minus <- round(three_prime_utr_length,0)
                #
              } else if (max(tmp_p$end) < (min(tmp1$start) - (five_prime_utr_length )/2) & sen == '+') {
                five_prime_utr_length_cor_plus <- round(five_prime_utr_length*0.5,0)
              } else if (max(tmp_p$end) < (min(tmp1$start) - (three_prime_utr_length )/2) & sen == '-') {
                three_prime_utr_length_cor_minus <- round(three_prime_utr_length*0.5,0)
                #
              } else if (max(tmp_p$end) < (min(tmp1$start) - (five_prime_utr_length)/3) & sen == '+') {
                five_prime_utr_length_cor_plus <- round(five_prime_utr_length*0.3,0)
              } else if (max(tmp_p$end) < (min(tmp1$start) - (three_prime_utr_length)/3) & sen == '-') {
                three_prime_utr_length_cor_minus <- round(three_prime_utr_length*0.3,0)
              } else if (max(tmp_p$end) < min(tmp1$start)) {
                
                if (sen == '+') {
                  five_prime_utr_length_cor_plus <- min(tmp1$start) - max(tmp_p$end) -2
                  if (five_prime_utr_length_cor_plus < 0) {
                    five_prime_utr_length_cor_plus <- 0
                  }
                } else {
                  three_prime_utr_length_cor_minus <- min(tmp1$start) - max(tmp_p$end) -2
                  if (three_prime_utr_length_cor_minus < 0) {
                    three_prime_utr_length_cor_minus <- 0
                  }
                  
                }
                
              }
              #################################################################################################
              
            } else if (max(tmp_p$end) < min(tmp1$start)) {
              
              
              if (max(tmp_p$end) < (min(tmp1$start) - (five_prime_utr_length)+2) & sen == '+') {
                five_prime_utr_length_cor_plus <- round(five_prime_utr_length,0)
              } else if (max(tmp_p$end) < (min(tmp1$start) - (three_prime_utr_length)+2) & sen == '-') {
                three_prime_utr_length_cor_minus <- round(three_prime_utr_length,0)
                #
              } else if (max(tmp_p$end) < (min(tmp1$start) - (five_prime_utr_length )/2) & sen == '+') {
                five_prime_utr_length_cor_plus <- round(five_prime_utr_length*0.5,0)
              } else if (max(tmp_p$end) < (min(tmp1$start) - (three_prime_utr_length )/2) & sen == '-') {
                three_prime_utr_length_cor_minus <- round(three_prime_utr_length*0.5,0)
                #
              } else if (max(tmp_p$end) < (min(tmp1$start) - (five_prime_utr_length)/3) & sen == '+') {
                five_prime_utr_length_cor_plus <- round(five_prime_utr_length*0.3,0)
              } else if (max(tmp_p$end) < (min(tmp1$start) - (three_prime_utr_length)/3) & sen == '-') {
                three_prime_utr_length_cor_minus <- round(three_prime_utr_length*0.3,0)
              } else if (max(tmp_p$end) < min(tmp1$start)) {
                
                if (sen == '+') {
                  five_prime_utr_length_cor_plus <- min(tmp1$start) - max(tmp_p$end) -2
                  if (five_prime_utr_length_cor_plus < 0) {
                    five_prime_utr_length_cor_plus <- 0
                  }
                } else {
                  three_prime_utr_length_cor_minus <- min(tmp1$start) - max(tmp_p$end) -2
                  if (three_prime_utr_length_cor_minus < 0) {
                    three_prime_utr_length_cor_minus <- 0
                  }
                  
                }
                
              }
              
              
            }
            
            
            tmp0 <- tmp1
            tmp0 <- tmp0[tmp0$start == min(tmp0$start),]
            tmp0 <- tmp0[tmp0$end == max(tmp0$end),]
            tmp0 <- tmp0[1,]
            
            tmp0_UTR5 <- tmp0
            tmp0_UTR5$source <- 'JBIO-predicted'
            tmp0_UTR5$annotationType <- 'five_prime_UTR'
            tmp0_UTR5$score <- '.'
            tmp0_UTR5$phase <- '.'
            
            
            tmp0_UTR3 <- tmp0
            tmp0_UTR3$source <- 'JBIO-predicted'
            tmp0_UTR3$annotationType <- 'three_prime_UTR'
            tmp0_UTR3$score <- '.'
            tmp0_UTR3$phase <- '.'
            
            tmp0_transcript <- tmp0
            tmp0_transcript$source <- 'JBIO-predicted'
            tmp0_transcript$annotationType <- 'transcript'
            tmp0_transcript$score <- '.'
            tmp0_transcript$phase <- '.'
            
            if (sen == '+') {
              
              
              #UTR5
              if (!is.na(five_prime_utr_length_cor_plus)) {
                
                tmp0_UTR5$start <- as.integer(tmp1_min$start - five_prime_utr_length_cor_plus)
                tmp0_UTR5$end <- as.integer(tmp1_min$start - 1)
                
                tmp_final <- rbind(tmp_final, tmp0_UTR5)

              }
              
              
              #UTR3
              if (!is.na(three_prime_utr_length_cor_plus)) {
                
                tmp0_UTR3$start <- as.integer(tmp1_max$end + 1)
                tmp0_UTR3$end <- as.integer(tmp1_max$end + three_prime_utr_length_cor_plus)
                
                tmp_final <- rbind(tmp_final, tmp0_UTR3)

              }
              
              
              #TRANSCRIPT
              if (!is.na(three_prime_utr_length_cor_plus) & !is.na(five_prime_utr_length_cor_plus)) {
                
                tmp0_transcript$start <- as.integer(tmp0_UTR5$start)
                tmp0_transcript$end <- as.integer(tmp0_UTR3$end)
                
                tmp_final <- rbind(tmp_final, tmp0_transcript)
                tmp <- rbind(tmp, tmp0_transcript)
                
              } else if (is.na(three_prime_utr_length_cor_plus) & !is.na(five_prime_utr_length_cor_plus)) {
                
                tmp0_transcript$start <- as.integer(tmp0_UTR5$start)
                tmp0_transcript$end <- as.integer(tmp1_max$end)
                
                tmp_final <- rbind(tmp_final, tmp0_transcript)
                tmp <- rbind(tmp, tmp0_transcript)
                
              } else if (!is.na(three_prime_utr_length_cor_plus) & is.na(five_prime_utr_length_cor_plus)) {
                
                tmp0_transcript$start <- as.integer(tmp1_min$start)
                tmp0_transcript$end <- as.integer(tmp0_UTR3$end)
                
                tmp_final <- rbind(tmp_final, tmp0_transcript)
                tmp <- rbind(tmp, tmp0_transcript)
                
              } else {
                
                tmp0_transcript$start <- as.integer(tmp1_min$start)
                tmp0_transcript$end <- as.integer(tmp1_max$end)
                
                tmp_final <- rbind(tmp_final, tmp0_transcript)
                tmp <- rbind(tmp, tmp0_transcript)
                
              }
              
              
              
            } else if (sen == '-') {
              
              #UTR3
              if (!is.na(three_prime_utr_length_cor_minus)) {
                
                tmp0_UTR3$start <- as.integer(tmp1_min$start - three_prime_utr_length_cor_minus)
                tmp0_UTR3$end <- as.integer(tmp1_min$start - 1)
                
                tmp_final <- rbind(tmp_final, tmp0_UTR3)

              }
              
              
              #UTR5
              if (!is.na(five_prime_utr_length_cor_minus)) {
                
                tmp0_UTR5$start <- as.integer(tmp1_max$end + 1)
                tmp0_UTR5$end <- as.integer(tmp1_max$end + five_prime_utr_length_cor_minus)
                
                tmp_final <- rbind(tmp_final, tmp0_UTR5)

              }
                
              
              #TRANSCRIPT
              if (!is.na(five_prime_utr_length_cor_minus) & !is.na(three_prime_utr_length_cor_minus)) {
                
                tmp0_transcript$start <- as.integer(tmp0_UTR3$start)
                tmp0_transcript$end <- as.integer(tmp0_UTR5$end)
                
                tmp_final <- rbind(tmp_final, tmp0_transcript)
                tmp <- rbind(tmp, tmp0_transcript)
              
              } else if (is.na(five_prime_utr_length_cor_minus) & !is.na(three_prime_utr_length_cor_minus)) {
                
                tmp0_transcript$start <- as.integer(tmp0_UTR3$start)
                tmp0_transcript$end <- as.integer(tmp1_max$end)
                
                tmp_final <- rbind(tmp_final, tmp0_transcript)
                tmp <- rbind(tmp, tmp0_transcript)
                
              } else if (!is.na(five_prime_utr_length_cor_minus) & is.na(three_prime_utr_length_cor_minus)) {
                
                tmp0_transcript$start <- as.integer(tmp1_min$start)
                tmp0_transcript$end <- as.integer(tmp0_UTR5$end)
                
                tmp_final <- rbind(tmp_final, tmp0_transcript)
                tmp <- rbind(tmp, tmp0_transcript)
                
              } else {
                
                tmp0_transcript$start <- as.integer(tmp1_min$start)
                tmp0_transcript$end <- as.integer(tmp1_max$end)
                
                tmp_final <- rbind(tmp_final, tmp0_transcript)
                tmp <- rbind(tmp, tmp0_transcript)
                
              }
              
            }
            
          } else if (i == length(gen_list)) {
            
            # final element
            
            tmp1 <- tmp[tmp$gene_name %in% gen_list[i],]
            
            tmp1_min <- tmp1[tmp1$start == min(tmp1$start),]
            tmp1_min <- tmp1_min[1,]
            tmp1_max <- tmp1[tmp1$end == max(tmp1$end),]
            tmp1_max <- tmp1_max[1,]
            
            tmp2 <- tmp[tmp$gene_name %in% gen_list[i-1],] 
            
            
            # End UTR
            
            if (max(tmp2$end) < max(tmp1$end)) {
              
              three_prime_utr_length_cor_plus <- round(three_prime_utr_length,0)
              five_prime_utr_length_cor_minus <- round(five_prime_utr_length,0)
              
            }
              
            
            # tmp2 -> tmp1 +
            if (max(tmp2$end) < (min(tmp1$start) - (five_prime_utr_length )+2) & sen == '+') {
              
                five_prime_utr_length_cor_plus <- round(five_prime_utr_length,0)
              
            } else if (max(tmp2$end) < (min(tmp1$start) - (five_prime_utr_length)/2) & sen == '+') {
              
                five_prime_utr_length_cor_plus <- round(five_prime_utr_length*0.5,0)
              
            } else if (max(tmp2$end) < (min(tmp1$start) - (five_prime_utr_length)/3) & sen == '+') {
              
                five_prime_utr_length_cor_plus <- round(five_prime_utr_length*0.3,0)
              
            } else if (max(tmp2$end) < min(tmp1$start) & sen == '+') {
              
                five_prime_utr_length_cor_plus <- round((min(tmp1$start) - max(tmp2$end) -2), 0)
                if (five_prime_utr_length_cor_plus < 0) {
                  five_prime_utr_length_cor_plus <- 0
                }
              
            }
            
            
            # tmp2 -> tmp1 -
            if (max(tmp2$end) < (min(tmp1$start) - (three_prime_utr_length)+2) & sen == '-') {
              
                three_prime_utr_length_cor_minus <- round(three_prime_utr_length,0)
              
            } else if (max(tmp2$end) < (min(tmp1$start) - (three_prime_utr_length)/2) & sen == '-') {
              
                three_prime_utr_length_cor_minus <- round(three_prime_utr_length*0.5,0)
              
            } else if (max(tmp2$end) < (min(tmp1$start) - (three_prime_utr_length)/3) & sen == '-') {
              
                three_prime_utr_length_cor_minus <- round(three_prime_utr_length*0.3,0)
              
            } else if (max(tmp2$end) < min(tmp1$start) & sen == '-') {
              
                three_prime_utr_length_cor_minus <- round((min(tmp1$start) - max(tmp2$end) -2) / 2, 0)
                if (three_prime_utr_length_cor_minus < 0) {
                  three_prime_utr_length_cor_minus <- 0
               
                }
              
            }
            
            
            
            tmp0 <- tmp1
            tmp0 <- tmp0[tmp0$start == min(tmp0$start),]
            tmp0 <- tmp0[tmp0$end == max(tmp0$end),]
            tmp0 <- tmp0[1,]
            
            tmp0_UTR5 <- tmp0
            tmp0_UTR5$source <- 'JBIO-predicted'
            tmp0_UTR5$annotationType <- 'five_prime_UTR'
            tmp0_UTR5$score <- '.'
            tmp0_UTR5$phase <- '.'
            
            
            tmp0_UTR3 <- tmp0
            tmp0_UTR3$source <- 'JBIO-predicted'
            tmp0_UTR3$annotationType <- 'three_prime_UTR'
            tmp0_UTR3$score <- '.'
            tmp0_UTR3$phase <- '.'
            
            tmp0_transcript <- tmp0
            tmp0_transcript$source <- 'JBIO-predicted'
            tmp0_transcript$annotationType <- 'transcript'
            tmp0_transcript$score <- '.'
            tmp0_transcript$phase <- '.'
            
            if (sen == '+') {
              
              
              #UTR5
              if (!is.na(five_prime_utr_length_cor_plus)) {
                
                utr5_start <- as.integer(tmp1_min$start - five_prime_utr_length_cor_plus)
                if (utr5_start < 0) {
                  utr5_start <- 0
                }
                
                tmp0_UTR5$start <- utr5_start
                tmp0_UTR5$end <- as.integer(tmp1_min$start - 1)
                
                tmp_final <- rbind(tmp_final, tmp0_UTR5)

              }
              
              
              
              #UTR3
              if (!is.na(three_prime_utr_length_cor_plus)) {
                
                tmp0_UTR3$start <- as.integer(tmp1_max$end + 1)
                tmp0_UTR3$end <- as.integer(tmp1_max$end + three_prime_utr_length_cor_plus)
                
                tmp_final <- rbind(tmp_final, tmp0_UTR3)

              }
              
              
              #TRANSCRIPT
              if (!is.na(five_prime_utr_length_cor_plus) & !is.na(three_prime_utr_length_cor_plus)) {
                
                tmp0_transcript$start <- as.integer(tmp0_UTR5$start)
                tmp0_transcript$end <- as.integer(tmp0_UTR3$end)
                
                tmp_final <- rbind(tmp_final, tmp0_transcript)
                tmp <- rbind(tmp, tmp0_transcript)
              
              } else if (is.na(five_prime_utr_length_cor_plus) & !is.na(three_prime_utr_length_cor_plus)) {
                
                tmp0_transcript$start <- as.integer(tmp1_min$start)
                tmp0_transcript$end <- as.integer(tmp0_UTR3$end)
                
                tmp_final <- rbind(tmp_final, tmp0_transcript)
                tmp <- rbind(tmp, tmp0_transcript)
                
              } else if (!is.na(five_prime_utr_length_cor_plus) & is.na(three_prime_utr_length_cor_plus)) {
                
                tmp0_transcript$start <- as.integer(tmp0_UTR5$start)
                tmp0_transcript$end <- as.integer(tmp1_max$end)
                
                tmp_final <- rbind(tmp_final, tmp0_transcript)
                tmp <- rbind(tmp, tmp0_transcript)
                
              } else {
                
                tmp0_transcript$start <- as.integer(tmp1_min$start)
                tmp0_transcript$end <- as.integer(tmp1_max$end)
                
                tmp_final <- rbind(tmp_final, tmp0_transcript)
                tmp <- rbind(tmp, tmp0_transcript)
                
              }
              
              
              
            } else if (sen == '-') {
              
              #UTR3
              if (!is.na(three_prime_utr_length_cor_minus)) {
                
                utr3_start <- as.integer(tmp1_min$start - three_prime_utr_length_cor_minus)
                if (utr3_start < 0) {
                  utr3_start <- 0
                }
              
                tmp0_UTR3$start <- utr3_start
                tmp0_UTR3$end <- as.integer(tmp1_min$start - 1)
                
                tmp_final <- rbind(tmp_final, tmp0_UTR3)

              }
                
              
              #UTR5
              if (!is.na(five_prime_utr_length_cor_minus)) {
                
                tmp0_UTR5$start <- as.integer(tmp1_max$end + 1)
                tmp0_UTR5$end <- as.integer(tmp1_max$end + five_prime_utr_length_cor_minus)
                
                tmp_final <- rbind(tmp_final, tmp0_UTR5)

              }
              
              
              #TRANSCRIPT
              if (!is.na(five_prime_utr_length_cor_minus) & !is.na(three_prime_utr_length_cor_minus)) {
                
                tmp0_transcript$start <- as.integer(tmp0_UTR3$start)
                tmp0_transcript$end <- as.integer(tmp0_UTR5$end)
                
                tmp_final <- rbind(tmp_final, tmp0_transcript)
                tmp <- rbind(tmp, tmp0_transcript)
              
              } else if (is.na(five_prime_utr_length_cor_minus) & !is.na(three_prime_utr_length_cor_minus)) {
                
                tmp0_transcript$start <- as.integer(tmp0_UTR3$start)
                tmp0_transcript$end <- as.integer(tmp1_max$end)
                
                tmp_final <- rbind(tmp_final, tmp0_transcript)
                tmp <- rbind(tmp, tmp0_transcript)
                
              } else if (!is.na(five_prime_utr_length_cor_minus) & is.na(three_prime_utr_length_cor_minus)) {
                
                tmp0_transcript$start <- as.integer(tmp1_min$start)
                tmp0_transcript$end <- as.integer(tmp0_UTR5$end)
                
                tmp_final <- rbind(tmp_final, tmp0_transcript)
                tmp <- rbind(tmp, tmp0_transcript)
                
              } else {
                
                tmp0_transcript$start <- as.integer(tmp1_min$start)
                tmp0_transcript$end <- as.integer(tmp1_max$end)
                
                tmp_final <- rbind(tmp_final, tmp0_transcript)
                tmp <- rbind(tmp, tmp0_transcript)
                
              }
                
            }
            
            
          } 
          

        } 
      }
    }
  }

  
  final <- rbind(final, tmp_final)
  
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
#
#' # Run the function
#' refflat_data <- refflat_create(input, geneName = 'gene_name', name = 'gene_id')
#'
#' @import stringr dplyr doSNOW foreach doParallel
#' @export
refflat_create <- function(input, geneName = 'gene_name', name = 'gene_id') {
  
  set.seed(123)
  
  options(scipen = 999)

  iterations <- length(unique(input$chr))
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  
  CPU <- detectCores() - 2
  
  cl <- makeCluster(CPU)

  
  
  registerDoParallel(cl)
  registerDoSNOW(cl)
  
  
  
  
  input <- input[,c('chr', 'start', 'end' ,'gene_name', 'gene_id', 'transcript_name', 'transcript_id', 'strand', 'annotationType')]
  input$annotationType[toupper(input$annotationType) %in% c('MRNA', 'GENE', 'TRANSCRIPT')] <- 'TRANSCRIPT'
  input <- distinct(input)
  
  input$sort_val <- 1:length(input$chr)
  input$diff <- input$end - input$start
  

  chromosomes <- unique(input$chr)

  
  df <- data.frame(matrix(ncol = 11, nrow = 0))
  colnames(df) <- c('geneName', 'name', 'chrom', 'strand', 'txStart','txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds')
  
  cat('\n\n Refflat creating ... \n\n')
  results <- foreach(chr = chromosomes, .packages = c('dplyr'), .options.snow = opts, .combine=rbind) %dopar% {
    
    tmp_primary <- input[input$chr %in% chr,]
    tmp_primary <- tmp_primary[order(tmp_primary$start),]
    tmp_primary$lp <- 1:length(tmp_primary$chr)
    
    
    for (sen in c('+','-')) {
      tmp <- tmp_primary[tmp_primary$strand %in% sen,]
      gen_list <- unique(tmp$gene_name)
      
      if (length(gen_list) == 1) {break}
      
      for (i in 1:length(gen_list)) {
        

          check_vector <- toupper(tmp$annotationType[tmp$gene_name %in% gen_list[i]])
          if ('EXON' %in% check_vector & 'TRANSCRIPT' %in% check_vector) {
          tmp_tx <- tmp[tmp$gene_name %in% gen_list[i],]
          tmp_tx <- tmp_tx[toupper(tmp_tx$annotationType) == 'TRANSCRIPT',]
          tmp_dif <- tmp_tx[tmp_tx$diff == max(tmp_tx$diff),]
          min_tr <- min(tmp_dif$start)
          max_tr <- max(tmp_dif$end)
          tmp_tx <- tmp_tx[tmp_tx$diff == max(tmp_tx$diff) | tmp_tx$start < min_tr | tmp_tx$end > max_tr ,]
          tmp_ex <- tmp[tmp$gene_name %in% gen_list[i] & toupper(tmp$annotationType) == 'EXON',]

          for (trancript in 1:nrow(tmp_tx)) {
            tmp_tx_tmp <- tmp_tx[trancript,]
            tmp_ex_tmp <- tmp_ex[tmp_ex$start >= tmp_tx_tmp$start & tmp_ex$end <= tmp_tx_tmp$end,]
            
            decission <- TRUE
            trn <- 0
            while (decission) {
              trn <- trn + 1
                dec <- c()
                start <- c()
                end <- c()
                before <- 0
                for(ex in 1:nrow(tmp_ex_tmp)) {
                  if (before <= tmp_ex_tmp$start[ex]) {
                    start <- c(start, as.integer(tmp_ex_tmp$start[ex]))
                    end <- c(end, as.integer(tmp_ex_tmp$end[ex]))
                    before  <-  as.integer(tmp_ex_tmp$end[ex])
                    
                   
      
                  } else if (before > as.integer(tmp_ex_tmp$start[ex])) {
                    
                    dec <- c(dec, ex-1)
                   
                  }
                }
                
                if (length(dec) > 0) {
                tmp_ex_tmp <- tmp_ex_tmp[-dec[1], ]
                decission <- TRUE
                
                } else if (length(dec) == 0) {
                  decission <- FALSE
                  
                }
                
                
                df[nrow(df) + 1,] <- c(as.character(gsub('"', '', tmp_tx_tmp[[geneName]][1])), as.character(gsub('"', '', tmp_ex_tmp[[name]][1])), as.character(tmp_ex_tmp$chr[1]), as.character(tmp_ex_tmp$strand[1]), as.integer(tmp_tx_tmp$start[1]), as.integer(tmp_tx_tmp$end[1]), as.integer(min(start)), as.integer(max(end)), as.integer(length(start)), sub('"', '', paste0(start, collapse = ',')), sub('"', '', paste0(end, collapse = ',')))
                df <- distinct(df)
                
                
                } 
                  
                 
                
                }
              
            
                }
          
          
              }
            }
                results <- df
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

  output <- input[,1:8]
  gene_id <- paste0('gene_id ',gsub(' ', '',input$gene_id), ';')
  gene_name <-  paste0('gene_name ', gsub(' ','',input$gene_name), ';')
  transcript_name <- paste0('transcript_name ', gsub(' ', '',input$transcript_name), ';')
  transcript_id <- paste0('transcript_id ', gsub(' ', '',input$transcript_id), ';')
  output$combine <- paste0(gene_id, gene_name, transcript_name, transcript_id)
  
  
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
    dplyr::select(chr,start,end,strand,transcript_name,transcript_id, gene_name, gene_id)
  output$annotationType <- input$annotationType
  output$transcriptType <- input$source
  
  
  return(output)  
  
}

