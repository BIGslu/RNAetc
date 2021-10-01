#' Clean alignment metrics into a data frame
#'
#' @param dir Character. Path to directory containing metric files
#' @param flagstat Logical. Include samtools flagstat metrics?
#' @param picard Logical. Include Picard RNAseq metrics?
#'
#' @return Data frame with alignment quality metrics
#' @export
#'

clean_metrics <- function(dir="data_raw/metrics/", flagstat=TRUE, picard=TRUE){
  LIBRARY<-READ_GROUP<-SAMPLE<-X1<-libID<-name<-value<-NULL

  #Read in flagstat results
  if(flagstat){
    files1 <- list.files(path=dir, pattern="flagstat", full.names = TRUE, recursive = TRUE)
    flag.df <- data.frame()

    for(f1 in files1){
      #Original alignments
      if(!grepl("filter", f1)){
        temp1 <- readr::read_tsv(f1, col_names = FALSE, show_col_types = FALSE) %>%
          #Separate data to column
          tidyr::separate(X1, into=c("value","name"), sep=" \\+ 0 ", fill="right") %>%
          dplyr::mutate(value = as.numeric(value)) %>%
          tidyr::drop_na(value) %>%
          #Recode variable names
          tidyr::separate(name, into=c("name"), sep="[(]", extra="drop") %>%
          dplyr::mutate(name = forcats::fct_recode(factor(name),
                                          to.align="in total ",
                                          align="mapped ",
                                          secondary.align="secondary",
                                          chimeric.align="supplementary",
                                          PCR.dups="duplicates",
                                          paired="paired in sequencing",
                                          R1.paired="read1", R2.paired="read2",
                                          both.align.paired= "with itself and mate mapped",
                                          one.align.paired="singletons " ,
                                          both.align.paired.diffCHR="with mate mapped to a different chr",
                                          both.align.paired.diffCHR.mapq="with mate mapped to a different chr ")) %>%
          #add libID
          dplyr::mutate(libID = gsub("_Aligned_flagstat.tsv","",basename(f1)))

        flag.df <- dplyr::bind_rows(flag.df, temp1)
      } else {
        #Filtered alignments
        temp2 <- readr::read_tsv(f1, col_names = FALSE, show_col_types = FALSE) %>%
          #Separate data to column
          tidyr::separate(X1, into=c("value","name"), sep=" \\+ 0 ", fill="right") %>%
          dplyr::mutate(value = as.numeric(value)) %>%
          #Keep filtered align total
          dplyr::filter(name == "in total (QC-passed reads + QC-failed reads)") %>%
          dplyr::mutate(name = "align.filtered") %>%
          #add libID
          dplyr::mutate(libID = gsub("_filter_flagstat.tsv","",basename(f1)))

        flag.df <- dplyr::bind_rows(flag.df, temp2)
      }
    }

    #format to wide
    flag.df.format <- flag.df %>%
      tidyr::pivot_wider() %>%
      dplyr::select(libID, dplyr::everything())
  } else { flag.df.format <- NULL }

  #Read in Picard results
  if(picard){
    files2 <- list.files(path=dir, pattern="picard", full.names = TRUE, recursive = TRUE)
    pic.df <- data.frame()

    for(f2 in files2){
      temp3 <- readr::read_tsv(f2, skip=6, n_max=1, show_col_types = FALSE) %>%
        dplyr::select(-c(SAMPLE,LIBRARY,READ_GROUP)) %>%
        #add libID
        dplyr::mutate(libID = gsub("_Aligned_picard.tsv","",basename(f2)))

      pic.df <- dplyr::bind_rows(pic.df, temp3)
    }

    pic.df.format <- pic.df %>%
      dplyr::select(libID, dplyr::everything())

  } else { pic.df.format <- NULL }

  metrics.all <- dplyr::full_join(flag.df.format, pic.df.format, by="libID")
  return(metrics.all)
}
