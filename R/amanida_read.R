
#' Import data
#' 
#' \code{amanida_read} imports the data and formats for \code{compute_amanida} 
#' or \code{amanida_vote} functions
#'
#' Note that \code{amanida_read} skips rows with missing values or NA. Negatives values for fold-change are transformed to positive (1/value). 
#'
#' Formats compatible are csv, xlsx, xls or txt.
#' 
#' @param file path to file
#' @param mode indicate if data will be quantitative or qualitative. Options are:
#' \itemize{
#'   \item "quan" for quantitative meta-analysis using p-value and fold-change
#'   \item "qual" for qualitative meta-analysis using trend label
#'   }
#' @param coln columns names to use. It has to be in order identification, p-values, fold-changes, sample size and reference.
#' @param separator the separator used on file
#' @return tibble table with data imported
#' 
#' @import dplyr
#' @import readr
#' @import readxl
#' @importFrom stats complete.cases
#' 
#' @examples
#' coln <-  c("Compound Name", "P-value", "Fold-change", "N total", "References")
#' input_file <- getsampleDB()
#' datafile <- amanida_read(input_file, mode = "quan", coln, separator=";")
#' 
#' @export

amanida_read <- function(file, mode, coln, separator=NULL) {
  . = NULL; foldchange = NULL; pvalue = NULL; N = NULL; 
  
  set.seed(123)
  
  VAR_NAMES <- c('id', 'pvalue', 'foldchange', 'N', 'ref')

  # Get file type
  ext <- tools::file_ext(file)
  
  # Read file using appropiate extension
  if (ext %in% c("csv", "tsv", "txt")) {
    stopifnot("Please, specify a separator."=!is.null(separator))
    
    datafile <- readr::read_delim(file, delim = separator, col_types = readr::cols()) %>%
      # In some (specially spanish) locales, when the delimiter is ";", the decimal
      # point is ","; let's make sure here this is correct
      mutate(
        across(coln[2:3], function(x) sub(",", ".", x, fixed = TRUE))
      )
  } else if (ext %in% c("xlsx", "xls")) {
    datafile <- readxl::read_excel(file)
  } else {
    stop("Format not compatible; try csv, tsv, excel or txt. Aborting.")
  }
  
  misrow <- sum(!complete.cases(datafile))
  
  if (mode == "quan") {
    datafile <- datafile %>%
      # Select columns with data needed
      select(all_of(coln)) %>%
      # Only complete cases and rename
      filter(complete.cases(.)) %>%
      rename_with(.cols = everything(), .fn = ~ VAR_NAMES) %>%
      mutate(
        # Make sure numeric things are numeric
        `foldchange` = as.numeric(`foldchange`),
        `pvalue` = as.numeric(`pvalue`),
        `N` = as.integer(`N`),
        # Inverted fold-change for negative values
        `foldchange` = case_when(
          `foldchange` < 0 ~ 1 / abs(`foldchange`),
          T ~ foldchange
        ),
        # Add trend column
        trend = case_when(
          `foldchange` < 1 ~ -1,
          `foldchange` == 1 ~ 0,
          T ~ 1
        )
      )
  } else if (mode == "qual") {
    VAR_NAMES <- c('id', 'trend', 'ref')
    
    datafile <- datafile %>%
      # Select columns with data needed
      select(all_of(coln)) %>%
      # Only complete cases and rename
      filter(complete.cases(.)) %>%
      rename_with(.cols = everything(), .fn = ~ VAR_NAMES) %>%
      mutate(trend = case_when(
        tolower(trend) == "down" ~ -1,
        T ~ 1
      ))
  } else {
    stop("Please, indicate mode: 'quan' for quantitative and 'qual' for qualitative")
  }
  
  message(paste("Loaded dataset with", nrow(datafile),"rows which contains", 
                length(unique(datafile$id)), "different identifiers. There are",
                misrow, "rows skipped from original dataset because contained NA values."))
  return(datafile)
}

