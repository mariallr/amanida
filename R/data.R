## Load data

data.read <- function(file, coln, separator=NULL) {
  
  #' Import data
  #' 
  #'\code{data.read} imports the data and formats for metamet function
  #'
  #'Note that \code{data.read} skips rows with missing values or NA. 
  #'
  #'Formats compatible are csv, xlsx, xls or txt.
  #' 
  #' @param file path to file
  #' @param coln columns names to use
  #' @param separator the separator used on file
  #' @return tibble table with data imported
  #' @examples 
  #' data.read(file, separator, coln)
 
  
  VAR_NAMES <- c("id", "pvalue", "foldchange", "N", "ref")

  # Get file type
  ext <- tools::file_ext(file)
  
  # Read file using appropiate extension
  if (ext %in% c("csv", "tsv", "txt")) {
    stopifnot("Please, specify a separator."=!is.null(separator))
    
    datafile <- readr::read_delim(file, delim = separator) %>%
      # In some (specially spanish) locales, when the delimiter is ";", the decimal
      # point is ","; let's make sure here this is correct
      mutate(
        across(c("P-value", "Fold-change"), function(x) sub(",", ".", x, fixed = TRUE))
      )
  } else if (ext %in% c("xlsx", "xls")) {
    datafile <- readxl::read_excel(file)
  } else {
    stop("Format not compatible; try csv, tsv, excel or txt. Aborting.")
  }
  

  datafile %>%
    # Select columns with data needed and rename
    select(all_of(coln)) %>%
    rename_with(.cols = everything(), .fn = ~ VAR_NAMES) %>%
    mutate(
      # Make sure numeric things are numeric
      foldchange = as.numeric(foldchange),
      pvalue = as.numeric(pvalue),
      N = as.integer(N),
      # Inverted fold-change for negative values
      foldchange = case_when(
        foldchange < 0 ~ 1 / abs(foldchange),
        T ~ foldchange
        ),
      # Add trend column
      trend = case_when(
        foldchange < 1 ~ -1,
        foldchange == 1 ~ 0,
        T ~ 1
        )
    )
  
}
