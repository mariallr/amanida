## Load data


VAR_NAMES <- c("id", "pvalue", "foldchange", "N", "ref")


data.read <- function(file, coln, separator=NULL) {
  # Function: read data files
  # Arguments: filename, separator symbol, columns to import
  # Return: tibble table with data imported
  
  library(dplyr)
  
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
