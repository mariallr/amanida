## Load data


VAR_NAMES <- c("id", "pvalue", "foldchange", "N", "ref")


data.read <- function(file, separator, coln) {
  # Function: read data files
  # Arguments: filename, separator symbol, columns to import
  # Return: tibble table with data imported
  
  library(dplyr)
  
  # Read file using appropiate extension
  if (length(grep("\\.csv|\\.txt", file)) > 0) {
    datafile <- readr::read_delim(file, sep = separator)
  } else if (length(grep("\\.xlsx|\\.xls", file)) > 0) {
    datafile <- readxl::read_excel(file)
  } else {
    stop("Format not compatible; try csv, excel or txt. Aborting.")
  }
  
  datafile %>%
    # Select columns with data needed and rename
    select(all_of(coln)) %>%
    rename_with(.cols = everything(), .fn = ~ VAR_NAMES) %>%
    mutate(
      # Inverted fold-change for negative values
      foldchange = case_when(foldchange < 0 ~ 1 / abs(foldchange),
                             T ~ foldchange),
      # Add trend column
      trend = case_when(foldchange < 1 ~ -1,
                        foldchange == 1 ~ 0,
                        T ~ 1)
    )
  
}
