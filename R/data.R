## Load data

data.read <- function(file, separator, coln) {
  
  # Function: read data files
  # Arguments: filename, separator symbol, columns to import
  # Return: tibble table with data imported
  
  library(dplyr)
  
  # Search file extension
  
  csv <- grep("\\.csv", file)
  excel <- grep("\\.xlsx|\\.xls", file)
  txt <- grep("\\.txt", file)
  
  # Read file using appropiate extension
  if (length(csv) == 1) {
    datafile <- read_csv(file, delim = separator)
  } else if (length(excel) == 1) {
    datafile <- readxl::read_excel(file)
  } else if (length(txt) == 1) {
    datafile <- read_delim(file, delim = separator)
  } else {
    print("Not compatible format, try csv, excel or txt")
  }
  
  # Select columns with data needed
  datafile <- datafile %>% select(all_of(coln)) 
  
  # Rename columns
  colnames(datafile) <- c("id", "pvalue", "foldchange", "N", "ref")
  
  # Correct fold-change values
  
  datafile <- datafile %>% mutate(foldchange = case_when(
    foldchange < 0 ~ 1/abs(foldchange), 
    TRUE ~ foldchange
  ))
  
  # Add trend column
  
  datafile <- datafile %>% mutate(trend = case_when (
    foldchange < 1 ~ -1,
    foldchange > 1 ~ 1,
    foldchange == 1 ~ 0
  ))
  
  
  return(datafile)
  
}
