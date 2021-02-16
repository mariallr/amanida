## Load data

data.read <- function(file, separator, coln) {
  
  #' Import data
  #' 
  #'\code{data.read} imports the data and formats for metamet function
  #'
  #'Note that \code{data.read} skips rows with missing values or NA. 
  #'
  #'Formats compatible are csv, xlsx, xls or txt.
  #' 
  #' @param file path to file
  #' @param separator the separator used on file
  #' @param coln columns names to use
  #' @return tibble table with data imported
  #' @examples 
  #' data.read(file, separator, coln)
 
  
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
