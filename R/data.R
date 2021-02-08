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
    datafile <- read.csv(file, sep = separator)
  } else if (length(excel) == 1) {
    datafile <- readxl::read_excel(file)
  } else if (length(txt) == 1) {
    datafile <- read.table(file, sep = separator)
  } else {
    print("Not compatible format, try csv, excel or txt")
  }
  
  # Convert data imported as tibble
  datafile <- dplyr::as_tibble(datafile)
  
  # Select columns with data needed
  datafile <- datafile %>% select(all_of(coln)) 
  
  # Rename columns
  colnames(datafile) <- c("id", "pvalue", "foldchange", "N", "ref")
  
  # Add trend column
  if (length(datafile$foldchange[datafile$foldchange < 0]) > 0) { 
    
    # Inverted fold-change for negative values 
    datafile$foldchange[datafile$foldchange < 0] <- 
      1 / abs(datafile$foldchange[datafile$foldchange < 0] )
    
    # Negative values
    datafile[datafile$foldchange < 1, "trend"] <- -1 
    
    # Positive values
    datafile[datafile$foldchange > 1, "trend"] <- 1 
    
    # No behaviour change  
    datafile[datafile$foldchange == 1, "trend"] <- 0
    
  }
  
  
  return(datafile)
  
}
