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
  
  
  if (length(coln) == 6) { #hi ha l'opció d'afegir la direcció dels canvis o no
    colnames(datafile) <- c("id", "pvalue", "foldchange", "N", "trend", "ref") #canviem el nom de les columnes
  } else { 
    colnames(datafile) <- c("id", "pvalue", "foldchange", "N", "ref")
    if (length(datafile$foldchange[datafile$foldchange < 0])) { #si tenim valors positius i negatius podem extreure el trend
      datafile[datafile$foldchange < 0, "trend"] <- -1 #valors negatius posem -1
      datafile[datafile$foldchange > 0, "trend"] <- 1 #valors positius 1
      datafile[datafile$foldchange == 0, "trend"] <- 0 #si no hi ha canvi 0
    } 
  }
  
  
  return(datafile)
  
}
