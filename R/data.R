## Load data

data.read <- function(file, separator, coln) {
  library(dplyr)
  csv <- grep("\\.csv", file) #mirem si es csv o excel per llegir el document
  excel <- grep("\\.xlsx|\\.xls", file)
  txt <- grep("\\.txt", file)
  if (length(csv) == 1) {
    datafile <- read.csv(file, sep = separator)
  } else if (length(excel) == 1) {
    datafile <- readxl::read_excel(file)
  } else if (length(txt) == 1){
    datafile <- read.table(file, sep = separator)
  } else {
    print("Not compatible format, try csv, excel or txt")
  }
  
  datafile <- dplyr::as_tibble(datafile)
  datafile <- datafile %>% select(all_of(coln)) #seleccionem només les columnes que tenen les dades que ens interesa. L'usuari ha de donar els noms de les columnes. 
  
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
