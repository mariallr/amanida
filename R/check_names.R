
#' Amanida harmonization
#' 
#' \code{check_names} check the names to harmonize them to a common nomenclature. Valid names are: 
#' chemical name, InChI, InChIKey and SMILES.
#'
#' Note that \code{check_names} depends on `webchem`package and it slows down the process. 
#'
#' Formats compatible are \code{amanida_read} output
#' 
#' @param data data imported using \code{amanida_read} function
#' @return tibble table with data imported with PubChem ID retrieved
#' 
#' @import dplyr
#' @import webchem
#' 
#' @examples
#' \dontrun{
#' coln <-  c("Compound Name", "P-value", "Fold-change", "N total", "References")
#' input_file <- getsampleDB()
#' datafile <- amanida_read(input_file, mode = "quan", coln, separator=";")
#' 
#' data_checked <- check_names(datafile) 
#' }
#' 
#' @export

check_names <- function(data) {
  id_mod = NULL;
  
  data <- data |> mutate(cid = NA, id_mod = NA)
  
  for(i in 1:length(data$id)){
    a <- get_cid(data$id[i], 
                 from = "name",
                 domain = c("compound", "substance", "assay"), 
                 match = "first")
    if(is.na(a$cid)){
      a <- get_cid(data$id[i], 
                   from = "inchikey",
                   domain = c("compound", "substance", "assay"),
                   match = "first")
      if(is.na(a$cid)){
        a <- get_cid(data$id[i], 
                     from = "inchi",
                     domain = c("compound", "substance", "assay"), 
                     match = "first")
        if(is.na(a$cid)){
          a <- get_cid(data$id[i], 
                       from = "smiles",
                       domain = c("compound", "substance", "assay"), 
                       match = "first")
        }
      } 
    } 
    
    data$cid[i] <- a$cid
    
    data$id_mod[i] <- pc_synonyms(a$cid, 
                                  from = "cid", 
                                  match = "first")
  }
  
  data <- data |> 
    mutate("id_mod" = ifelse(is.na(id_mod), id, id_mod))
  
  return(data)
}
