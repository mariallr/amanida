
#' Amanida harmonization
#' 
#' \code{check_names} check the names to harmonize them to a common nomenclature. Valid names are: 
#' chemical name, InChI, InChIKey and SMILES.
#'
#' Note that \code{check_names} depends on `webchem`package and it slows down the process. 
#'
#' Formats compatible are \code{amanida_read} output
#' 
#' @param datafile data imported using \code{amanida_read} function
#' @return tibble table with data imported
#' 
#' @import dplyr
#' @import webchem
#' 
#' @examples
#' coln <-  c("Compound Name", "P-value", "Fold-change", "N total", "References")
#' input_file <- getsampleDB()
#' datafile <- amanida_read(input_file, mode = "quan", coln, separator=";")
#' 
#' datafile <- check_names(datafile)
#' 
#' @export
#' 
check_names <- function(datafile){
  
  datafile <- datafile |> mutate(cid = NA, id_mod = NA)
  
  for(i in 1:length(datafile$id)){
    a <- get_cid(datafile$id[i], 
                 from = "name",
                 domain = c("compound", "substance", "assay"), 
                 match = "first")
    if(is.na(a$cid)){
      a <- get_cid(datafile$id[i], 
                   from = "inchikey",
                   domain = c("compound", "substance", "assay"),
                   match = "first")
      if(is.na(a$cid)){
        a <- get_cid(datafile$id[i], 
                     from = "inchi",
                     domain = c("compound", "substance", "assay"), 
                     match = "first")
        if(is.na(a$cid)){
          a <- get_cid(datafile$id[i], 
                       from = "smiles",
                       domain = c("compound", "substance", "assay"), 
                       match = "first")
        }
      } 
    } 
    
    datafile$cid[i] <- a$cid
    
    datafile$id_mod[i] <- pc_synonyms(a$cid, 
                                      from = "cid", 
                                      match = "first")
  }
  
  datafile <- datafile |> 
    mutate("id_mod" = ifelse(is.na(id_mod), id, id_mod))
}
