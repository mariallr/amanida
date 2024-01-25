
#' Qualitative meta-analysis
#' 
#' \code{amanida_vote} performs vote-counting on qualitative data. 
#' 
#' Vote-counting is computed without trend division. Punctuation of entries is based on trend, up-regulation gives 1, down-regulation give -1 and equal behavior gives 0. Total sum is divided then by the total number of entries on each compound. Compound combination is made with PubChem CID when is available. 
#'
#' Note that \code{amanida_vote} skips rows with missing values or NA. 
#'
#' Formats compatible are csv, xlsx, xls or txt.
#' 
#' @param data data imported using \code{amanida_read} function w/o names checked by \code{check_names}
#' @return METAtable S4 object with vote-counting for each compound on @slot vote
#' @export
#' 
#' @import dplyr
#' @import webchem
#' 
#' @examples
#' \dontrun{
#' coln = c("Compound Name", "Behaviour", "References")
#' input_file <- system.file("extdata", "dataset2.csv", package = "amanida")
#' data_votes <- amanida_read(input_file, mode = "qual", coln, separator = ";")
#' 
#' vote_result <- amanida_vote(data_votes)
#' }
#' 
#' @export
#' 

amanida_vote <- function(data) {
    . = NULL; votes = NULL; articles = NULL; vote_counting = NULL; trend = NULL;
    id_mod = NULL;query = NULL; cid = NULL; CID = NULL; KEGG = NULL; 
    ChEBI = NULL; HMDB = NULL; Drugbank = NULL;
    
    set.seed(123)
    
  if ("cid" %in% colnames(data)) { 
    
    vote <- data |>
      dplyr::group_by(`id_mod`) |>
      summarize(
        # Votes per compound
        votes = sum(trend),
        # Number of reports
        articles = n(),
        # Vote-counting
        vote_counting = votes/articles, 
        cid = unique(cid)
      )
    
      b <- pc_prop(vote$cid, properties = c("MolecularFormula", "MolecularWeight", "InChIKey", "CanonicalSMILES"))
        
      vote <- vote |> mutate(cid = as.integer(cid)) |>
          full_join(b, by = c("cid" = "CID"), relationship = "many-to-many") |>
          distinct() 
      if(requireNamespace("metaboliteIDmapping", quietly = TRUE)) {
        extra <- NULL
        for (i in 1:nrow(vote)){
          b <- metaboliteIDmapping::metabolitesMapping |> 
            mutate(CID = as.character(CID)) |>
            dplyr::filter(CID %in% vote$cid[i]) |> 
            slice(1) |>
            select(c(CID, KEGG, ChEBI, HMDB, Drugbank))
          extra <- extra |> bind_rows(b)
        }
        vote <- vote |> mutate(cid = as.character(cid), 
                               id = tolower(unlist(id_mod))) |>
          full_join(extra, by = c("cid" = "CID"), relationship = "many-to-many") |>
          distinct() |>
          rename(PubChem_CID = cid)
      } else {
        msg <- c("metaboliteIDmapping is not installed. amanida can operate without metaboliteIDmapping, unless you want the complete information using comp.inf = F")
        warning(msg)
        
        vote <- vote |> mutate(cid = as.character(cid),
                               id = tolower(unlist(id_mod))) |>
          distinct() |>
          rename(PubChem_CID = cid)
      }
  } else {
    
    vote <- data |>
      dplyr::group_by(`id`) |>
      summarize(
        # Votes per compound
        votes = sum(trend),
        # Number of reports
        articles = n(),
        # Vote-counting
        vote_counting = votes/articles
      )
  }
  
  
 
   sta = tibble()
  
  # Save results in S4 object and return
  METAtables(stat=sta, vote=vote)

}

