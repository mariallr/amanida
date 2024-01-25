#' Combine statistical results and compute vote-counting
#' 
#' \code{compute_amanida} Combines for the same entry or metabolite the statistical values of p-value and fold-change. Also is computed a vote-counting for each compound. Compound combination is made with PubChem CID when is available. 
#' 
#' Entries corresponding to metabolites are combined as follows:
#' \itemize{
#'  \item P-values are combined using Fisher method weighted by N and gamma distribution
#'  \item Fold-change are combined by weighted mean. Transformation works with fold-change transformed to log scale with base 2. 
#' }
#' 
#' Vote-counting is computed based on votes. Punctuation of entries is based on trend, up-regulation gives 1, down-regulation give -1 and equal behavior gives 0. Total sum is divided then by the total number of entries on each compound.
#'   
#' @param datafile data imported using \code{amanida_read} function w/o names checked by \code{check_names}
#' @param comp.inf include compounds IDs from PubChem, InChIKey, SMILES, KEGG, ChEBI, HMDB, Drugbank, Molecular Mass and Molecular Formula
#' @return METAtable S4 object with p-value combined, fold-change combined and vote-counting for each compound
#' 
#' @examples
#' \dontrun{
#' data("sample_data")
#' 
#' compute_amanida(sample_data)
#' }
#' 
#' @import dplyr
#' @import webchem
#' @importFrom stats qgamma pgamma
#' 
#' @export

compute_amanida <- function(datafile, comp.inf = F) {
  
  pvalue = NULL; foldchange = NULL; ratio = NULL; df = NULL; G = NULL;  
  pval = NULL; fc = NULL; N_total = NULL; reference = NULL; N = NULL;
  votes = NULL; articles = NULL; vote_counting = NULL; ref = NULL; trend = NULL;
  query = NULL; cid = NULL; CID = NULL; KEGG = NULL; ChEBI = NULL; HMDB = NULL;
  Drugbank = NULL; id_mod = NULL;
  
  set.seed(123)
  
  if(ncol(datafile) == 3) {
    stop("Compute_amanida needs quantitative data with p-value and fold-change. To import it use amanida_read in 'quan' mode.")
  }
  
  if (hasArg(comp.inf)) { 
     compinf <- comp.inf
    
  } else {
    compinf = FALSE
  }
  
  
  if ("cid" %in% colnames(datafile)) { 
  
    # Statistics grouping by compound identifier
    sta <- datafile |> group_by(id_mod) |>
      mutate(ratio = N /sum(N),
             df = n()*ratio,
             G = qgamma(pvalue, shape = df, scale = 2, lower.tail = F)) |>
      reframe(
        # Weigthed P-value combination
        pval = pgamma(sum(G), shape = n(), scale = 2, lower.tail = F),
        # Weighted average of fold-change
        fc = 2^(sum(log2(foldchange) * `N`) / sum(`N`)), N_total = sum(N),
        reference = paste(`ref`, collapse = "; "),
        id = unique(id_mod),
        cid = unique(cid)) |>
      group_by(id_mod) |>
      mutate(trend = case_when(fc < 1 ~ -1, T ~ 1)) |>
      select(c(`id`, `trend`, `pval`, `fc`, `N_total`, `reference`, `cid`))
    
    if(compinf == T) {
      b <- pc_prop(sta$cid, properties = c("MolecularFormula", "MolecularWeight", 
                                           "InChIKey", "CanonicalSMILES"))
      
      sta <- sta |> mutate(cid = as.integer(cid)) |>
        full_join(b, by = c("cid" = "CID"), relationship = "many-to-many") |>
        distinct() 
      
      if(requireNamespace("metaboliteIDmapping", quietly = TRUE)) {
        extra <- NULL
        for (i in 1:nrow(sta)){
          b <- metaboliteIDmapping::metabolitesMapping |> 
            mutate(CID = as.character(CID)) |>
            dplyr::filter(CID %in% sta$cid[i]) |> 
            slice(1) |>
            select(c(CID, KEGG, ChEBI, HMDB, Drugbank))
          extra <- extra |> bind_rows(b)
        }
        
        sta <- sta |> mutate(cid = as.character(cid),
                             id = tolower(unlist(id_mod))) |>
          full_join(extra, by = c("cid" = "CID"), relationship = "many-to-many") |>
          distinct() |>
          rename(PubChem_CID = cid) |>
          select(-reference)
        
      } else {
        msg <- c("metaboliteIDmapping is not installed. amanida can operate without metaboliteIDmapping, unless you want the complete information using comp.inf = F")
        warning(msg)
        
        sta <- sta |> mutate(cid = as.character(cid),
                             id = tolower(unlist(id_mod))) |>
          distinct() |>
          rename(PubChem_CID = cid) |>
          select(-reference)
      }
    } else {
      sta <- sta |> mutate(cid = as.character(cid),
                           id = tolower(unlist(id_mod))) |>
        distinct() |>
        rename(PubChem_CID = cid) |>
        select(-reference)
    }
    
  
  ## Vote-counting per each compound id
    vote <- datafile |>
      dplyr::group_by(`id_mod`) |> 
      reframe(
        id = tolower(unlist(id_mod)),
        # Votes per compound
        votes = sum(`trend`),
        # Number of reports
        articles = n(),
        # Vote-counting
        vote_counting = `votes`/`articles`
      ) |> group_by(id_mod) |> distinct()
  
  } else {
    sta <- datafile %>% group_by(id) |>
      mutate(ratio = N /sum(N),
             df = n()*ratio,
             G = qgamma(pvalue, shape = df, scale = 2, lower.tail = F)) |>
      reframe(
        # Weigthed P-value combination
        pval = pgamma(sum(G), shape = n(), scale = 2, lower.tail = F),
        # Weighted average of fold-change
        fc = 2^(sum(log2(foldchange) * `N`) / sum(`N`)), N_total = sum(N),
        reference = paste(`ref`, collapse = "; "),
        id = unique(id)) |>
      group_by(id) |>
      mutate(trend = case_when(fc < 1 ~ -1, T ~ 1)) |>
      select(c(`id`, `trend`, `pval`, `fc`, `N_total`, `reference`))
  
    ## Vote-counting per each compound id
    vote <- datafile |>
      dplyr::group_by(`id`) |> 
      reframe(
        id = id,
        # Votes per compound
        votes = sum(`trend`),
        # Number of reports
        articles = n(),
        # Vote-counting
        vote_counting = `votes`/`articles`
      ) |> 
      group_by(id) |> 
      distinct()
  }

  # Save results in S4 object and return
  METAtables(stat=sta, vote=vote)
}

