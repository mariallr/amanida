#' Combine statistical results and compute vote-counting
#' 
#' \code{compute_amanida} Combines for the same entry or metabolite the statistical values of p-value and fold-change. Also is computed a vote-counting for each compound. 
#' 
#' Entries corresponding to metabolites are combined as follows:
#' \itemize{
#'  \item P-values are combined using Fisher method weighted by N and gamma distribution
#'  \item Fold-change are combined by weighted mean. Transformation works with fold-change transformed to log scale with base 2. 
#' }
#' 
#' Vote-counting is computed based on votes. Punctuation of entries is based on trend, up-regulation gives 1, down-regulation give -1 and equal behavior gives 0. Total sum is divided then by the total number of entries on each compound.
#'   
#' @param datafile data imported using amanida_read function
#' @return METAtable S4 object with p-value combined, fold-change combined and vote-counting for each compound
#' 
#' @examples
#' data("sample_data")
#' 
#' compute_amanida(sample_data)
#' 
#' @import dplyr
#' @importFrom stats qgamma pgamma
#' 
#' @export

compute_amanida <- function(datafile) {
  
  pvalue = NULL; foldchange = NULL; ratio = NULL; df = NULL; G = NULL;  
  pval = NULL; fc = NULL; N_total = NULL; reference = NULL; N = NULL;
  votes = NULL; articles = NULL; vote_counting = NULL; ref = NULL; trend = NULL;
  
  set.seed(123)
  
  if(ncol(datafile) == 3) {
    stop("Compute_amanida needs quantitative data with p-value and fold-change. To import it use amanida_read in 'quan' mode.")
  }
    # Statistics grouping by compound identifier
    sta <- datafile %>% group_by(id) %>%
      mutate(ratio = N /sum(N),
             df = n()*ratio,
             G = qgamma(pvalue, shape = df, scale = 2, lower.tail = F)) %>%
      summarise(
        # Weigthed P-value combination
        pval = pgamma(sum(G), shape = n(), scale = 2, lower.tail = F),
        # Weighted average of fold-change
        fc = 2^(sum(log2(foldchange) * `N`) / sum(`N`)), N_total = sum(N),
                reference = paste(`ref`, collapse = "; ")) %>%
      mutate(trend = case_when(fc < 1 ~ -1, T ~ 1)) %>%
      select(c(`id`, `trend`, `pval`, `fc`, `N_total`, `reference`))
      
    ## Vote-counting per each compound id
    vote <- datafile %>% 
      dplyr::group_by(`id`) %>% 
      summarize(
      # Votes per compound
      votes = sum(`trend`),
      # Number of reports
      articles = n(),
      # Vote-counting
      vote_counting = `votes`/`articles`
    )
    
  # Save results in S4 object and return
  METAtables(stat=sta, vote=vote)
}
