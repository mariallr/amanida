#' Combine statistical results and compute vote-counting
#' 
#' \code{compute_amanida} Combines for the same entry or metabolite the statistical values of p-value and fold-change. Also is computed a vote-counting for each compound. 
#' 
#' Entries corresponding to metabolites are divided by trend and then combined as follows:
#' \itemize{
#'  \item P-values are combined using Fisher method weighted by N and chi-squared distribution
#'  \item Fold-change are combined by weighted mean
#' }
#' 
#' Vote-counting is computed without trend division. Punctuation of entries is based on trend, up-regulation gives 1, down-regulation give -1 and equal behavior gives 0. Total sum is divided then by the total number of entries on each compound.
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
#' 
#' @export

compute_amanida <- function(datafile) {
  
  pvalue = NULL; foldchange = NULL; trend = NULL; N = NULL; logp = NULL; 
  chisq = NULL; logfc = NULL; ref = NULL; pval = NULL; fc = NULL; N_total = NULL;
  reference = NULL; votes = NULL; articles = NULL; vote_counting = NULL;

    # Statistics grouping by compound identifier and trend
    sta <- datafile %>% 
      mutate(logp = log10(`pvalue`),
             logfc = log2(`foldchange`)) %>%
      group_by(`id`, `trend`) %>% 
      summarize(
      # Combine p-values using Fisher's method weighted by number of individuals
      chisq = (-2/sum(`N`))*sum(`logp` * `N`),
      # P-values comparison using Chi-squared distribution
      pval = pchisq(`chisq`, 2*n(), lower.tail = F),
      # Wheigthed mean for combining fold-change values
      fc = 2^(sum(`logfc` * `N`) / sum(`N`)),
      # Sum of total individuals
      N_total = sum(`N`),
      # References
      reference = paste(`ref`, collapse = "; ")
      ) %>%
      select(c(`id`, `trend`, `pval`, `fc`, `N_total`, `reference`)) %>%
      ungroup()
    
    ## Vote-counting per each compound id
    vote <- datafile %>% 
      group_by(`id`) %>% 
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
