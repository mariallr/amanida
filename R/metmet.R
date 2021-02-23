
## Statistics

#' Combine statistical results and compute vote-counting
#' 
#' \code{metamet} Combines for the same entry or metabolite the statistical values of p-value and fold-change. Also is computed a vote-counting for each compound. 
#' 
#' Entries corresponding to metabolites are divided by trend and then combined as follows:
#' \itemize{
#'  \item P-values are combined using Fisher method weighted by N and chi-squared distribution
#'  \item Fold-change are combined by weighted mean
#' }
#' 
#' Vote-counting is computed without trend division. Punctuation of entries is based on trend, up-regulation gives 1, down-regulation give -1 and equal behavior gives 0. Total sum is divided then by the total number of entries on each compound.
#'   
#' @param datafile data imported using data.read function
#' @return METAtable S4 object with p-value combined, fold-change combined and vote-counting for each compound
#' @examples 
#' coln = c("Compound Name", "P-value", "Fold-change", "N total", "References")
#' datafile <- data.read("~/OneDrive - URV/metaanalysis_pack/dataset2.xlsx", coln)
#' res.met <- metmet(datafile)
#' 
#' @import dplyr
#' @importFrom magrittr %>%
#' 
#' @export

metmet <- function(datafile) {

    # Statistics grouping by compound identifier and trend
    sta <- datafile %>% 
      mutate(logp = log10(pvalue),
             logfc = log(foldchange)) %>%
      group_by(id, trend) %>% 
      summarize(
      # Combine p-values using Fisher's method weighted by number of individuals
      chisq = (-2/sum(N))*sum(logp * N),
      # P-values comparison using Chi-squared distribution
      pval = pchisq(chisq, 2*n(), lower.tail = F),
      # Wheigthed mean for combining fold-change values
      fc = 2^(sum(logfc * N) / sum(N)),
      # Sum of total individuals
      N_total = sum(N),
      # References
      reference = paste(ref, collapse = "; ")
      ) %>%
      select(c(id, trend, pval, fc, N_total, reference))
    
    ## Vote-counting per each compound id
    vote <- datafile %>% 
      group_by(id) %>% 
      summarize(
      # Votes per compound
      votec = sum(trend),
      # Number of reports
      articles = n(),
      # Vote-counting
      VC = votec/articles
    )
    
  # Save results in S4 object
  #mets <- initialize("METAtables", stat = stat,  vote = vote)
  mets <- METAtables(sta, vote)
  
  return(mets)
}
