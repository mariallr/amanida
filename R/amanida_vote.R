
#' Semi-quantitative meta-analysis
#' 
#' \code{amanida_vote} performs vote-counting on qualitative data
#'
#' Note that \code{amanida_vote} skips rows with missing values or NA. 
#'
#' Formats compatible are csv, xlsx, xls or txt.
#' 
#' @param data datafile data imported using amanida_read function
#' @return METAtable S4 object with vote-counting for each compound on @slot vote
#' 
#' @import dplyr
#' 
#' @examples
#' coln = c("Compound Name", "Behaviour", "References")
#' input_file <- system.file("extdata", "dataset2.csv", package = "amanida")
#' data_votes <- amanida_read(input_file, mode = "qual", coln, separator = ";")
#' 
#' vote_result <- amanida_vote(data_votes)
#' 
#' @export
#' 

amanida_vote <- function(data) {
    . = NULL; votes = NULL; articles = NULL; vote_counting = NULL; trend = NULL;
  
  vote <-  data %>%
    group_by(`id`) %>% 
    summarize(
      # Votes per compound
      votes = sum(trend),
      # Number of reports
      articles = n(),
      # Vote-counting
      vote_counting = votes/articles
      )
 
   sta = tibble()
  
  # Save results in S4 object and return
  METAtables(stat=sta, vote=vote)

}
