
#' Report
#' 
#' \code{amanida_report} creates a report for data using amanida functions
#'
#' In case of quantitative analysis \code{amanida_report} uses the functions \code{amanida_read} and \code{compute_amanida} for analyse the input data. Then the results are showed using \code{volcano_plot}, \code{explore_plot} and \code{vote_plot}.
#' 
#' @param input_file path to file
#' @param separator indicate the separator used in the input_file parameter
#' @param column_id columns names to use. It have to be in order identification, p-values, fold-changes, sample size and reference. 
#' @param analysis_type indicate if data will be quantitative or qualitative. Options are:
#' \itemize{
#'   \item "quan" for quantitative meta-analysis using p-value and fold-change
#'   \item "qual" for qualitative meta-analysis using trend label
#'   }
#' @param pvalue_cutoff numeric value to consider statistical significance
#' @param fc_cutoff numeric value to consider significance for effect size
#' @param votecount_lim minimum numeric value for vote-counting visualization
#' @return an html document
#' 
#' @import rmarkdown
#' 
#' @examples
#' column_id = c("Compound Name", "P-value", "Fold-change", "N total", "References")
#' input_file <- system.file("extdata", "dataset2.csv", package = "amanida")
#' 
#' amanida_report(input_file, separator = ";", column_id, analysis_type = "quan", pvalue_cutoff = 0.05, 
#' fc_cutoff = 4, votecount_lim = 2)
#' 
#' @export


amanida_report <- function(input_file, separator = NULL, analysis_type, column_id, 
                           pvalue_cutoff, fc_cutoff, votecount_lim) {
  
  rmarkdown::render(
  input = system.file("rmd", "amanida_report.Rmd", package = "amanida"),
  output_file = "Amanida_report.html",
  output_dir = getwd(),
  params = list(
    file_name = input_file,
    separator = separator,
    analysis_type = analysis_type,
    column_id = column_id,
    pvalue_cutoff = pvalue_cutoff,
    fc_cutoff = fc_cutoff,
    votecount_lim = votecount_lim,
    show_code = FALSE
    )
  )
}
