
#' Report
#' 
#' \code{amanida_report} creates a report from the data using amanida functions
#'
#' This function uses directly the dataset to create a report with the meta-analysis results. In case of quantitative analysis \code{amanida_report} uses the functions \code{amanida_read} and \code{compute_amanida} for analyse the input data. Then the results are showed using \code{volcano_plot}, \code{explore_plot} and \code{vote_plot}.
#' 
#' @param input_file path to the original dataset in xlsx, xls, csv or txt format
#' @param separator indicate the separator used in the input_file parameter
#' @param column_id vector containing columns names to use. It has to be in order identification, p-values, fold-changes, sample size and reference. 
#' @param path path to the directory where html file is created, otherwise the file will be saved in a temporal folder
#' @param analysis_type indicate if data will be quantitative, qualitative or both. Options are:
#' \itemize{
#'   \item "quan-qual" for quantitative and qualitative meta-analysis
#'   \item "quan" for quantitative meta-analysis using p-value and fold-change
#'   \item "qual" for qualitative meta-analysis using trend label
#'   }
#' @param pvalue_cutoff numeric value to consider statistical significance
#' @param fc_cutoff numeric value to consider significance for effect size
#' @param votecount_lim minimum numeric value for vote-counting visualization
#' @param comp_inf name checking using information from public databases
#' @return an html document saved in the working directory
#' 
#' @import rmarkdown
#' @importFrom kableExtra kbl kable_styling scroll_box footnote
#' @import knitr
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' 
#' @examples
#' \dontrun{
#' column_id = c("Compound Name", "P-value", "Fold-change", "N total", "References")
#' input_file <- getsampleDB()
#' 
#' amanida_report(input_file, separator = ";", column_id, analysis_type = "quan", 
#'                 pvalue_cutoff = 0.05, fc_cutoff = 4, votecount_lim = 2, 
#'                 comp_inf = F)
#' }
#' 
#' @export

amanida_report <- function(input_file, separator = NULL, analysis_type = NULL, 
                           column_id, pvalue_cutoff = NULL, fc_cutoff = NULL, 
                           votecount_lim, path = NULL, comp_inf = NULL) {
  analysis = NULL;
  
  Sys.setlocale("LC_TIME", "C")
  
  if (hasArg(analysis_type)) { 
    analysis <- analysis_type
    
  } else {
    analysis <- "quan-qual"
  }
  
  if (hasArg(path)) {
    path <- path
  } else {
    path <- tempdir()
  }
  
  if(analysis == "quan-qual") {
    try()
    rmarkdown::render(
      input = system.file("rmd", "amanida_report_quanqual.Rmd", package = "amanida"),
      output_file = "Amanida_report.html",
      output_dir = path,
      params = list(
        file_name = input_file,
        separator = separator,
        analysis_type = analysis,
        column_id = column_id,
        pvalue_cutoff = pvalue_cutoff,
        fc_cutoff = fc_cutoff,
        votecount_lim = votecount_lim,
        comp_inf = comp_inf,
        show_code = FALSE
      ), 
      quiet = T
    )
  } else if(analysis == "quan") {
    rmarkdown::render(
      input = system.file("rmd", "amanida_report_quan.Rmd", package = "amanida"),
      output_file = "Amanida_report.html",
      output_dir = path,
      params = list(
        file_name = input_file,
        separator = separator,
        analysis_type = analysis,
        column_id = column_id,
        pvalue_cutoff = pvalue_cutoff,
        fc_cutoff = fc_cutoff,
        votecount_lim = votecount_lim,
        comp_inf = comp_inf,
        show_code = FALSE
      ), 
      quiet = T
    )
    
  } else if (analysis == "qual") {
    rmarkdown::render(
      input = system.file("rmd", "amanida_report_qual.Rmd", package = "amanida"),
      output_file = "Amanida_report_qualitative.html",
      output_dir = path,
      params = list(
        file_name = input_file,
        separator = separator,
        analysis_type = analysis,
        column_id = column_id,
        votecount_lim = votecount_lim,
        comp_inf = comp_inf,
        show_code = FALSE
      ), 
      quiet = T
    )
  } else {
    message("Please indicate analysis type, 'quan-qual for quantitative and qualitative, 'quan' for quantitative or 'qual' for qualitative")
  }
}
