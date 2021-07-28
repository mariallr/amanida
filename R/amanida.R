
#' amanida: A package for Meta-Analysis with non-integral data
#' 
#' @title amanida
#' 
#' @name amanida
#' 
#' @author Maria Llambrich, Eudald Correig and Raquel Cumeras
#' 
#' Results combination for meta-analysis using only significance and effect size. 
#' \itemize{
#'  \item P-values and fold-change are combined to obtain a global significance on each metabolite.
#'  \item Produces a volcano plot summarizing the relevant results from meta-analysis.
#'  \item Qualitative meta-analysis for metabolites
#'  \item Graphical representation of qualitative analysis by bar plot
#'  \item Trend explore plot to detect discrepancies between studies at a first glance
#' }
#' 
#' 
#' @import dplyr
#' @import ggplot2
#' @import ggrepel
#' @import readr
#' @import readxl
#' @import tidyr
#' @importFrom methods hasArg new
#' @importFrom stats qgamma pgamma reorder
#' @importFrom magrittr %>%
#' @importFrom kableExtra kbl kable_styling scroll_box footnote
#' 
#' @docType package
#' 

NULL
