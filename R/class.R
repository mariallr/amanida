#' An S4 class to return results from metamet function
#' @slot stat results for statistics combining p-values and fold-changes
#' @slot vote vote-counting for metabolites
#' 
#' @import tibble
#' @importFrom methods callNextMethod
#' 
#' @export

METAtables <- setClass("METAtables",
                       slots = c(stat = "tbl_df",
                                 vote = "tbl_df"))

setMethod(f = "initialize", "METAtables",
          function(.Object, stat, vote) {
            .Object <- callNextMethod()
            .Object@stat <- stat
            .Object@vote <- vote
            .Object
          })
