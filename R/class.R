
# S4 object to save the data

#' An S4 class to return results from metamet function
#' @slot stat results for statistics combining p-values and fold-changes
#' @slot vote vote-counting for metabolites
#' 
#' 
#' @export

#setOldClass(c("grouped_df", "tbl_df", "tbl"))
.METAtables <- setClass("METAtables",
                        slots = list(
                          stat = "tbl_df",
                          vote = "tbl_df"
                        ),
                        contains = class(tibble())
)

METAtables <- function(sta, vote) {
  .METAtables(stat = as_tibble(sta), vote = as_tibble(vote))
}


#setMethod(
#  f = "show",
#  signature = "METAtables",
#  definition = function(object) {
#    cat("An object of class", class(object), "\n", sep = "")
#    invisible(NULL)
#  }
#)

