

source("R/data.R")
source("R/plots.R")

## Statistics

# S4 object to save the data

#' An S4 class to return results from metamet function
#' @slot up results for up regulated metabolites
#' @slot down results for down regulated metabolites
#' @slot equal results for no trend provided or no change
#' @slot vote vote-counting for metabolites
#' 
setClass("METAtables",
         slots = c(
           up = "matrix",
           down = "matrix",
           equal = "matrix",
           vote = "matrix"
         ))
setMethod(
  f = "initialize",
  signature = "METAtables",
  definition = function(.Object, Aup, Adown, Aequal, Avote) {
    .Object@up <- Aup
    .Object@down <- Adown
    .Object@equal <- Aequal
    .Object@vote <- Avote
    return(.Object)
  }
)



metmet <- function(datafile) {
  
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
  #' 
  #' @return METAtable S4 object with p-value combined, fold-change combined and vote-counting for each compound
  #' 
  #' @example 
  #' res.met <- metmet(datafile)
  
  vote <- NULL
  
  up <- matrix()
  down <- matrix()
  equal <- matrix()
  
  # Split data using trend variable
    
  trend <- levels(as.factor(datafile$trend)) 
    
  # One table per time
  for (t in 1:length(trend)) {
      
    met <- NULL
      
    # Subset of data
    data.t <- datafile[datafile$trend == trend[t], ]
      
    # Compound name as identifier
    unid <- unique(data.t$id) 
    metres <- data.frame()
      
    # Search data for each compound
    for (i in 1:length(unid)) {
       datacomp <- data.t[data.t$id == unid[i], ] 
        
      # Keep the identifier
      metres[i, "id"] <- unid[i]
        
      ## P-value
        
      logp <- log10(datacomp$pvalue) 
        
      # P-value weigthed for study number of individuals
      ponder <- logp * datacomp$N 
        
      # Combine p-values using Fisher's method      
      chi_sq <- (-2 / sum(datacomp$N)) * sum(ponder) 
        
      # Degrees of freedom
      degree <- 2 * nrow(datacomp) 
        
      # P-values comparison using Chi-squared distribution
      chisq <- pchisq(chi_sq, degree, lower.tail = F) 
        
      # Save result obtained
      metres[i, "pvalue combined"] <- chisq 
        
      ## Fold-change
        
      # Logarithmic transformation of fold-change
      logfc <- log2(datacomp$foldchange) 
        
      # Wheigthed mean for combining values
      mean_fc <- sum(logfc * datacomp$N) / sum(datacomp$N) 
        
      # Save result reversing logarithm
      metres[i, "foldchange combined"] <- 2 ^ mean_fc 
        
      # Save trend of result
      metres[i, "trend"] <- trend[t] 
        
      ## Extra information
        
      # Sum of total number of participants
      metres[i, "N total"] <- sum(datacomp$N)
        
      # Save references
      metres[i, "Reference"] <- paste(datacomp$ref, collapse = ";")
        
    }
      
      # Keep all results for compounds
      met <- rbind(met, metres)
      
      # Rename object using trend
      if (all(met$trend == "1" |
              tolower(met$trend) == "up")) {
        up <- as.matrix(met)
      } else if (all(met$trend == "-1" |
                     tolower(met$trend) == "down")) {
        down <- as.matrix(met)
      } else {
        equal <- as.matrix(met)
      }
      
    }
    
    ## Vote-counting
    
    # Unique compound as identifier using all trends
    comp <- unique(datafile$id)  
    
    # Create a data.frame
    vote <- data.frame(id = comp)
    
    for (c in 1:length(comp)) {
      
      # Subset data for each compound in all trends
      data.c <- datafile[datafile$id == comp[c],] 
      
      # Initialize vote count
      votec <- 0
      
      for (d in 1:nrow(data.c)) {
        
        # If is found up will sum 1
        if ("up" %in% tolower(data.c[d, "trend"]) |
            "1" %in% data.c[d, "trend"]) {
          votec <- votec + 1
          
          # If found down will sum -1
        } else if ("down" %in% tolower(data.c[d, "trend"]) |
                   "-1" %in% data.c[d, "trend"]) {
          votec <- votec - 1
        }
        
      # Equal trend not make any change  
      }
      
      # Save vote-counting obtained
      vote[vote$id == comp[c], "votes"] <- votec #guardem els vots
      
      # Save number of articles of each compound
      vote[vote$id == comp[c], "articles"] <- nrow(data.c)
      
      # Vote-counting/number of articles
      vote[vote$id == comp[c], "VC"] <- votec / nrow(data.c) 
    }
    
  # Save results in S4 object
  mets <- new("METAtables", up, down, equal, as.matrix(vote))
  
  return(mets)
}
