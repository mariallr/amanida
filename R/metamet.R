

source("R/data.R")
source("R/plots.R")

## Statistics

# S4 object to save the data
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
  
  # Function: Combines a) p-values using Fisher methods weighted by N  and 
  #           chi-squared distribution
  #           b) Fold-change by weighted mean 
  #           Calculates vote-counting for each compound
  # Arguments: data obtained with data.read function
  # Returns: data-frame with p-value combined, fold-change combined and
  #          vote-counting for each compound
  
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
