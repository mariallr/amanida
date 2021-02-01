

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
  vote <- NULL
  
  up <- matrix()
  down <- matrix()
  equal <- matrix()
  
  if ("trend" %in% colnames(datafile)) {
    #primer mirem si s'ha facilitat la direcció dels canvis
    # split data in up and down
    trend <- levels(as.factor(datafile$trend)) #convertim en factor
    
    for (t in 1:length(trend)) {
      #dividim les dades en up i down
      met <- NULL
      
      data.t <- datafile[datafile$trend == trend[t], ]
      
      unid <- unique(data.t$id) # per cada compost en les dades
      metres <- data.frame()
      for (i in 1:length(unid)) {
        datacomp <-
          data.t[data.t$id == unid[i], ] #agafem les dades per aquest compost
        metres[i, "id"] <- unid[i]
        
        #p-value
        logp <-
          log10(datacomp$pvalue) #calculem logaritme del p-value
        ponder <-
          logp * datacomp$N #ponderació logaritme(p-valor) * N de cada estudi
        
        chi_sq <- (-2 / sum(datacomp$N)) * sum(ponder) # combinació de p-valors utilitizant el métode Fisher
        degree <- 2 * nrow(datacomp) #graus de llibertat
        chisq <-
          pchisq(chi_sq, degree, lower.tail = F) #comparació dels p-valors amb una distribució chi-quadrat
        metres[i, "pvalue combined"] <-
          chisq #guardem els valors en una taula
        
        #fold-change
        
        logfc <-
          log2(datacomp$foldchange) #transformació logaritmica del fold-change
        mean_fc <-
          sum(logfc * datacomp$N) / sum(datacomp$N) #mitjana dels valors
        metres[i, "foldchange combined"] <-
          2 ^ mean_fc #treiem el logaritme
        metres[i, "trend"] <- trend[t] #indiquem si és up o down
        
        # extra info
        
        metres[i, "N total"] <- sum(datacomp$N)
        metres[i, "Reference"] <-
          paste(datacomp$ref, collapse = ";")
        
      }
      
      met <- rbind(met, metres)
      
      if (all(met$trend == "1" |
              tolower(met$trend) == "up")) {
        #guardem la taula al nivell que correspongui
        up <- as.matrix(met)
      } else if (all(met$trend == "-1" |
                     tolower(met$trend) == "down")) {
        down <- as.matrix(met)
      } else {
        equal <- as.matrix(met)
      }
      
    }
    
    #vote-counting
    
    comp <-
      unique(datafile$id)  #calculem el vote-counting per cada compost
    
    vote <- data.frame(id = comp)
    
    for (c in 1:length(comp)) {
      data.c <-
        datafile[datafile$id == comp[c],] #ara agafem les dades sense separa per up/down
      
      votec <- 0
      
      for (d in 1:nrow(data.c)) {
        if ("up" %in% tolower(data.c[d, "trend"]) |
            "1" %in% data.c[d, "trend"]) {
          #si es up sumem 1
          votec <- votec + 1
        } else if ("down" %in% tolower(data.c[d, "trend"]) |
                   "-1" %in% data.c[d, "trend"]) {
          #si es down sumem -1
          votec <- votec - 1
        }
      }
      
      vote[vote$id == comp[c], "votes"] <- votec #guardem els vots
      
      vote[vote$id == comp[c], "articles"] <- nrow(data.c)
      
      vote[vote$id == comp[c], "VC"] <-
        votec / nrow(data.c) #dividim els vots pel total de articles reportant el compost
      
    }
    
  } else {
    # en cas que no donin trend considerem tots en la mateixa direcció
    met <- NULL
    data.t <- datafile
    unid <- unique(data.t$id)
    metres <- data.frame()
    for (i in 1:length(unid)) {
      datacomp <- data.t[data.t$id == unid[i], ]
      metres[i, "id"] <- unid[i]
      
      #p-value
      logp <- log10(datacomp$pvalue)
      ponder <- logp * datacomp$N
      
      chi_sq <- (-2 / sum(datacomp$N)) * sum(ponder)
      degree <- 2 * nrow(datacomp)
      chisq <- pchisq(chi_sq, degree, lower.tail = F)
      metres[i, "pvalue combined"] <- -log10(chisq)
      
      #fold-change
      
      logfc <-
        log2(datacomp$foldchange) #transformació logaritmica del fold-change
      mean_fc <-
        sum(logfc * datacomp$N) / sum(datacomp$N) #mitjana dels valors
      metres[i, "foldchange combined"] <-
        2 ^ mean_fc #treiem el logaritme
      metres[i, "trend"] <- trend[t] #indiquem si és up o down
      
      # extra info
      
      metres[i, "N total"] <- sum(datacomp$N)
      metres[i, "Reference"] <- paste(datacomp$ref, collapse = ";")
      
    }
    met <- rbind(met, metres)
    
    equal <- as.matrix(met)
  }
  
  
  mets <- new("METAtables", up, down, equal, as.matrix(vote))
  return(mets)
}
