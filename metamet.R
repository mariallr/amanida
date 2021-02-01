

## Load data

data.read <- function(file, separator, coln) {
  library(dplyr)
  csv <- grep("\\.csv", file) #mirem si es csv o excel per llegir el document
  excel <- grep("\\.xlsx|\\.xls", file)
  txt <- grep("\\.txt", file)
  if (length(csv) == 1) {
    datafile <- read.csv(file, sep = separator)
  } else if (length(excel) == 1) {
    datafile <- readxl::read_excel(file)
  } else if (length(txt) == 1){
    datafile <- read.table(file, sep = separator)
  } else {
    print("Not compatible format, try csv, excel or txt")
  }
  
  datafile <- dplyr::as_tibble(datafile)
  datafile <- datafile %>% select(all_of(coln)) #seleccionem només les columnes que tenen les dades que ens interesa. L'usuari ha de donar els noms de les columnes. 
  
  if (length(coln) == 6) { #hi ha l'opció d'afegir la direcció dels canvis o no
    colnames(datafile) <- c("id", "pvalue", "foldchange", "N", "trend", "ref") #canviem el nom de les columnes
  } else { 
    colnames(datafile) <- c("id", "pvalue", "foldchange", "N", "ref")
    if (length(datafile$foldchange[datafile$foldchange < 0])) { #si tenim valors positius i negatius podem extreure el trend
      datafile[datafile$foldchange < 0, "trend"] <- -1 #valors negatius posem -1
      datafile[datafile$foldchange > 0, "trend"] <- 1 #valors positius 1
      datafile[datafile$foldchange == 0, "trend"] <- 0 #si no hi ha canvi 0
    } 
  }
  
  
  return(datafile)

}


## Statistics

# S4 object to save the data
setClass("METAtables", 
         slots = c(up = "matrix", 
                   down = "matrix", 
                   equal = "matrix", 
                   vote = "matrix"))
setMethod( f = "initialize", signature = "METAtables",
           definition = function (.Object, Aup, Adown, Aequal, Avote) {
             .Object@up <- Aup
             .Object@down <- Adown
             .Object@equal <- Aequal
             .Object@vote <- Avote
             return(.Object)
           })



metmet <- function(datafile) {
  
  vote <- NULL
  
  up <- matrix()
  down <- matrix()
  equal <- matrix()
  
  if ("trend" %in% colnames(datafile)) { #primer mirem si s'ha facilitat la direcció dels canvis
    # split data in up and down
    trend <- levels(as.factor(datafile$trend)) #convertim en factor
    
    for (t in 1:length(trend)) { #dividim les dades en up i down
      met <- NULL
      
      data.t <- datafile[datafile$trend == trend[t],]
      
      unid <- unique(data.t$id) # per cada compost en les dades
      metres <- data.frame()
      for (i in 1:length(unid)) { 
        datacomp <- data.t[data.t$id == unid[i],] #agafem les dades per aquest compost
        metres[i,"id"] <- unid[i]
        
        #p-value
        logp <- log10(datacomp$pvalue) #calculem logaritme del p-value
        ponder <- logp * datacomp$N #ponderació logaritme(p-valor) * N de cada estudi
        
        chi_sq <- (-2/sum(datacomp$N)) * sum(ponder) # combinació de p-valors utilitizant el métode Fisher
        degree <- 2*nrow(datacomp) #graus de llibertat
        chisq <- pchisq(chi_sq, degree, lower.tail = F) #comparació dels p-valors amb una distribució chi-quadrat
        metres[i, "pvalue combined"] <- chisq #guardem els valors en una taula
        
        #fold-change
        
        logfc <- log2(datacomp$foldchange) #transformació logaritmica del fold-change
        mean_fc <- sum(logfc *datacomp$N) / sum(datacomp$N) #mitjana dels valors
        metres[i,"foldchange combined"] <- 2^mean_fc #treiem el logaritme
        metres[i, "trend"] <- trend[t] #indiquem si és up o down
        
        # extra info
        
        metres[i, "N total"] <- sum(datacomp$N)
        metres[i, "Reference"] <- paste(datacomp$ref, collapse = ";")
        
      }
      
      met <- rbind(met, metres)
      
      if (all(met$trend == "1" | tolower(met$trend) == "up" )) { #guardem la taula al nivell que correspongui
        up <- as.matrix(met)
      } else if (all(met$trend == "-1" | tolower(met$trend) == "down" )) {
        down <- as.matrix(met)
      } else {
        equal <- as.matrix(met)
      }
      
    }
    
    #vote-counting
    
    comp <- unique(datafile$id)  #calculem el vote-counting per cada compost
    
    vote <- data.frame(id = comp)
    
    for (c in 1:length(comp)) {
      data.c <- datafile[datafile$id == comp[c], ] #ara agafem les dades sense separa per up/down
      
      votec <- 0 
      
      for (d in 1:nrow(data.c)) {
        if ("up" %in% tolower(data.c[d, "trend"]) | "1" %in% data.c[d, "trend"]){ #si es up sumem 1
          votec <- votec + 1
        } else if ("down" %in% tolower(data.c[d, "trend"]) | "-1" %in% data.c[d, "trend"]) { #si es down sumem -1
          votec <- votec - 1
        }
      }
      
      vote[vote$id == comp[c], "votes"] <- votec #guardem els vots
      
      vote[vote$id == comp[c], "articles"] <- nrow(data.c)
      
      vote[vote$id == comp[c], "VC"] <- votec/nrow(data.c) #dividim els vots pel total de articles reportant el compost
      
    }
      
  } else { # en cas que no donin trend considerem tots en la mateixa direcció
    met <- NULL
    data.t <- datafile
    unid <- unique(data.t$id)
    metres <- data.frame()
    for (i in 1:length(unid)) {
      datacomp <- data.t[data.t$id == unid[i],]
      metres[i,"id"] <- unid[i]
      
      #p-value
      logp <- log10(datacomp$pvalue)
      ponder <- logp * datacomp$N
      
      chi_sq <- (-2/sum(datacomp$N)) * sum(ponder)
      degree <- 2*nrow(datacomp)
      chisq <- pchisq(chi_sq, degree, lower.tail = F)
      metres[i, "pvalue combined"] <- -log10(chisq)
      
      #fold-change
      
      logfc <- log2(datacomp$foldchange) #transformació logaritmica del fold-change
      mean_fc <- sum(logfc *datacomp$N) / sum(datacomp$N) #mitjana dels valors
      metres[i,"foldchange combined"] <- 2^mean_fc #treiem el logaritme
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


## Volcano Plot

metaplot <- function(mets, method, cutoff) {
  library(ggplot2)
  library(plotly)
  
  if (hasArg(cutoff)) { #busca l'argument de cutoff per guardar els valors
    cuts <- cutoff
    if (length(cuts) != 2) {
      print("2 cut-off needed")
    } else {
      cut_pval <- cuts[1]
      cut_fc <- cuts[2]
    }
  } else {
    cut_pval <- 0.05
    cut_fc <- 2.83
  }

  
  if (method == "alltogether") { #si volem el gràfic amb totes les dades de trend juntes
    #preparem dades per al plot
    
    nms <- slotNames(mets)
    
    metall <- NULL
    for (s in 1:3) {
      
      meta <- data.frame(slot(mets, nms[s]))
      
      if(length(meta) == 1) {
        next
      }
      
      metall <- rbind(metall, meta)
      
    }
    
    metall$id <- as.character(metall$id)
    
    metall$ppval <- -log10(as.numeric(as.character(metall$pvalue.combined))) #calculem el log del p-valor 
    
    metall$pfc <- log(as.numeric(as.character(metall$foldchange.combined)),2) #tornem a afegir el logaritme al fold-change
    
    if (all(metall$trend != "NA")){ #arreglem els noms per a que es vegi bé
      metall$trend <- as.factor(metall$trend)
      
      for (i in 1:length(levels(metall$trend))) {
        if ("up" %in% tolower(levels(metall$trend)[i]) | 1 %in% levels(metall$trend)[i] | "1" %in% levels(metall$trend)[i]) {
          levels(metall$trend)[i] <- "Up-regulated"
        } else if ("down" %in% tolower(levels(metall$trend)[i]) | -1 %in% levels(metall$trend)[i] | "-1" %in% levels(metall$trend)[i]) {
          levels(metall$trend)[i] <- "Down-regulated"
        } else {
          levels(metall$trend)[i] <- "equal"
        }
      }
      
    }
      
      # filtrem taula per veure els significatius
      
      hgdown <- metall[metall$pfc < -log(cut_fc, 2) & metall$ppval > -log10(cut_pval),]
      hgup <- metall[metall$pfc > log(cut_fc, 2) & metall$ppval > -log10(cut_pval),]
      hgmix <- metall[metall$ppval > -log10(cut_pval) & abs(metall$pfc) < log(cut_fc, 2),]
      
      cont <- data.frame(slot(mets, nms[4]))
      cont$articles <- as.numeric(as.character(cont$articles))
      cont$votes <- as.numeric(as.character(cont$votes))
      cont_id <- cont[cont$articles >= 2,]
      
      hgcont <- metall[metall$id %in% cont_id$id,]
      
      rx <- range(metall$pfc)
      
      a <- list(x = metall$pfc[metall$pfc > log(cut_fc, 2)| metall$pfc < -log(cut_fc, 2)],
                y = metall$ppval[metall$ppval > -log10(cut_pval)], 
                text = ifelse((metall$pfc > log(cut_fc, 2) | metall$pfc < -log(cut_fc, 2)) & metall$ppval > -log10(cut_pval), metall$id, ''),
                xref = "x",
                yref = "y", 
                showarrow = T,
                arrowhead = 7,
                arrowsize = .3, 
                ax = 20, 
                ay = -40)
      
      #fem el volcano plot
      
      d <- ggplot(metall, ggplot2::aes(pfc, ppval, label = id)) 
      e <- d  + ggplot2::theme_minimal() + 
        geom_point(color = "gray67") + #paletta color-blind-friendly
        geom_point(data = hgmix, aes(pfc, ppval), color = "#E69F00") + 
        geom_point(data = hgdown, aes(pfc, ppval), color = "#56B4E9") +
        geom_point(data = hgup, aes(pfc, ppval), color = "#009E73") +
        geom_point(data = hgcont, aes(pfc, ppval), shape = 8) +
        xlab( "log2(Fold-change)") + ylab("-log10(p-value)") + ggtitle("CRC and advanced adenomas + pre-surgery vs control + post-surgery")+
        scale_x_continuous(breaks = seq(round(-max(abs(rx)),0)-1, 
                                                 round(max(abs(rx)),0)+1, 1), limits = c(-max(abs(rx)), max(abs(rx)))) + labs (colour = "") +
        geom_hline(yintercept = -log10(cut_pval), colour = "black", linetype = "dashed") + geom_vline(xintercept = c(log(cut_fc, 2), -log(cut_fc, 2)), colour = "black", linetype = "dashed") +
        theme(legend.position = "bottom") 
      f <- ggplotly(e) %>% layout(annotations = a)
      print(f)
    
  } else if (method == "onebyone") { #si volem un plot per cada trend
    #preparem dades per al plot
    
    nms <- slotNames(mets)
    
    for (s in 1:3) {
      
      meta <- data.frame(slot(mets, nms[s]))
      
      if(length(meta) == 1) {
        next
      }
      
      meta$id <- as.character(meta$id)
      
      meta$ppval <- -log10(as.numeric(as.character(meta$pvalue.combined))) #calculem el log del p-valor 
      
      meta$pfc <- log(as.numeric(as.character(meta$foldchange.combined)),2) #tornem a afegir el logaritme al fold-change
      
      if (all(meta$trend != "NA")){ #arreglem els noms per a que es vegi bé
        meta$trend <- as.factor(meta$trend)
        
        for (i in 1:length(levels(meta$trend))) {
          if ("up" %in% tolower(levels(meta$trend)[i]) | 1 %in% levels(meta$trend)[i] | "1" %in% levels(meta$trend)[i]) {
            levels(meta$trend)[i] <- "Up-regulated"
          } else if ("down" %in% tolower(levels(meta$trend)[i]) | -1 %in% levels(meta$trend)[i] | "-1" %in% levels(meta$trend)[i]) {
            levels(meta$trend)[i] <- "Down-regulated"
          } else {
            levels(meta$trend)[i] <- "equal"
          }
        }
        
      }
        
        #fem el volcano plot
        
        # filtrem taula per veure els significatius
        
        hgdown <- meta[meta$pfc < -log(cut_fc, 2) & meta$ppval > -log10(cut_pval),]
        hgup <- meta[meta$pfc > log(cut_fc, 2) & meta$ppval > -log10(cut_pval),]
        hgmix <- meta[meta$ppval > -log10(cut_pval) & abs(meta$pfc) < log(cut_fc, 2),]
        
        cont <- data.frame(slot(mets, nms[4]))
        cont$articles <- as.numeric(as.character(cont$articles))
        cont$votes <- as.numeric(as.character(cont$votes))
        cont_id <- cont[cont$articles >= 2 & cont$articles == abs(cont$votes),]
        
        hgcont <- meta[meta$id %in% cont_id$id,]
        
        rx <- range(meta$pfc)
        
        #fem el volcano plot
        
       
        
        d <- ggplot(meta, ggplot2::aes(pfc, ppval, label = id)) 
        e <- d  + theme_minimal() + 
          geom_point(color = "gray67") +  #paletta color-blind-friendly
          geom_point(data = hgmix, aes(pfc, ppval), color = "#E69F00") + 
          geom_point(data = hgdown, aes(pfc, ppval), color = "#56B4E9") +
          geom_point(data = hgup, aes(pfc, ppval), color = "#009E73") +
          geom_point(data = hgcont, aes(pfc, ppval), shape = 8) +
          ggrepel::geom_text_repel(label = ifelse((meta$pfc > log(cut_fc, 2) | meta$pfc < -log(cut_fc, 2)) & meta$ppval > -log10(cut_pval), meta$id, ''), 
                                   size = 3.5, fontface = "bold", segment.size = 0.4, point.padding = (unit(0.3, "lines")),
                                   box.padding = unit(0.3, "lines")) + ggtitle(nms[s]) +
          xlab( "log2(Fold-change)") + ylab("-log10(p-value)") +
          scale_x_continuous(breaks = seq(round(-max(abs(rx)),0)-1, 
                                                   round(max(abs(rx)),0)+1, 0.5), limits = c(-max(abs(rx)), max(abs(rx)))) + labs (colour = "") +
          geom_hline(yintercept = -log10(cut_pval), colour = "black", linetype = "dashed") + geom_vline(xintercept = c(log(cut_fc, 2), -log(cut_fc, 2)), colour = "black", linetype = "dashed") +
          theme(legend.position = "bottom")  
        print(e)
    }
  }
}

#Plot articles nº

voteplot <- function(mets) {
  library(ggplot2)
  library(plotly)
  data <- data.frame(mets@vote)
  
  data <- data[order(data$articles),]
  
  ggplot(data, aes(id, articles)) + geom_bar(stat = "identity", fill=rgb(0.2,0.4,0.7,0.6) ) +
    theme(legend.position = "none") + theme_minimal() + coord_flip() +
    xlim(rev(levels(data$id))) #ara ordenat alfabéticament, millor per nº articles?
}



# Exemple 

filename <- "~/OneDrive - URV/metaanalysis_pack/dataset2.xlsx"

datafile <- data.read(filename, coln = c("Compound Name", "P-value", "Fold-change", "N total", "Behaviour",  "References"))

res.met <- metmet(datafile)

write.csv(res.met@down, "~/OneDrive - URV/Review celia/down.csv")

metaplot(res.met, method = "onebyone", cutoff = c(0.05, 4))

voteplot(res.met)

filename <- "~/OneDrive - URV/Review celia/data_prepost.xlsx"

datafile <- data.read(filename, coln = c("Common Name", "p-value", "fold-change", "N", "trend", "Reference"))

res.met <- metmet(datafile)

metaplot(res.met, method = "alltogether", cutoff = c(0.05, 4))

write.csv(res.met@up, "~/OneDrive - URV/Review celia/up_all.csv")

filename <- "~/OneDrive - URV/Review celia/volcano_all.xlsx"

datafile <- data.read(filename, coln = c("Compound Name", "P-value", "Fold-change", "N total", "Behaviour",  "References"))

res.met <- metmet(datafile)

metaplot(res.met, method = "alltogether", cutoff = c(0.05, 4))












geom_GeomTextRepel(label = ifelse((metall$pfc > log(cut_fc, 2) | metall$pfc < -log(cut_fc, 2)) & metall$ppval > -log10(cut_pval), metall$id, ''), 
                   size = 3.5, fontface = "bold", segment.size = 0.4, point.padding = (unit(0.3, "lines")),
                   box.padding = unit(0.3, "lines")) +
