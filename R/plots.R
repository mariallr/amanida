
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

