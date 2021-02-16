
## Volcano Plot

metaplot <- function(mets, method, cutoff) {
  library(ggplot2)
  library(plotly)
  
  # Function: creates a volcano plot for data combined using metmet function
  # Arguments: S4 METAtables object, cutoff for statistical significance (optional)
  # Return: volcano plot
  
  # Search for cutoff argument
  if (hasArg(cutoff)) { 
    cuts <- cutoff
    if (length(cuts) != 2) {
      print("Two cut-off values needed")
    } else {
      cut_pval <- cuts[1]
      cut_fc <- cuts[2]
    }
    
    # If not cutoff argument convention values are established: 
  } else {
    # Alpha < 0.05 
    cut_pval <- 0.05
    
    # Log(fold-change) = 1.5
    
    cut_fc <- 2.83
  }
  
    # Join all data in S4 object
    nms <- slotNames(mets)
    
    metall <- NULL
    
    for (s in 1:3) {
      
      meta <- data.frame(slot(mets, nms[s]))
      
      if (length(meta) == 1) {
        next
      }
      
      metall <- rbind(metall, meta)
      
    }
    
    # Format data needed
    metall$id <- as.character(metall$id)
    
    # Negative logarithm of p-value for plot
    metall$ppval <- -log10(as.numeric(as.character(metall$pvalue.combined)))
    
    # Logarithm of fold-change
    metall$pfc <- log(as.numeric(as.character(metall$foldchange.combined)),2) 
    
    # Subset statistically significant values based on cutoff
    
    # Subset compounds with significant fold-change and p-value (down-regulated)
    hgdown <- metall[metall$pfc < -log(cut_fc, 2) & 
                       metall$ppval > -log10(cut_pval),]
    
    # Subset compounds with significant fold-change and p-value (up-regulated)
    hgup <- metall[metall$pfc > log(cut_fc, 2) & 
                     metall$ppval > -log10(cut_pval),]
    
    # Subset compounds with significant p-value and NOT fold-change
    hgmix <- metall[metall$ppval > -log10(cut_pval) & 
                      abs(metall$pfc) < log(cut_fc, 2),]
    
    # Vote-counting and articles information
    cont <- data.frame(slot(mets, nms[4]))
    
    # Number of articles found
    cont$articles <- as.numeric(as.character(cont$articles))
    
    # Votes obtained
    cont$votes <- as.numeric(as.character(cont$votes))
    
    # Filter compounds found in more than 2 articles
    cont_id <- cont[cont$articles >= 2,]
    
    hgcont <- metall[metall$id %in% cont_id$id,]
    
    # Range of logarithmic fold-change for plot aesthetics
    rx <- range(metall$pfc)
    
    # Volcano plot
    
    # Scatter plot for logarithmic fold-change vs. -logarithmic p-value
    d <- ggplot(metall, ggplot2::aes(pfc, ppval, label = id)) 
    
    e <- d  + ggplot2::theme_minimal() + 
      
      # Color-blind friendly colors
      geom_point(color = "gray67") + 
      
      # Colors for diferent statistical significant values
      geom_point(data = hgmix, aes(pfc, ppval), color = "#E69F00") + 
      geom_point(data = hgdown, aes(pfc, ppval), color = "#56B4E9") +
      geom_point(data = hgup, aes(pfc, ppval), color = "#009E73") +
      geom_point(data = hgcont, aes(pfc, ppval), shape = 8) +
      
      # Label compounds significant for fold-change and p-value
      ggrepel::geom_text_repel(label = 
                                 ifelse((metall$pfc > log(cut_fc, 2) | 
                                           metall$pfc < -log(cut_fc, 2)) &
                                                metall$ppval > -log10(cut_pval),
                                              metall$id, ''), 
                               size = 3.5, 
                               fontface = "bold", 
                               segment.size = 0.4, 
                               point.padding = (unit(0.3, "lines")),
                               box.padding = unit(0.3, "lines")) +
      # Axis titles
      xlab( "log2(Fold-change)") + ylab("-log10(p-value)") + labs(colour = "") +
      
      # X axis breaks
      scale_x_continuous(breaks = seq(round(-max(abs(rx)),0) - 1, 
                                      round(max(abs(rx)),0) + 1, 1), 
                         limits = c(-max(abs(rx)), max(abs(rx)))) + 
      
      # Cutoff marks
      geom_hline(yintercept = -log10(cut_pval), 
                 colour = "black", 
                 linetype = "dashed") + 
      geom_vline(xintercept = c(log(cut_fc, 2), -log(cut_fc, 2)),
                 colour = "black", 
                 linetype = "dashed") 
    
    # Show plot
    print(e)
}

# Plot articles number

voteplot <- function(mets) {
  library(ggplot2)
  library(plotly)
  
  # Function: create a bar-plot for the number of studies found on each compound
  # Arguments: S4 METAtables object
  # Return: Bar-plot
  
  # Subset vote-couting data
  data <- data.frame(mets@vote)
  
  # Sort data decreasing
  data <- data[order(data$articles),]
  
  # Bar-plot
  ggplot(data, aes(id, articles)) + 
    geom_bar(stat = "identity", fill = rgb(0.2,0.4,0.7,0.6) ) +
    theme(legend.position = "none") + theme_minimal() + coord_flip() +
    xlim(rev(levels(data$id))) 
}

