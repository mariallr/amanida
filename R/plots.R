
## Volcano Plot

metaplot <- function(mets, cutoff = NULL) {
  library(ggplot2)
  
  #' Volcano plot of combined results 
  #' 
  #' \code{metaplot} returns a volcano plot of the combined results on each metabolite obtained by metmet function
  #' 
  #' Results are presented as -log10 for p-value and log2 for fold-change. 
  #' Values over the cut off are labeled. If not cutoff is provided will be used alpha 0.05 for p-value and 1.5 for logarithmic fold-change.
  #'  
  #'  @param mets an S4 METAtable object
  #'  @param cutoff values for p-value and fold-change significance
  #'  
  #'  @return plot of results
  #'  
  #'  @example 
  #'  metaplot(res.met, c(0.01, 2))
  #'  
  
  # Search for cutoff argument
  if (hasArg(cutoff)) { 
    cuts <- cutoff
    
    if(length(cuts) != 2) {
      stop( "Please indicate one cut-off for p-value and one for fold-change")
    }
    cut_pval <- -log10(cuts[1])
    cut_fc <- log2(cuts[2])
      
    # If not cutoff argument convention values are established: 
  } else {
    # Alpha < 0.05 
    cut_pval <- -log10(0.05)
    
    # Log(fold-change) = 1.5
    
    cut_fc <- log2(2.83)
  }
  
  
  data <- as_tibble(mets@stat)
  
  data <- data %>% mutate(pval,
                  # Format data needed
                  pval = as.numeric(pval), 
                  fc = as.numeric(fc),
                  # Negative logarithm of p-value for plot              
                  lpval = -log10(pval),
                  # Logarithm of fold-change
                  lfc = log2(fc)
                  )
  
      
    # Compounds with 2 or more reports
    
    cont <- as_tibble(mets@vote)
    
    cont <- cont %>% mutate(articles,
                            articles = as.numeric(articles)) %>% 
      filter(articles >= 2)
    
    # Volcano plot
    
    # Scatter plot for logarithmic fold-change vs. -logarithmic p-value
    d <- ggplot(data, ggplot2::aes(lfc, lpval, label = id)) 
    
    e <- d  + ggplot2::theme_minimal() + 
      
      # Color-blind friendly colors
      geom_point(color = "gray67") +
      
      # Colors for diferent statistical significant values
      geom_point(data = data %>% 
                   filter(lpval > cut_pval &
                            abs(lfc) < cut_fc), 
                 aes(lfc, lpval), color = "#E69F00") + 
      geom_point(data = data %>% filter(lfc < -cut_fc & 
                                          lpval > cut_pval), 
                 aes(lfc, lpval), color = "#56B4E9") +
      geom_point(data = data %>% filter(lfc > cut_fc & 
                                          lpval > cut_pval),
                 aes(lfc, lpval), color = "#009E73") +
      geom_point(data = data %>% filter(id %in% cont$id), 
                 aes(lfc, lpval), shape = 8) +
      
      # Label compounds significant for fold-change and p-value
      ggrepel::geom_text_repel(label = 
                                  data %>% mutate(id = case_when(
                                   (abs(lfc) < cut_fc & lpval < cut_pval) ~ "",
                                   (abs(lfc) > cut_fc & lpval < cut_pval) ~ "",
                                   (abs(lfc) < cut_fc & lpval > cut_pval) ~ "",
                                   T ~ id)) %>% .$id,
                               size = 3.5, 
                               fontface = "bold", 
                               segment.size = 0.4, 
                               point.padding = (unit(0.3, "lines")),
                               box.padding = unit(0.3, "lines")) + 
      # Axis titles
      xlab( "log2(Fold-change)") + ylab("-log10(p-value)") + labs(colour = "") +
      
      # X axis breaks
      scale_x_continuous(breaks = seq(round(-max(abs(data$lfc)),0) - 1, 
                                      round(max(abs(data$lfc)),0) + 1, 1), 
                         limits = c(-max(abs(data$lfc)), max(abs(data$lfc)))) + 
      
      # Cutoff marks
      geom_hline(yintercept = cut_pval, 
                 colour = "black", 
                 linetype = "dashed") + 
      geom_vline(xintercept = c(cut_fc, -cut_fc),
                 colour = "black", 
                 linetype = "dashed") 
    
    # Show plot
    print(e)
}

# Plot articles number

voteplot <- function(mets) {
  library(ggplot2)
  
  #' Bar-plot for compouends number of reports
  #' 
  #' \code{voteplot} creates a bar-plot showing the number of entries for each compound. 
  #' 
  #' Independently of the compound trend, the total number of entries on each compound are plotted.
  #'  
  #'  @params mets an S4 METAmet object obtained by \code{metmet}
  #'  
  #'  @return bar-plot of results
  #'  @example 
  #'  voteplot(res.met)
  
  # Subset vote-couting data
  data <- as_tibble(mets@vote)
  
  # Sort data decreasing
  data <- data %>% arrange(desc(articles))
  
  # Bar-plot
  ggplot(data, aes(id, articles)) + 
    geom_bar(stat = "identity", fill = rgb(0.2,0.4,0.7,0.6) ) +
    theme(legend.position = "none") + theme_classic() + coord_flip() +
    xlim(data$id) + xlab("id") + ylab("Nº of articles")
}

