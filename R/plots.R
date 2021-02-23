
amanida_palette <- function() {
  
  #' Get nice colour-blind colours
  #' 
  #' @return vector of colours
  #' 
  
  c("gray67","#56B4E9", "#E69F00", "#009E73")
}


metaplot <- function(mets, cutoff = NULL) {
  
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
  
  col_palette <- amanida_palette()
  
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
  
 
  
  # Compounds with 2 or more reports
  
  cont <- as_tibble(mets@vote) %>% 
    mutate(articles = as.numeric(articles)) %>% 
    filter(articles >= 2)
  
  cont_ids <- cont %>% pull(id)
  
  # # Volcano plot
  
  # Scatter plot for logarithmic fold-change vs. -logarithmic p-value
  
  as_tibble(mets@stat) %>% 
    mutate( 
      # Format data needed
      across(c(pval,fc), as.numeric),
      # Negative logarithm of p-value for plot              
      lpval = -log10(pval),
      # Logarithm of fold-change
      lfc = log2(fc)) %>% 
    mutate(sig = case_when(
    (lpval > cut_pval & abs(lfc) < cut_fc) ~ "sig_p",
    (lfc < -cut_fc & lpval > cut_pval) ~ "sig_down",
    (lfc > cut_fc & lpval > cut_pval) ~ "sig_up",
    T ~ "no_sig"),
   label = case_when(
     sig == "sig_p" ~ "",
     sig == "no_sig" ~ "",
     T ~ id),
   reports = case_when(
     id %in% cont_ids ~ "y",
     T ~ "n" )) %>% 
    group_by(sig) %>% 
    {
      ggplot(., aes(lfc, lpval, label = label, colour = sig)) +
    geom_point() + 
    theme_minimal() +
    ggrepel::geom_text_repel(size = 3.5, 
                             fontface = "bold", 
                             segment.size = 0.4, 
                             point.padding = (unit(0.3, "lines")),
                             box.padding = unit(0.3, "lines"),
                             colour = "black") + 
    # Axis titles
    xlab( "log2(Fold-change)") + 
    ylab("-log10(p-value)") + 
    labs(colour = "") +
    
    # X axis breaks
    scale_x_continuous(breaks = seq(round(-max(abs(.$lfc)),0) - 1, 
                                    round(max(abs(.$lfc)),0) + 1, 1), 
                       limits = c(-max(abs(.$lfc)), max(abs(.$lfc)))) + 
    
    # Cutoff marks
    geom_hline(yintercept = cut_pval, 
               colour = "black", 
               linetype = "dashed") + 
    geom_vline(xintercept = c(cut_fc, -cut_fc),
               colour = "black", 
               linetype = "dashed") + 
    theme(legend.position = "none") +
    scale_color_manual(values = col_palette)
    }
    
}

# Plot articles number

voteplot <- function(mets) {
  
  #' Bar-plot for compouends number of reports
  #' 
  #' \code{voteplot} creates a bar-plot showing the number of entries for each compound. 
  #' 
  #' Independently of the compound trend, the total number of entries on each compound are plotted.
  #'  
  #' @param mets an S4 METAmet object obtained by \code{metmet}
  #'  
  #' @return bar-plot of results
  #' @example 
  #' voteplot(res.met)
  
  col_palette <- amanida_palette()
  
  # Subset vote-couting data
  as_tibble(mets@vote) %>% 
    mutate(
    articles = as.numeric(articles)) %>%
    ggplot(aes(reorder(id, -articles), articles)) + 
    geom_bar(stat = "identity", fill = col_palette[1] ) +
    theme(legend.position = "none") + 
    theme_classic() + 
    coord_flip() +
    xlab("id") + 
    ylab("Number of articles")
}

