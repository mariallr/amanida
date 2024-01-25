
amanida_palette <- function() {
  
  #' Get nice colour-blind colours
  #' 
  #' @return vector of colours
  #' 
  set.seed(123) 
  
  c("#F4A460", "#87CEEB", "#CD5C5C", "#A9A9A9", "#FFEFD5")
}

volcano_plot <- function(mets, cutoff = NULL) {
  
  #' Volcano plot of combined results 
  #' 
  #' \code{volcano_plot} returns a volcano plot of the combined results on each metabolite obtained by \code{compute_amanida} function
  #' 
  #' Results are presented as -log10 for p-value and log2 for fold-change. 
  #' Values over the cut off are labeled. If not cutoff is provided will be used alpha 0.05 for p-value and 1.5 for logarithmic fold-change.
  #'  
  #' @param mets an S4 METAtables object
  #' @param cutoff values for p-value and fold-change significance
  #'  
  #' @return plot of results
  #'  
  #' @examples 
  #' \dontrun{
  #'    data("sample_data")
  #'    
  #'    amanida_result <- compute_amanida(sample_data)
  #'    volcano_plot(amanida_result)
  #' }
  #'
  #' @import ggplot2
  #' @import dplyr
  #' @export
  #' 
  
  articles = NULL; pval = NULL; fc = NULL; lfc = NULL; lpval = NULL; . = NULL; 
  label = NULL; sig = NULL; reports = NULL;
  set.seed(123)
  
  col_palette <- amanida_palette()
  
  # Search for cutoff argument
  if (hasArg(cutoff)) { 
    cuts <- cutoff
    
    if (length(cuts) != 2) {
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
  
  message(paste("The cut-off used are "), (10^-cut_pval), " for p-value (", round(cut_pval,2), 
          " in log10 scale) and ", 2^cut_fc, 
          " for fold-change (", round(cut_fc,2), " in log2 scale).", sep = "")
  
  # Compounds with 2 or more reports
  cont <- as_tibble(mets@vote) |> 
    mutate(articles = as.numeric(articles)) |> 
    filter(articles >= 2)
  
  cont_ids <- cont |> pull(id)
  
  # Function for labels
  case_character_type <- function(lfc, lpval) {
    case_when(
      (lfc < -cut_fc & lpval > cut_pval) ~ paste("p-value < ", 10^-cut_pval, 
                                                 "& fold-change < ", -2^cut_fc),
      (lpval > cut_pval & abs(lfc) < cut_fc) ~ paste("p-value < ", 10^-cut_pval),
      (lfc > cut_fc & lpval > cut_pval) ~ paste("p-value <", 10^-cut_pval, 
                                                "& fold-change >", 2^cut_fc),
      T ~ "under cut-offs"
    )}
  
  ## Volcano plot
  
  # Scatter plot for logarithmic fold-change vs. -logarithmic p-value
  met <- as_tibble(mets@stat) |>
    mutate( 
      # Format data needed
      across(c(pval,fc), as.numeric),
      # Negative logarithm of p-value for plot              
      lpval = -log10(pval),
      # Logarithm of fold-change
      lfc = log2(fc)) |>
    mutate(sig = case_character_type(lfc, lpval),
    label = case_when(
     sig == paste("p-value < ", 10^-cut_pval) ~ "",
     sig == "under cut-offs" ~ "",
     T ~ id),
   reports = case_when(
     id %in% cont_ids ~ "> 1 report",
     T ~ "single report" )) |>
    dplyr::group_by('sig')
  
    ggplot(met, aes(lfc, lpval, label = label, colour = sig)) +
    geom_point(aes(shape = reports), size = 2.5) + 
    scale_shape_manual(values = c(8, 16), name = "") +
    theme_minimal() +
    ggrepel::geom_text_repel(size = 4,
                             fontface = "bold",
                             segment.size = 0.4,
                             point.padding = (unit(0.3, "lines")),
                             box.padding = unit(0.3, "lines"),
                             colour = "black",
                             max.overlaps = Inf) +
    # Axis titles
    xlab( "log2(Fold-change)") + 
    ylab(expression(paste("-log10(", italic(p), "-value)"))) + 
    labs(colour = "", tag = "Created with amanida") +
    
    # X axis breaks
    # scale_x_continuous(breaks = seq(round(-max(abs(lfc)),0) - 1, 
    #                                 round(max(abs(lfc)),0) + 1, 1), 
    #                    limits = c(-max(abs(lfc)), max(abs(lfc)))) + 
    
    # Cutoff marks
    geom_hline(yintercept = cut_pval, 
               colour = "black", 
               linetype = "dashed") + 
    geom_vline(xintercept = c(cut_fc, -cut_fc),
               colour = "black", 
               linetype = "dashed") + 
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5),
          legend.text = element_text(size = 10),
          plot.tag = element_text(size = 9, colour = "grey"),
          plot.tag.position = "bottomright") +
    guides(col = guide_legend(nrow = 2, byrow = T)) + 
    guides(shape = guide_legend(nrow = 2, byrow = T)) +
    scale_color_manual(values = col_palette) +
    ggtitle("Volcano plot of adapated meta-analysis results")
}


# Plot for vote-counting
vote_plot <- function(mets, counts = NULL) {
  
  #' Bar-plot for compounds vote-counting
  #' 
  #' \code{vote_plot} creates a bar-plot showing the vote-count for each compound. 
  #' 
  #' Vote-couting is the sum of number of reports up-regulated and the substraction of reports down-regulated. 
  #'  
  #' @param mets an S4 METAtables object obtained by \code{compute_amanida} or \code{amanida_vote}.
  #' @param counts value of vote-counting cut-off. Will be only displayed data over the cut-off.
  #'  
  #' @return a ggplot bar-plot showing the vote-count per compound
  #' 
  #' @importFrom stats reorder
  #' @examples 
  #' \dontrun{
  #'     data("sample_data")
  #'     result <- compute_amanida(sample_data)
  #'     vote_plot(result)
  #' }
  #' 
  #' @import ggplot2
  #' @import dplyr
  #' @export
  #' 
  
  votes = NULL; . = NULL;
  set.seed(123)
  
  col_palette <- amanida_palette()
  
  if (hasArg(counts)) { 
    cuts <- counts
    
    if (length(counts) != 1) {
      stop( "Please indicate one cut-off only")
    }
  } else {
    cuts <- 1
  }
  
  message("Cut-off for votes is ", cuts, ".", sep = "")
  
  # Subset vote-couting data
  tb <- as_tibble(mets@vote) |> 
    mutate(
    votes = as.numeric(votes),
    test = case_when(
      votes > 0 ~ 0,
      T ~ 1)) |>
    filter (abs(votes) >= cuts)
  
  if(nrow(tb) > 30) {
    message("Too much values, only showing 30 highest values. Please check counts parameter.")
    
   tb <- tb  |>
      slice_max(abs(votes), n = 30, with_ties = FALSE) 
  }
  
  min_p <- min(tb$votes)
  max_p <- max(tb$votes)
  
   ggplot(tb, aes(reorder(id, votes), votes, fill = votes)) + 
    geom_bar(stat = "identity", show.legend = F, width = .5
             ) +
    geom_text(aes(label = reorder(id, votes)), vjust = 0.2, size = 3.5, 
              hjust = tb$test) +
    scale_fill_gradient(low = col_palette[3], high = col_palette[5]) +
    theme_light() + 
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(size = 10),
          axis.title = element_text(size = 10),
          axis.ticks.y = element_blank(), 
          plot.title = element_text(size = 12), 
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(), 
          panel.border = element_blank(), 
          panel.grid.major.x = element_line(linetype = "dashed"),
          plot.tag = element_text(size = 9, colour = "grey"),
          plot.tag.position = "bottomright") +
    coord_flip() +
    ylab("Vote-counting") +
    xlab('')+
    labs(tag = "Created with amanida") +
    ggtitle("Total vote count of compounds behaviour") +
    scale_y_continuous(expand = c(0.6, 0), 
                       breaks = seq(min_p, max_p, by = 1),
                       limits = c(min_p - 1,
                                  max_p + 1)
                       )
  
}


# Piramid plot
explore_plot <- function(data, type = "all", counts = NULL) {
  
  #' Plot for compounds divergence in reports
  #' 
  #' \code{explore_plot} creates a bar-plot showing the votes divided in up-regulated and down-regulated and the global result for each compound. 
  #' 
  #' Sum of votes divided by trend are plotted, then is obtained the total result by compound summing both trends.  
  #'  
  #' @param data an tibble obtained by \code{amanida_read} w/o names checked by \code{check_names}
  #' @param type select the subset of data to plot. Options are: 
  #' \itemize{
  #'    \item "all": all data will be displayed
  #'    \item "sub": only data over counts value will be displayed. Need counts value.
  #'    \item "mix": will display data over count value and elements with reports in both trends.Need counts value.
  #'    }
  #' @param counts value of vote-counting cut-off. Will be only displayed data over the cut-off.  
  #'  
  #' @return a ggplot bar-plot showing the sum of votes for each compound divided by the trend
  #' 
  #' @importFrom stats reorder
  #' @import ggplot2
  #' @import dplyr
  #' 
  #' @examples 
  #' data("sample_data")
  #' 
  #' explore_plot(sample_data, type = "mix", counts = 1)
  #' 
  #' @export
  #' 

  trend = NULL; trend_l = NULL; N = NULL; vc = NULL; . = NULL; 
  cont = NULL; lab = NULL;
  set.seed(123)
  
  col_palette <- amanida_palette()
  
  if (hasArg(counts)) { 
    cuts <- counts
    
    if (length(cuts) != 1) {
      stop( "Please indicate one cut-off")
    }
  } else {
    stop("Function needs counts parameter")
  } 
  
  message("Cut-off for votes is ", cuts, ".", sep = "")
  
  if("cid" %in% colnames(data)){
    
    data$id <- tolower(unlist(data$id_mod))
  }
  
  if (type == "all") {
    dt <- data |>
      mutate(
        trend_l = case_when(
          trend == -1 ~ "Down-regulated", 
          T ~ "Up-regulated"
        )
      ) |> dplyr::group_by(id) |> 
      mutate(vc = sum(trend)) |>
      dplyr::group_by(id, trend_l) |>
      summarise(
        cont = n(),
        total_N = sum(N),
        vc = unique(vc),
        lab = c("Vote-counting")
      ) |>
      mutate(cont = case_when(
        trend_l == "Down-regulated" ~ cont*-1,
        T ~ cont*1
      ))
    
  } else if (type == "sub") {
    dt <- data |>
      mutate(
        trend_l = case_when(
          trend == -1 ~ "Down-regulated", 
          T ~ "Up-regulated"
        )
      ) |> dplyr::group_by(id) |>
      mutate(vc = sum(trend)) |>
      dplyr::group_by(id, trend_l) |>
      summarise(
        cont = n(),
        total_N = sum(N),
        vc = unique(vc),
        lab = c("Vote-counting")
      ) |>
      mutate(cont = case_when(
        trend_l == "Down-regulated" ~ cont*-1,
        T ~ cont*1
      )) |>
      filter(vc > cuts | vc < -1*cuts)
    
  } else if (type == "mix") {
    dt <- data |>
      mutate(
        trend_l = case_when(
          trend == -1 ~ "Down-regulated", 
          T ~ "Up-regulated"
        )
      ) |> dplyr::group_by(id) |> 
      mutate(vc = sum(trend)) |>
      dplyr::group_by(id, trend_l) |>
      summarise(
        cont = n(),
        total_N = sum(N),
        vc = unique(vc),
        lab = c("Vote-counting")
      ) |>
      mutate(cont = case_when(
        trend_l == "Down-regulated" ~ cont*-1,
        T ~ cont*1
      )) |>
      filter(vc > cuts | vc < -1*cuts |
               vc != cont)
  }
  
  # Prepare data for plot
  
  if(nrow(dt) > 25) {
    message("Too much values, only showing 30 first values. Please check counts and/or type parameters.")

    dt <- dt |>
      ungroup() |>
      arrange(id) |>
      slice(1:30)
  }
  
  max_p <- max(abs(dt$cont))
  
      ggplot(dt, mapping = aes(x = cont,
                              y = reorder(id, vc), fill = trend_l)) +
        geom_bar(aes(x = ifelse(test = trend_l == "Up-regulated",
                                yes = cont, no = cont)), stat = "identity",
                 alpha = .3) +
        scale_fill_manual(values = col_palette[2:3]) +
        geom_segment(aes(x = 0, xend = vc, 
                         y = id, yend = id, linetype = lab),
                     linewidth = 0.4, alpha = 0.9, 
                     arrow = arrow(length = unit(0.1, "cm")), 
                     lineend = "round", linejoin = "round") +
        scale_x_continuous(labels = abs, limits = max_p * c(-1,1) * 1.01) +
        scale_color_manual(values = c(col_palette[2], col_palette[3])) +
        theme_minimal() +
        xlab("Counts by trend") + 
        ylab("") +
        labs(fill = "Counts by trend", tag = "Created with amanida") +
        ggtitle("Qualitative compounds trend plot") +
        theme(legend.position = "bottom", legend.title = element_blank(),
              axis.text.y = element_text(size = 14),
              plot.tag = element_text(size = 9, colour = "grey"),
              plot.tag.position = "bottomright") +
        guides(col = guide_legend(nrow = 2, byrow = T)) + 
        guides(shape = guide_legend(nrow = 2, byrow = T)) 
}


