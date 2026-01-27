require(tidyverse)
require(ggfortify)
require(ggpubr)
require(ggrepel)
require(ggforce)
require(RColorBrewer)
require(patchwork)

PvA_corr.plot <- function(
    data, save = NULL, 
    x = "Actual_s_dss", y = "Mean_Combo_s_dss",
    xlab = "Observed Mean Combo DSS", 
    ylab = "Predicted Mean Combo DSS",
    xlim = c(-0.1,1), ylim = c(0,1), add = "reg.line",
    color = "cor", colorname = "Spearman\nCorr Coef",
    size = "n", sizename = "Sample size", text.size = 4.5,
    label1 = "Combination", label2 = "venaza",
    palette = c("green","black", "red"),
    ...
){
  
  if (color!="black"){
    midpoint = range(data[,color])
    midpoint = (0-midpoint[1])/(midpoint[2]-midpoint[1]) 
  } else {midpoint=0.5}
  p <- data %>%
    ggscatter(., x = x, y = y,
              add = add, 
              cor.coef.size = 4.9,
              cor.coef = TRUE, cor.method = "spearman",
              xlab = xlab, ylab = ylab,
              color = color, 
              size = size, alpha=0.7) + 
    scale_size_continuous(range = c(2,7), name = sizename) +
    scale_color_gradientn(colours = palette,
                          values = c(0, midpoint-0.1,midpoint+0.1, 1),
                          n.break = 5,
                          name = colorname) +
    xlim(xlim) + ylim(ylim)
  
  if (!is.na(label1)){
    p <- p + geom_text_repel(
      data = data, aes(x = get(x), y = get(y), label = get(label1)),
      box.padding = unit(1.2, "lines"), size=text.size,
      point.padding = unit(0.5, "lines"),
      segment.color = "black", min.segment.length = 0.1,
      force = 1, max.overlaps = 20, seed = 123)
    }

  if (!is.na(label2)){
    p <- p + geom_text_repel(
      data = data, aes(x = get(x), y = get(y), label = get(label2)),
      color = "blue", fontface="bold",
      box.padding = unit(4, "lines"), size=text.size,
      point.padding = unit(0.5, "lines"),
      segment.color = "black", min.segment.length = 0.03,
      force = 1, max.overlaps = 20, nudge_x = -0.4)
    }

  p <- p +
    theme_classic(base_size = 14)  +
    theme(legend.position = "right",
          axis.line = element_line(colour = "black"),
          legend.title = element_text(size=12),
          legend.text = element_text(size=12),
          aspect.ratio = 1
    )
  
  if (!is.null(save)){
    ggsave(save, plot=p, units="cm", ...)
  }
  return (p)
}


#' Plot top N drug frequency
plot_top_drugs <- function(data, top_n = 20, subtitle = "", title = "") {
  # title <- title %||% paste("Top", top_n, "Drug Partners")
  
  max_val <- ceiling(max(data$Proportion) / 5) * 5
  step <- ifelse(max_val >= 50, 10, 5)
  
  ggplot(data[1:top_n, ], aes(Proportion, reorder(Drug, Proportion))) +
    geom_col(aes(fill = Drug_Class), width = 0.8) +
    scale_x_continuous(
      limits = c(0, max_val),
      breaks = seq(0, max_val, by = step),
      expand = c(0, 0),
      position = "bottom"
    ) +
    scale_y_discrete(expand = expansion(add = c(0, 0.5))) +
    scale_fill_manual(values = drugclass, name = "Drug Class") +
    labs(x = "Proportion (%) of the Top Predicted Combinations", y = "Top Drug Partners", title = title, subtitle = subtitle) +
    theme_pubr() +
    theme(
      aspect.ratio = 0.8,
      legend.position = c(.75, .25),
      legend.text = element_text(size = 11),
      legend.key.height = unit(0.2, "cm"),
      panel.grid.major.x = element_line(color = "#A8BAC4", size = 0.3)
    )
}

