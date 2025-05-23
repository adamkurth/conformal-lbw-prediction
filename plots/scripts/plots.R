# rm(list=ls())
library(rpart)
library(rpart.plot)
year <- 2021

# Determine type based on data path
use_without_2_5kg <- FALSE # Set to TRUE for type 2 (without 2.5kg), FALSE for type 1
type <- ifelse(use_without_2_5kg, 2, 1)

# Set appropriate data path based on type
if(type == 1) {
  data.pwd <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/birthweight_data/rebin/"
  plots.pwd <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/plots/type1"
  results.pwd <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/results/results1"
  load(sprintf("%s/quantile_cutpoints_%d.RData", data.pwd, 2020))
} else {
  data.pwd <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/birthweight_data/rebin_without_2.5kg/"
  plots.pwd <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/plots/type2"
  results.pwd <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/results/results2"
  load(sprintf("%s/quantile_cutpoints_%d.RData", data.pwd, 2020))
}

# data.pwd <- sprintf("/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/birthweight_data/rebin/")
# load(sprintf("%s/birthweight_data_rebin_%d.RData", data.pwd, year))
# load(sprintf("%s/informed_prior_%d.RData", data.pwd, 2020))
# load(sprintf("%s/bootstrap_results_%d.RData", data.pwd, year))
load(sprintf("%s/bootstrap_tree_results_%d.RData", results.pwd, year))
# load(sprintf("%s/dm_tree_%d.RData", data.pwd, year))


pdf(sprintf("%s/dm_tree_%d.pdf", plots.pwd, type), width = 10, height = 8)
plot(dm.tree, main="Dirichlet−Multinomial Decision Tree", uniform=TRUE, margin=0.1)
text(dm.tree, use.n=TRUE, cex=0.8, all=TRUE, pretty=0, xpd=NA, digits=3, font=3)
dev.off()



# K <- length(alphavec)
# barplot(alphavec,
#         names.arg = paste("cat", 1:K),
#         col = "skyblue",
#         border = "blue",
#         main = "Informed Dirichlet Prior (alphavec)",
#         xlab = "Birth Weight Categories",
#         ylab = "Prior Probability",
#         ylim = c(0, max(alphavec) * 1.2))
# 
# plots.pwd <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/plots"
# pdf(sprintf("%s/alphavec_plot_%d.pdf", plots.pwd, 2020), width = 10, height = 8)
# barplot(alphavec,
#         names.arg = paste0("C", 1:length(alphavec)),
#         col = "skyblue",
#         border = "blue",
#         main = "Informed Dirichlet Prior",
#         xlab = "Birth Weight Categories",
#         ylab = "Prior Probability",
#         ylim = c(0, max(alphavec) * 1.2))
# dev.off()



library(ggplot2)
library(gridExtra)
K <- length(alphavec)
df <- data.frame(
  Category = factor(paste0("C", 1:length(alphavec)), levels = paste0("C", 1:length(alphavec))), 
  Value = alphavec
 )

ggplot(df, aes(x = Category, y = Value)) +
  # Add bars with defined colors
  geom_col(fill = "skyblue", color = "blue", width = 0.7) +
  
  # Add labels and title with better typography
  labs(
    title = "Informed Dirichlet Prior",
    x = "Birth Weight Categories",
    y = "Prior Probability",
    caption = "Note: Values represent prior probabilities for each birth weight category"
  ) +
  
  # Set y-axis limits
  scale_y_continuous(
    limits = c(0, max(alphavec) * 1.2),
    expand = expansion(mult = c(0, 0.05))  # Remove padding at bottom, add a bit at top
  ) +
  
  # Apply a clean, publication-ready theme
  theme_minimal() +
  theme(
    # Title formatting
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    
    # Axis titles
    axis.title = element_text(size = 12, face = "bold"),
    
    # Axis text
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 0, vjust = 0.5),
    
    # Grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    # Add a subtle border
    panel.border = element_rect(color = "gray90", fill = NA, linewidth = 0.5),
    
    # Caption formatting
    plot.caption = element_text(size = 9, color = "gray30", hjust = 0),
    
    # Overall margins
    plot.margin = margin(20, 20, 20, 20)
  )

# To save the plot:
ggsave(sprintf("%s/dirichlet_prior_%d.pdf", boot.pwd, year), 
       width = 8, height = 6, units = "in", dpi = 300)


library(ggplot2)
library(patchwork)
library(dplyr)
library(grid)        # textGrob()

load(sprintf("%s/quantile_cutpoints_%d.RData", data.pwd, 2020))

cp   <- as.numeric(cut.points)
cats <- paste0("C", 1:K)

lb   <- cp[seq_len(min(K, length(cp) - 1))]
ub   <- cp[seq_len(min(K, length(cp) - 1)) + 1]

cat_labels <- paste0(cats[seq_along(lb)], ": ", lb, "–", ub, " g")
if (K > length(cp) - 1)           # optional extra bin
  cat_labels <- c(cat_labels,
                  paste0(cats[K], ": > ", cp[length(cp)], " g"))

## ------------------------------------------------------------------------
## 1.   build plotting data frames
## ------------------------------------------------------------------------
high_risk_df <- data.frame(
  Category    = factor(cats, levels = cats),
  Probability = high.risk.probs,
  Lower       = high.risk.lwr,
  Upper       = high.risk.upr
)
low_risk_df  <- data.frame(
  Category    = factor(cats, levels = cats),
  Probability = low.risk.probs,
  Lower       = low.risk.lwr,
  Upper       = low.risk.upr
)
combined_df  <- bind_rows(high_risk_df, low_risk_df)

## ── 2. helper styling ────────────────────────────────────────────────────
red_bar   <- "#FF6666"; red_border  <- "#990000"
blue_bar  <- "#377EB8"; blue_border <- "#00008B"

base_theme <- theme_minimal() +
  theme(plot.title   = element_text(size = 14, face = "bold", hjust = .5),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.text    = element_text(size = 10),
        plot.margin  = margin(5, 10, 5, 10))

## reset any half-closed graphic devices that can trigger the error
while (!is.null(dev.list())) dev.off()

## ------------------------------------------------------------------------
## 2.   global theme & bar colours
## ------------------------------------------------------------------------
red_bar  <- "#FF6666"; red_border  <- "#990000"
blue_bar <- "#377EB8"; blue_border <- "#00008B"

base_theme <- theme_minimal() +
  theme(axis.title.y = element_text(size = 12, face = "bold"),
        axis.text    = element_text(size = 10),
        plot.title   = element_text(size = 14, face = "bold", hjust = .5),
        plot.margin  = margin(5, 10, 5, 10))

## ------------------------------------------------------------------------
## 3.   helper to make a panel
## ------------------------------------------------------------------------
y_max <- max(combined_df$Upper, na.rm = TRUE) * 1.1   # same y-axis for both

make_panel <- function(df, fill, border, title) {
  ggplot(df, aes(Category, Probability)) +
    geom_col(fill = fill, colour = border, width = .7) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper),
                  width = .2, colour = "black", linewidth = .8) +
    scale_y_continuous(limits = c(0, y_max),
                       expand  = expansion(mult = c(0, .05))) +
    labs(title = title, x = "", y = "Probability") +
    base_theme
}

if (type == 1) {# ── TYPE 1 ──────────────────────────
  p1 <- make_panel(high_risk_df, red_bar,  red_border,  "High-Risk Subgroup")
  p2 <- make_panel(low_risk_df,  blue_bar, blue_border, "Low-Risk Subgroup")
  
} else {# ── TYPE 2 ──────────────────────────
  # (identical layout but y-limit based on *all* Upper values)
  p1 <- make_panel(high_risk_df, red_bar,  red_border,  "High-Risk Subgroup")
  p2 <- make_panel(low_risk_df,  blue_bar, blue_border, "Low-Risk Subgroup")
}

## 4. build a tidy legend grob --------------------------------------------
legend_grob <- grobTree(
  rectGrob(            # light background box
    x = 0, y = 1, width = 1, height = 1, hjust = 0, vjust = 1,
    gp = gpar(fill = "white", col = NA, alpha = 0.8)
  ),
  textGrob(
    paste(cat_labels, collapse = "\n"),
    x = 0.02, y = 0.98, hjust = 0, vjust = 1,
    gp = gpar(cex = 0.72, fontfamily = "sans")
  )
)

# turn the grob into its own tiny plot (1 × 1 “panel”)
legend_plot <- ggplot() + theme_void() +
  annotation_custom(legend_grob, xmin = -Inf, xmax = Inf,
                    ymin = -Inf, ymax = Inf)

## 5. combine panels + legend as a 2-column layout -------------------------
combined_plot <- (p1 / p2) | legend_plot   # patchwork “|” puts it in a column
combined_plot <- combined_plot +
  plot_layout(widths = c(1, 0.25)) +      # 25 % width for the legend column
  plot_annotation(
    title = "Birth-weight Category Bootstrap Probabilities",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = .5),
      plot.margin = margin(10, 10, 10, 10)
    )
  )

## 6. export / view --------------------------------------------------------
ggsave(sprintf("%s/high_low_risk_pred_%d.pdf", boot.pwd, year),
       combined_plot, width = 9, height = 8, units = "in", dpi = 300)

combined_plot   # show in Viewer

## ── 5. table for the console ─────────────────────────────────────────────
cat("\nMean probabilities for all categories:\n")
results.data <- data.frame(cat = cats,
                 high.risk.prob = high.risk.probs,
                 high.risk.lwr  = high.risk.lwr,
                 high.risk.upr  = high.risk.upr,
                 low.risk.prob  = low.risk.probs,
                 low.risk.lwr   = low.risk.lwr,
                 low.risk.upr   = low.risk.upr)





library(ggplot2)
library(patchwork)
library(dplyr)
library(grid)

## ── 0. load cut-points & build labels ───────────────────────────────────
load(sprintf("%s/quantile_cutpoints_%d.RData", data.pwd, 2020))

cp   <- as.numeric(cut.points)
K    <- length(cp)                      # total possible bins
plot_idx <- seq_len(min(10, K))         # -> C1 … C10 only

cats <- paste0("C", 1:K)[plot_idx]      # keep first 10 labels

lb <- cp[plot_idx]
ub <- cp[plot_idx + 1]

cat_labels <- paste0(cats, ": ", lb, "–", ub, " g")

## ── 1. build trimmed plotting data frames ───────────────────────────────
high_risk_df <- data.frame(
  Category    = factor(cats, levels = cats),
  Probability = high.risk.probs[plot_idx],
  Lower       = high.risk.lwr [plot_idx],
  Upper       = high.risk.upr [plot_idx]
)
low_risk_df <- data.frame(
  Category    = factor(cats, levels = cats),
  Probability = low.risk.probs[plot_idx],
  Lower       = low.risk.lwr [plot_idx],
  Upper       = low.risk.upr [plot_idx]
)
combined_df <- bind_rows(high_risk_df, low_risk_df)

## ── 2. palettes & theme (unchanged) ─────────────────────────────────────
red_bar <- "#FF6666"; red_border <- "#990000"
blue_bar <- "#377EB8"; blue_border <- "#00008B"

base_theme <- theme_minimal() +
  theme(axis.title.y = element_text(size = 12, face = "bold"),
        axis.text    = element_text(size = 10),
        plot.title   = element_text(size = 14, face = "bold", hjust = .5),
        plot.margin  = margin(5, 10, 5, 10))

## ── 3. helper to make a panel ───────────────────────────────────────────
y_max <- max(combined_df$Upper, na.rm = TRUE) * 1.1
make_panel <- function(df, fill, border, title) {
  ggplot(df, aes(Category, Probability)) +
    geom_col(fill = fill, colour = border, width = .7) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper),
                  width = .2, colour = "black", linewidth = .8) +
    scale_y_continuous(limits = c(0, y_max),
                       expand = expansion(mult = c(0, .05))) +
    labs(title = title, x = "", y = "Probability") +
    base_theme
}
p1 <- make_panel(high_risk_df, red_bar,  red_border,  "High-Risk Subgroup")
p2 <- make_panel(low_risk_df,  blue_bar, blue_border, "Low-Risk Subgroup")

## ── 4. legend grob (only 10 labels) ─────────────────────────────────────
legend_grob <- grobTree(
  rectGrob(gp = gpar(fill = "white", col = NA, alpha = 0.9)),
  textGrob(paste(cat_labels, collapse = "\n"),
           x = 0.02, y = 0.98, hjust = 0, vjust = 1,
           gp = gpar(cex = 0.72))
)
legend_plot <- wrap_elements(full = legend_grob) +
  theme(plot.margin = margin(5, 5, 5, 5))

## small spacer so legend hugs the top-right
legend_col <- plot_spacer() / legend_plot +
  plot_layout(heights = c(0.05, 0.95))

## ── 5. assemble & save ──────────────────────────────────────────────────

combined_plot <- (p1 / p2) | legend_plot   # patchwork “|” puts it in a column
combined_plot <- combined_plot +
  plot_layout(widths = c(1, 0.25)) +      # 25 % width for the legend column
  plot_annotation(
    title = "Birth-weight Category Bootstrap Probabilities",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = .5),
      plot.margin = margin(10, 10, 10, 10)
    )
  )

ggsave(sprintf("%s/high_low_risk_pred_%d.pdf", boot.pwd, year),
       combined_plot, width = 9, height = 8, units = "in", dpi = 300)

combined_plot



# result.table -> latex table
library(xtable)







# pdf(sprintf("%s/dm_tree_%d.pdf", plots.pwd, type), width = 10, height = 8)
# plot(dm.tree, uniform=TRUE, margin=0.1)
# text(dm.tree, use.n=TRUE, cex=0.5)
# dev.off()




# for depth comparisons

# Create a higher resolution PDF with larger text and spacing
pdf(sprintf("%s/decision_tree_comparison_%d.pdf", results.pwd, year), 
    width = 16, height = 10, pointsize = 14)

# Set up a 2x2 layout with more generous margins
par(mfrow = c(2, 2), mar = c(1, 1, 3, 1), oma = c(1, 1, 1, 1))

# Depth 2 - Maximized for clarity
plot(trees$`2`, 
     main = "Depth = 2", 
     uniform = TRUE, 
     margin = 0.2,       # Increased margin for better spacing
     cex = 2.2)          # Very large node text for readability
text(trees$`2`, 
     use.n = FALSE,       # Show number of observations
     all = TRUE,         # Label all nodes 
     fancy = FALSE,      # Simple labels
     cex = 1.8)          # Very large text

# Depth 3
plot(trees$`3`, 
     main = "Depth = 3", 
     uniform = TRUE, 
     margin = 0.2, 
     cex = 2.0)
text(trees$`3`, 
     use.n = FALSE, 
     all = TRUE, 
     fancy = FALSE, 
     cex = 1.6)

# Depth 4
plot(trees$`4`, 
     main = "Depth = 4", 
     uniform = TRUE, 
     margin = 0.2, 
     cex = 1.8)
text(trees$`4`, 
     use.n = FALSE, 
     all = TRUE, 
     fancy = FALSE, 
     cex = 1.4)

# Depth 5
plot(trees$`5`, 
     main = "Depth = 5", 
     uniform = TRUE, 
     margin = 0.2, 
     
     cex = 1.6)
text(trees$`5`, 
     use.n = FALSE, 
     all = TRUE, 
     fancy = FALSE, 
     cex = 1.2)

# Add a main title
title(main = "Decision Trees at Different Maximum Depths", 
      outer = TRUE, 
      cex.main = 2, 
      line = -1)

# Reset to standard plotting parameters
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

# Close the PDF device
dev.off()


for (depth in depths) {
  # Create a very high resolution PDF for each tree
  pdf(sprintf("%s/decision_tree_depth_%d_%d_large.pdf", results.pwd, depth, year), 
      width = 14, height = 10, pointsize = 16)
  
  # Set margins for maximum clarity
  par(mar = c(1, 1, 3, 1))
  
  # Plot the tree with exaggerated text size and spacing
  plot(trees[[as.character(depth)]], 
       main = sprintf("Decision Tree (Depth = %d)", depth),
       uniform = TRUE,
       margin = 0.05,     # Very generous margin
       cex = 2.5)        # Extremely large node text
  
  # Add text with maximum readability
  text(trees[[as.character(depth)]], 
       cex = 1.7)        # Very large text
  
  # Reset margins
  par(mar = c(5, 4, 4, 2) + 0.1)
  
  # Close the PDF device
  dev.off()
}







