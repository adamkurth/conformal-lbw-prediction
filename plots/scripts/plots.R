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
} else {
  data.pwd <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/birthweight_data/rebin_without_2.5kg/"
  plots.pwd <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/plots/type2"
  results.pwd <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/results/results2"
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

# 
# library(ggplot2)
# library(plotly)
# library(data.tree)
# library(visNetwork)
# library(rpart.plot)
# library(DiagrammeR)
# library(igraph)
# rpart_to_dataframe <- function(model) {
#   if (!inherits(model, "rpart")) {
#     stop("Model must be an rpart object")
#   }
#   
#   # Extract the frame information
#   frame <- model$frame
#   
#   # Get splits information
#   splits <- model$splits
#   
#   # Create a dataframe for nodes
#   nodes_df <- data.frame(
#     id = rownames(frame),
#     label = row.names(frame),
#     var = frame$var,
#     n = frame$n,
#     dev = frame$dev,
#     yval = sapply(1:nrow(frame), function(i) {
#       if (is.leaf(frame[i,])) {
#         return(paste("Node", rownames(frame)[i], "n =", frame$n[i]))
#       } else {
#         return(paste(frame$var[i]))
#       }
#     }),
#     terminal = is.leaf(frame),
#     stringsAsFactors = FALSE
#   )
#   
#   # Create a dataframe for edges
#   # In rpart, if node i is not a leaf, then its children are 2i and 2i+1
#   edges_df <- data.frame()
#   for (i in nodes_df$id) {
#     i_num <- as.numeric(i)
#     left_child <- 2 * i_num
#     right_child <- 2 * i_num + 1
#     
#     # Check if children exist in the nodes
#     if (as.character(left_child) %in% nodes_df$id) {
#       # Get the split value
#       var_name <- nodes_df$var[nodes_df$id == i]
#       
#       # If not a leaf
#       if (var_name != "<leaf>") {
#         # Left child edge
#         edges_df <- rbind(edges_df, data.frame(
#           from = i,
#           to = as.character(left_child),
#           label = "yes",
#           stringsAsFactors = FALSE
#         ))
#         
#         # Right child edge
#         edges_df <- rbind(edges_df, data.frame(
#           from = i,
#           to = as.character(right_child),
#           label = "no",
#           stringsAsFactors = FALSE
#         ))
#       }
#     }
#   }
#   
#   return(list(nodes = nodes_df, edges = edges_df))
# }
# plot_interactive_tree <- function(model) {
#   tree_data <- rpart_to_dataframe(model)
#   
#   # Create a list of nodes with customized appearance
#   nodes <- data.frame(
#     id = tree_data$nodes$id,
#     label = sapply(1:nrow(tree_data$nodes), function(i) {
#       if (tree_data$nodes$terminal[i]) {
#         return(paste("Node", tree_data$nodes$id[i], "\nn =", tree_data$nodes$n[i]))
#       } else {
#         return(paste(tree_data$nodes$var[i]))
#       }
#     }),
#     shape = ifelse(tree_data$nodes$terminal, "box", "ellipse"),
#     color = ifelse(tree_data$nodes$terminal, 
#                    list(background = "#E8F8F5", border = "#1ABC9C"), 
#                    list(background = "#D4E6F1", border = "#3498DB")),
#     font = list(size = 14),
#     size = 25
#   )
#   
#   # Create a list of edges with customized appearance
#   edges <- data.frame(
#     from = tree_data$edges$from,
#     to = tree_data$edges$to,
#     label = tree_data$edges$label,
#     arrows = "to",
#     color = "#85929E",
#     width = 2,
#     font = list(size = 12)
#   )
#   
#   # Create the network visualization
#   visNetwork(nodes, edges, width = "100%", height = "800px") %>%
#     visLayout(hierarchical = list(
#       enabled = TRUE,
#       direction = "UD",
#       sortMethod = "directed",
#       levelSeparation = 100
#     )) %>%
#     visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
#                nodesIdSelection = TRUE) %>%
#     visInteraction(navigationButtons = TRUE)
# }
# plot_enhanced_rpart <- function(model) {
#   # Set up a high-resolution PNG device for better quality
#   # png("enhanced_tree_plot.png", width = 2000, height = 1600, res = 300)
#   
#   # Custom palette for better visualization
#   custom_palette <- colorRampPalette(c("#D4E6F1", "#2E86C1", "#1A5276"))(20)
#   
#   rpart.plot(model, 
#              type = 4,            # Type of plot (4 shows the splits)
#              extra = 101,         # Display the number of observations at each node
#              box.palette = custom_palette, # Custom colors
#              shadow.col = "gray", # Add shadows for depth
#              nn = TRUE,           # Display node numbers
#              fallen.leaves = TRUE, # Align terminal nodes
#              branch = 0.5,        # Set branch length
#              branch.lty = 1,      # Branch line type
#              branch.col = "gray30", # Branch color
#              split.cex = 1.2,     # Size of split text
#              split.box.col = "lightgray", # Split box color
#              split.border.col = "darkgray", # Split box border
#              split.round = 0.5,   # Round corners of split boxes
#              cex = 0.8)           # Overall text size
#   
#   # dev.off()
# }
# plot_ggplot_tree <- function(model) {
#   # Extract data from the model
#   rpart_frame <- model$frame
#   rpart_split <- model$splits
#   
#   # Convert rownames to numeric for easier manipulation
#   node_ids <- as.numeric(rownames(rpart_frame))
#   
#   # Create a clean nodes dataframe (excluding leaf nodes)
#   decision_nodes <- rpart_frame[rpart_frame$var != "<leaf>", ]
#   decision_node_ids <- as.numeric(rownames(decision_nodes))
#   
#   # Extract variable names and create labels
#   var_names <- as.character(decision_nodes$var)
#   
#   # Create a nodes dataframe with improved labels
#   nodes <- data.frame(
#     id = decision_node_ids,
#     var = var_names,
#     n = decision_nodes$n,
#     stringsAsFactors = FALSE
#   )
#   
#   # Create simplified tree structure
#   tree_depth <- ceiling(log2(max(node_ids) + 1))
#   max_nodes_at_depth <- 2^(tree_depth - 1)
#   
#   # Calculate x-coordinates based on node position in full binary tree
#   nodes$depth <- floor(log2(nodes$id))
#   
#   # Get the count of nodes at each depth
#   depth_counts <- table(nodes$depth)
#   max_depth <- max(nodes$depth)
#   
#   # Create layout matrix for nodes
#   x_coords <- list()
#   for (d in 0:max_depth) {
#     nodes_at_depth <- nodes[nodes$depth == d, ]
#     if (nrow(nodes_at_depth) > 0) {
#       n_nodes <- nrow(nodes_at_depth)
#       x_coords[[d+1]] <- seq(0, 1, length.out = n_nodes + 2)[2:(n_nodes+1)]
#     }
#   }
#   
#   # Assign x and y coordinates
#   nodes$y <- -nodes$depth  # Negative to put root at top
#   
#   # Assign x coordinates for each node based on its position at its depth
#   for (d in 0:max_depth) {
#     depth_nodes <- nodes[nodes$depth == d, ]
#     if (nrow(depth_nodes) > 0) {
#       # Sort nodes by ID to maintain binary tree structure
#       depth_nodes <- depth_nodes[order(depth_nodes$id), ]
#       # Assign coordinates
#       coords <- x_coords[[d+1]]
#       nodes$x[nodes$depth == d] <- coords
#     }
#   }
#   
#   # Create connections between nodes
#   edges <- data.frame()
#   for (i in nodes$id) {
#     left_child <- 2 * i
#     right_child <- 2 * i + 1
#     
#     # Check if children exist in the nodes
#     if (left_child %in% nodes$id) {
#       # Left child edge
#       edges <- rbind(edges, data.frame(
#         from = i,
#         to = left_child,
#         direction = "left",
#         stringsAsFactors = FALSE
#       ))
#     }
#     
#     if (right_child %in% nodes$id) {
#       # Right child edge
#       edges <- rbind(edges, data.frame(
#         from = i,
#         to = right_child,
#         direction = "right",
#         stringsAsFactors = FALSE
#       ))
#     }
#   }
#   
#   # Add coordinates to edges
#   if (nrow(edges) > 0) {
#     edges$x <- nodes$x[match(edges$from, nodes$id)]
#     edges$y <- nodes$y[match(edges$from, nodes$id)]
#     edges$xend <- nodes$x[match(edges$to, nodes$id)]
#     edges$yend <- nodes$y[match(edges$to, nodes$id)]
#   }
#   
#   # Create the ggplot
#   p <- ggplot() +
#     # Add nice grid for reference (optional)
#     theme_minimal() +
#     
#     # Draw edges with different colors for left/right branches
#     geom_segment(data = edges, 
#                  aes(x = x, y = y, xend = xend, yend = yend,
#                      color = direction, linetype = direction),
#                  size = 1) +
#     
#     # Draw nodes as rectangles
#     geom_rect(data = nodes,
#               aes(xmin = x - 0.03, xmax = x + 0.03,
#                   ymin = y - 0.15, ymax = y + 0.15),
#               fill = "#D6EAF8", color = "#2E86C1", alpha = 0.8) +
#     
#     # Add variable names inside nodes
#     geom_text(data = nodes,
#               aes(x = x, y = y, label = var),
#               size = 3.5, fontface = "bold") +
#     
#     # Add sample size under each node
#     geom_text(data = nodes,
#               aes(x = x, y = y + 0.25, 
#                   label = paste("n =", n)),
#               size = 2.8, color = "#566573") +
#     
#     # Add branch labels (yes/no or condition)
#     geom_text(data = edges,
#               aes(x = (x + xend) / 2, 
#                   y = (y + yend) / 2,
#                   label = ifelse(direction == "left", "yes", "no")),
#               size = 2.5, color = "#566573", vjust = 1.5) +
#     
#     # Set colors for the branches
#     scale_color_manual(values = c("left" = "#3498DB", "right" = "#E74C3C")) +
#     scale_linetype_manual(values = c("left" = "solid", "right" = "solid")) +
#     
#     # Clean up the plot
#     theme(
#       panel.grid = element_blank(),
#       axis.text = element_blank(),
#       axis.title = element_blank(),
#       axis.ticks = element_blank(),
#       legend.position = "none"
#     ) +
#     
#     # Set proper limits to avoid clipping
#     xlim(-0.05, 1.05) +
#     
#     # Add title
#     labs(title = "Dirichlet-Multinomial Decision Tree",
#          subtitle = "Decision Nodes Only")
#   
#   return(p)
# }
# 
# # df <- rpart_to_dataframe(dm.tree)
# # plot_enhanced_rpart(dm.tree)
# # plot_interactive_tree(dm.tree)
# plot_ggplot_tree(dm.tree)



# 
# 
# plot_bootstrap_results <- function(Yhat, lower, upper) {
#   # Create a dataframe for the heatmap
#   bootstrap_df <- data.frame(
#     Observation = 1:nrow(Yhat),
#     Mean_LBW_Prob = rowMeans(Yhat),
#     Lower_CI = lower,
#     Upper_CI = upper
#   )
#   
#   # Calculate width of confidence intervals
#   bootstrap_df$CI_Width <- bootstrap_df$Upper_CI - bootstrap_df$Lower_CI
#   
#   # Sort by mean LBW probability
#   bootstrap_df <- bootstrap_df[order(bootstrap_df$Mean_LBW_Prob, decreasing = TRUE), ]
#   
#   # Create the heatmap plot
#   p <- ggplot(bootstrap_df, aes(x = reorder(Observation, -Mean_LBW_Prob), y = 1)) +
#     geom_tile(aes(fill = Mean_LBW_Prob)) +
#     geom_errorbar(aes(ymin = 0.7, ymax = 1.3, width = 0.7, 
#                       color = CI_Width), alpha = 0.7) +
#     scale_fill_gradient2(low = "white", mid = "steelblue", high = "darkblue",
#                          midpoint = 0.5, name = "LBW Probability") +
#     scale_color_gradient(low = "green", high = "red", name = "CI Width") +
#     theme_minimal() +
#     theme(axis.text.x = element_blank(),
#           axis.ticks.x = element_blank(),
#           panel.grid = element_blank()) +
#     labs(title = "Bootstrap Estimates of Low Birth Weight Probabilities",
#          subtitle = "Ordered by Mean Probability with 95% Confidence Intervals",
#          x = "Observations (Ordered by LBW Probability)",
#          y = "")
#   
#   return(p)
# }
# 
# plot_bootstrap_results(Yhat, lower, upper)



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
#         names.arg = paste("Cat", 1:length(alphavec)),
#         col = "skyblue",
#         border = "blue",
#         main = "Informed Dirichlet Prior (alphavec)",
#         xlab = "Birth Weight Categories",
#         ylab = "Prior Probability",
#         ylim = c(0, max(alphavec) * 1.2))
# dev.off()

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









#-------------------------------------------------------------
# Visualize Bootstrap Results for Publication
#-------------------------------------------------------------
library(ggplot2)
library(tidyr)
library(dplyr)
library(viridis)      # For better color palettes
library(gridExtra)    # For arranging multiple plots
library(ggthemes)     # For publication-ready themes
library(scales)       # For better axis formatting

# Set theme for publication-quality plots
publication_theme <- theme_minimal(base_size = 12, base_family = "Helvetica") +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 15)),
    plot.caption = element_text(size = 9, hjust = 0, margin = margin(t = 10)),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "bottom",
    legend.box = "horizontal",
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95"),
    panel.border = element_rect(fill = NA, color = "gray80", size = 0.5),
    strip.text = element_text(size = 11, face = "bold")
  )

# Set default saving parameters
save_plot <- function(filename, plot, width = 7, height = 5, dpi = 300) {
  # Save as PDF
  ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi
  )
  
  # Also save as PNG for easy viewing
  ggsave(
    filename = gsub("\\.pdf$", ".png", filename),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi
  )
}

# 0. Internal consistency check (if necessary)
if (!exists("coverage.check")) {
  coverage.check <- matrix(FALSE, nrow = nrow(pred.results), ncol = B)
  for(b in 1:B) {
    # for each bootstrap sample, check if predictions are within CI bounds
    coverage.check[, b] <- (Yhat.all[, b] >= pred.results$lower.ci) & 
      (Yhat.all[, b] <= pred.results$upper.ci)
  }
  empirical.coverage <- rowMeans(coverage.check)
  pred.results$empirical.coverage <- empirical.coverage
}


# 2. Create combined plot with both intervals
interval_plot <- ggplot(pred.results, aes(x = 1:nrow(pred.results))) +
  # Add confidence intervals (narrower band)
  geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), 
              fill = "#3366CC", alpha = 0.4) +
  # Add prediction intervals (wider band)
  geom_ribbon(aes(ymin = lower.pi, ymax = upper.pi), 
              fill = "#CC3333", alpha = 0.2) +
  # Add the point predictions
  geom_point(aes(y = bagged.pred), size = 1.5, color = "black") +
  # Add horizontal reference line at mean prediction
  geom_hline(yintercept = mean(pred.results$bagged.pred), 
             linetype = "dashed", color = "gray50") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "95% Confidence and Prediction Intervals",
    subtitle = paste("Based on", B, "bootstrap samples"),
    x = "Observation Index", 
    y = "Probability",
    caption = "Blue band: 95% confidence intervals | Red band: 95% prediction intervals"
  )

# 3. Empirical Coverage Plot
coverage_plot <- ggplot(pred.results, aes(x = empirical.coverage)) +
  geom_histogram(bins = 50, fill = "#3366CC", color = "white", alpha = 0.8) +
  geom_vline(xintercept = 0.95, color = "#CC3333", linetype = "dashed", size = 1) +
  geom_density(aes(y = ..count.. * 0.8), color = "#000066", size = 0.8, alpha = 0.2) +
  scale_x_continuous(limits = c(0.90, 1.0), labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Empirical Coverage Rates of 95% Confidence Intervals",
    subtitle = paste("Based on", B, "bootstrap samples"),
    x = "Coverage Rate", 
    y = "Count",
    caption = "Red dashed line indicates nominal 95% coverage level"
  )

# 4. Combine plots using patchwork or gridExtra
if (requireNamespace("patchwork", quietly = TRUE)) {
  combined_plot <- interval_plot / coverage_plot
  save_plot(sprintf("%s/intervals_coverage_combined_pub.pdf", save.path), combined_plot, 
            width = 8, height = 10)
} else {
  # Save plots individually
  save_plot(sprintf("%s/prediction_confidence_intervals_pub.pdf", save.path), interval_plot)
  save_plot(sprintf("%s/empirical_coverage_pub.pdf", plots.pwd), coverage_plot)
}






ci_interval_plot <- ggplot(pred.results, aes(x = 1:nrow(pred.results))) +
  # Add confidence intervals band
  geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), 
              fill = "#3366CC", alpha = 0.4) +
  # Add the point predictions
  geom_point(aes(y = bagged.pred), size = 1.5, color = "black") +
  # Add horizontal reference line at mean prediction
  geom_hline(yintercept = mean(pred.results$bagged.pred), 
             linetype = "dashed", color = "gray50") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "95% Confidence Intervals for Low Birth Weight Predictions",
    subtitle = paste("Based on", B, "bootstrap samples"),
    x = "Observation Index", 
    y = "Probability",
    caption = "Blue band: 95% confidence intervals | Black points: Bagged predictions"
  ) +
  theme_minimal(base_size = 12, base_family = "Helvetica") +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 15)),
    plot.caption = element_text(size = 9, hjust = 0, margin = margin(t = 10)),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95"),
    panel.border = element_rect(fill = NA, color = "gray80", size = 0.5)
  )

ci_interval_sorted_plot <- ggplot(pred.results[order(pred.results$bagged.pred),], 
                                  aes(x = 1:nrow(pred.results))) +
  geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), 
              fill = "#3366CC", alpha = 0.4) +
  geom_point(aes(y = bagged.pred), size = 1.5, color = "black") +
  geom_hline(yintercept = mean(pred.results$bagged.pred), 
             linetype = "dashed", color = "gray50") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "95% Confidence Intervals (Sorted by Prediction Value)",
    subtitle = paste("Based on", B, "bootstrap samples"),
    x = "Sorted Observation Index", 
    y = "Probability",
    caption = "Blue band: 95% confidence intervals | Observations sorted by increasing prediction value"
  ) +
  theme_minimal(base_size = 12, base_family = "Helvetica") +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 15)),
    plot.caption = element_text(size = 9, hjust = 0, margin = margin(t = 10)),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95"),
    panel.border = element_rect(fill = NA, color = "gray80", size = 0.5)
  )


save_plot(sprintf("%s/ci_intervals_pub.pdf", results.pwd), ci_interval_plot, 
          width = 8, height = 5)
save_plot(sprintf("%s/ci_intervals_sorted_pub.pdf", plots.pwd), ci_interval_sorted_plot, 
          width = 8, height = 5)

combined_plot <- ci_interval_plot / ci_interval_sorted_plot
save_plot(sprintf("%s/ci_intervals_combined_pub.pdf", plots.pwd), combined_plot, 
          width = 8, height = 10)





# 2. Calibration Plot (Bagged vs. OOB)
calibration_plot <- ggplot(calibration.data, aes(x = bagged, y = oob)) +
  # Background panel without geom_rect
  theme_classic() +
  theme(panel.background = element_rect(fill = "gray95", color = NA)) +
  
  # Simple scatter plot with semi-transparency for density
  geom_point(alpha = 0.4, color = "#000066", size = 1.5) +
  
  # Perfect calibration line
  geom_abline(intercept = 0, slope = 1, color = "#CC3333", 
              linetype = "dashed", size = 1) +
  
  # Loess smooth line
  geom_smooth(method = "loess", color = "#000066", fill = "#3366CC", 
              alpha = 0.2, size = 1.2) +
  
  # Add marginal density at the axes
  geom_rug(sides = "b", alpha = 0.2, color = "#000066") +
  geom_rug(sides = "l", alpha = 0.2, color = "#000066") +
  
  # Format axes properly
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = c(0.01, 0.01)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = c(0.01, 0.01)) +
  
  # Equal coordinates for proper 45-degree line interpretation
  coord_equal() +
  
  # Labels
  labs(
    title = "Calibration Plot: Bagged vs. Out-of-Bag Predictions",
    subtitle = "Assessing prediction consistency across bootstrap samples",
    x = "Bagged Predictions (Probability)",
    y = "Out-of-Bag Predictions (Probability)",
    caption = "Red dashed line shows perfect calibration"
  ) +
  
  # Custom theme elements for a more traditional look
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 15)),
    plot.caption = element_text(size = 9, hjust = 0, margin = margin(t = 10)),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    panel.border = element_rect(fill = NA, color = "black", size = 0.5),
    axis.line = element_line(color = "black")
  )

save_plot(sprintf("%s/calibration_plot_pub.pdf", save.path), calibration_plot, width = 8, height = 7)


# 3. Prediction Intervals Plot
# Create a subset for clearer visualization if needed
if (nrow(pred.results) > 50) {
  # Sort by bagged prediction and select representative nodes
  # Fix for the slice() error - use explicit indices instead of seq()
  pred_sorted <- pred.results %>% arrange(bagged.pred)
  indices <- round(seq(1, nrow(pred_sorted), length.out = 50))
  pred_subset <- pred_sorted[indices, ]
} else {
  pred_subset <- pred.results
}

# Rearrange data for plotting
pred_long <- pred_subset %>%
  mutate(node_id = as.factor(node)) %>%
  arrange(bagged.pred)

# Create custom node ordering
pred_long$node_id <- factor(pred_long$node_id, levels = pred_long$node_id[order(pred_long$bagged.pred)])

pred_interval_plot <- ggplot(pred_long, aes(x = node_id, y = bagged.pred)) +
  geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), fill = "#3366CC", alpha = 0.2) +
  geom_errorbar(aes(ymin = lower.ci, ymax = upper.ci), width = 0.3, color = "#3366CC", alpha = 0.7) +
  geom_point(size = 3, color = "#CC3333") +
  geom_line(aes(group = 1), color = "#CC3333", alpha = 0.5) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Bootstrap Prediction Intervals for Low Birth Weight Probability",
    subtitle = paste("95% confidence intervals based on", B, "bootstrap samples"),
    x = "Node (ordered by prediction value)", 
    y = "Predicted LBW Probability",
    caption = "Red points show bagged predictions; blue bands show 95% prediction intervals"
  ) +
  publication_theme +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

save_plot(sprintf("%s/pred_interval_pub.pdf", save.path), pred_interval_plot, width = 10, height = 7)




# 4. Uncertainty vs. Prediction Values
uncertainty_plot <- ggplot(pred.results, aes(x = bagged.pred, y = width)) +
  geom_point(alpha = 0.6, color = "#3366CC", size = 2) +
  geom_smooth(method = "loess", color = "#CC3333", fill = "#CC3333", alpha = 0.2) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Relationship Between Prediction Value and Uncertainty",
    subtitle = "Width of confidence interval vs. predicted probability",
    x = "Predicted LBW Probability", 
    y = "Width of 95% Prediction Interval",
    caption = "Red curve shows smoothed relationship trend"
  ) +
  publication_theme

save_plot(sprintf("%s/uncertainty_vs_prediction_pub.pdf", save.path), uncertainty_plot)

hist_data <- hist(pred.results$bagged.pred, plot = FALSE, breaks = 40)
max_count <- max(hist_data$counts)
density_data <- density(pred.results$bagged.pred)
max_density <- max(density_data$y)
scale_factor <- max_count / max_density


pred_dist_plot <- ggplot(pred.results, aes(x = bagged.pred)) + 
  # Single, clear histogram
  geom_histogram(bins = 30, 
                 fill = "#3366CC", 
                 color = "white", 
                 alpha = 0.8) +
  
  # Clean density curve with proper scaling
  geom_density(aes(y = after_stat(count)), 
               color = "#CC3333", 
               size = 1.2,
               fill = NA) +
  
  # Format x-axis as percentage
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  
  # Labels
  labs(
    title = "Distribution of Predicted Low Birth Weight Probabilities",
    x = "Predicted LBW Probability",
    y = "Count",
    caption = "Red line shows probability density curve"
  ) +
  
  # Simplified theme for better readability
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    plot.caption = element_text(size = 9, hjust = 0, margin = margin(t = 10)),
    panel.border = element_rect(fill = NA, color = "black", size = 1),
    axis.line = element_line(color = "black")
  )

# save_plot(sprintf("%s/pred_dist_pub.pdf", save.path), pred_dist_plot, width = 8, height = 6)


# 6. Bootstrap Distribution for Selected Nodes
# Select representative nodes based on prediction percentiles
if (exists("Yhat.all")) {
  quantiles <- quantile(pred.results$bagged.pred, probs = c(0.1, 0.3, 0.7, 0.9))
  
  # Find nodes closest to these quantiles
  nodes.to.plot <- sapply(quantiles, function(q) {
    which.min(abs(pred.results$bagged.pred - q))
  })
  
  # Create data frame for plotting
  node.dists <- data.frame()
  for (i in seq_along(nodes.to.plot)) {
    node <- nodes.to.plot[i]
    node_data <- data.frame(
      node = paste0("Node ", node, " (p = ", round(pred.results$bagged.pred[node] * 100, 1), "%)"),
      prediction = Yhat.all[node, ]
    )
    node.dists <- rbind(node.dists, node_data)
  }
  
  # Order factor levels by prediction value
  node.dists$node <- factor(
    node.dists$node,
    levels = unique(node.dists$node)[order(sapply(nodes.to.plot, function(n) pred.results$bagged.pred[n]))]
  )
  
  bootstrap_dist_plot <- ggplot(node.dists, aes(x = prediction, fill = node)) +
    geom_density(alpha = 0.7, color = NA) +
    facet_wrap(~ node, ncol = 2) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_fill_viridis_d(option = "plasma") +
    labs(
      title = "Bootstrap Distribution of Predictions for Selected Nodes",
      subtitle = "Density plots for representative nodes across prediction range",
      x = "Predicted LBW Probability", 
      y = "Density"
    ) +
    publication_theme +
    theme(legend.position = "none")
  
  save_plot(sprintf("%s/bootstrap_dist_nodes_pub.pdf", save.path), bootstrap_dist_plot, width = 8, height = 8)
}



# 7. Combined Risk Factors Plot
# Select key risk factors from the data
risk_factors <- c("sex", "mrace15", "dmar", "mager")
risk_factors <- risk_factors[risk_factors %in% colnames(pred.results)]

if (length(risk_factors) > 0) {
  # Create long-format data for plotting
  risk_long <- pred.results %>%
    select(bagged.pred, width, all_of(risk_factors)) %>%
    pivot_longer(
      cols = all_of(risk_factors),
      names_to = "risk_factor",
      values_to = "value"
    )
  
  # Create factor labels - corrected for proper interpretation
  factor_labels <- c(
    "sex" = "Sex",
    "mrace15" = "Race",
    "dmar" = "Marital Status",
    "mager" = "Maternal Age"
  )
  
  # Define value labels for each factor
  value_labels <- list(
    "sex" = c("0" = "Female", "1" = "Male"),
    "mrace15" = c("0" = "Non-Black", "1" = "Black"),
    "dmar" = c("0" = "Not Married", "1" = "Married"),
    "mager" = c("0" = "Age ≤ 33", "1" = "Age > 33")
  )
  
  # Create custom x-axis labeling function
  custom_labeller <- function(variable, value) {
    if (variable %in% names(value_labels)) {
      return(value_labels[[variable]][as.character(value)])
    }
    return(value)
  }
  
  # Create plot with corrected labels
  risk_plot <- ggplot(risk_long, aes(x = as.factor(value), y = bagged.pred)) +
    geom_boxplot(fill = "#3366CC", alpha = 0.6, outlier.shape = NA) +
    geom_jitter(width = 0.2, height = 0, alpha = 0.3, size = 1) +
    # Use labeller for facet labels, and custom labels for x-axis values
    facet_wrap(~ risk_factor, scales = "free_x", labeller = labeller(risk_factor = factor_labels)) +
    scale_x_discrete(labels = function(x) {
      # Get the current facet's risk factor
      current_panel <- ggplot_build(risk_plot)$layout$panel_layout$risk_factor[panel_select()]
      return(sapply(x, function(val) custom_labeller(current_panel, val)))
    }) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      title = "Predicted Low Birth Weight Probability by Risk Factors",
      subtitle = "Distribution of predictions across key demographic variables",
      x = "Risk Factor Value", 
      y = "Predicted LBW Probability",
      caption = "0/1 values represent: Female/Male, Non-Black/Black, Not Married/Married, Age ≤33/Age >33"
    ) +
    theme_minimal(base_size = 12, base_family = "Helvetica") +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 15)),
      plot.caption = element_text(size = 9, hjust = 0, margin = margin(t = 10)),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      strip.text = element_text(size = 11, face = "bold"),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95"),
      panel.border = element_rect(fill = NA, color = "gray80", size = 0.5)
    )
  
  # Simpler approach that works more reliably
  # Create a modified dataset with human-readable factor levels
  risk_long_labeled <- risk_long %>%
    mutate(
      value_label = case_when(
        risk_factor == "sex" & value == 0 ~ "Female",
        risk_factor == "sex" & value == 1 ~ "Male",
        risk_factor == "mrace15" & value == 0 ~ "Non-Black",
        risk_factor == "mrace15" & value == 1 ~ "Black",
        risk_factor == "dmar" & value == 0 ~ "Not Married",
        risk_factor == "dmar" & value == 1 ~ "Married",
        risk_factor == "mager" & value == 0 ~ "Age ≤ 33",
        risk_factor == "mager" & value == 1 ~ "Age > 33",
        TRUE ~ as.character(value)
      )
    )
  
  # Create plot with corrected labels using the modified dataset
  risk_plot_labeled <- ggplot(risk_long_labeled, aes(x = value_label, y = bagged.pred)) +
    geom_boxplot(fill = "#3366CC", alpha = 0.6, outlier.shape = NA) +
    geom_jitter(width = 0.2, height = 0, alpha = 0.3, size = 1) +
    facet_wrap(~ risk_factor, scales = "free_x", labeller = labeller(risk_factor = factor_labels)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      title = "Predicted Low Birth Weight Probability by Risk Factors",
      subtitle = "Distribution of predictions across key demographic variables",
      x = "Risk Factor Value", 
      y = "Predicted LBW Probability"
    ) +
    theme_minimal(base_size = 12, base_family = "Helvetica") +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 15)),
      plot.caption = element_text(size = 9, hjust = 0, margin = margin(t = 10)),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      strip.text = element_text(size = 11, face = "bold"),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95"),
      panel.border = element_rect(fill = NA, color = "gray80", size = 0.5)
    )
  
  # Save the plot with properly labeled factors
  save_plot(sprintf("%s/risk_factors_labeled_pub.pdf", plots.pwd), risk_plot_labeled, width = 10, height = 8)
}

# # Create LaTeX code for the figure
# latex_code <- paste("% Risk factors and predicted low birth weight probabilities",
#                     "\\begin{figure}[htbp]",
#                     "    \\centering",
#                     "    \\includegraphics[width=0.85\\textwidth]{chapters/chapter3/figures/boot/risk_factors_labeled_pub.pdf}",
#                     "    \\caption{Predicted Low Birth Weight Probability by Risk Factors. This figure shows the distribution of predicted probabilities across key demographic variables. Black mothers (mrace15=1) show notably higher predicted risk of low birth weight compared to non-Black mothers. Similarly, unmarried mothers (dmar=0) exhibit higher predicted risk than married mothers. The effect of sex appears more modest, while maternal age shows that mothers over 33 (mager=1) have slightly different risk distributions than younger mothers.}",
#                     "    \\label{fig:boot-risk-factors}",
#                     "\\end{figure}", sep = "\n")
# 
# # Print LaTeX code to console
# cat(latex_code)





# 8. Create Enhanced Waterfall Plot of Node-specific Predictions
# This plot shows progression of predictions across nodes with uncertainty
waterfall_data <- pred_subset %>%
  arrange(bagged.pred) %>%
  mutate(
    node_id = factor(node, levels = node[order(bagged.pred)]),
    lower_diff = bagged.pred - lower.ci,
    upper_diff = upper.ci - bagged.pred
  )

waterfall_plot <- ggplot(waterfall_data, aes(x = node_id, y = bagged.pred)) +
  geom_linerange(aes(ymin = lower.ci, ymax = upper.ci), 
                 color = "#3366CC", size = 0.8, alpha = 0.5) +
  geom_point(size = 2.5, color = "#CC3333") +
  geom_hline(yintercept = median(waterfall_data$bagged.pred), 
             linetype = "dashed", color = "darkgray") +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Waterfall Plot of Node Predictions with Uncertainty",
    subtitle = "Nodes ordered by increasing predicted probability",
    x = "Node", 
    y = "Predicted LBW Probability",
    caption = "Red points show point estimates; blue lines show 95% prediction intervals"
  ) +
  publication_theme +
  theme(axis.text.y = element_blank())

save_plot(sprintf("%s/waterfall_plot_pub.pdf", save.path), waterfall_plot, width = 8, height = 10)

# Print completion message
cat("Successfully created publication-quality visualizations in", save.path, "\n")



# library(plotly, quietly=TRUE)
# plot_ly(
#   data=dm.tree$frame,
#   x=~var,
#   y=~n,
#   type="bar",
#   marker=list(color="blue")
# ) %>%
#   layout(
#     title="DM Tree Structure",
#     xaxis=list(title="Variable"),
#     yaxis=list(title="Sample Size")
#   )
# )
# 
# library(ggplot2)
# library(tree)
# library(ggdendro)
# 
# data(cpus, package = "MASS")
# model <- tree(log10(perf) ~ syct + mmin + mmax + cach + chmin + chmax, data = cpus)
# tree_data <- dendro_data(model)
# p <- ggplot(segment(tree_data)) +
#   geom_segment(aes(x = x, y = y, xend = xend, yend = yend, size = n), 
#                colour = "blue", alpha = 0.5) +
#   scale_size("n") +
#   geom_text(data = label(tree_data), 
#             aes(x = x, y = y, label = label), vjust = -0.5, size = 3) +
#   geom_text(data = leaf_label(tree_data), 
#             aes(x = x, y = y, label = label), vjust = 0.5, size = 2) +
#   theme_dendro()
# 
# ggplotly(p)



# # Alternative simplified visualization function without rpart.plot dependency
# visualize_dm_tree_simple <- function(tree) {
#   # If tree has splits
#   if (nrow(tree$frame) > 1) {
#     # Basic graphics setup
#     par(mfrow=c(1,1), mar=c(1,1,2,1))
#     
#     # Use base R plot for the tree structure
#     plot(tree, uniform=TRUE, margin=0.1, main="Dirichlet-Multinomial Tree")
#     text(tree, use.n=TRUE, cex=0.8, fancy=TRUE)
#     
#     # Return success message
#     return("Tree visualization completed with base R plot")
#   } else {
#     cat("Tree has no splits - visualization skipped\n")
#     return(NULL)
#   }
# }
# 
# # More comprehensive visualization for publications
# visualize_dm_publication <- function(tree) {
#   par(mfrow=c(2,2))
#   
#   # 1. Tree structure with base R plot
#   plot(tree, uniform=TRUE, branch=0.5, compress=TRUE, 
#        margin=0.1, main="DM Tree Structure")
#   text(tree, use.n=TRUE, cex=0.7)
#   
#   # 2. Terminal nodes analysis
#   terminal_idx <- which(tree$frame$var == "<leaf>")
#   if (length(terminal_idx) > 0) {
#     terminal_nodes <- data.frame(
#       node_id = as.numeric(rownames(tree$frame)[terminal_idx]),
#       n = tree$frame$n[terminal_idx]
#     )
#     
#     # Sort by sample size
#     terminal_nodes <- terminal_nodes[order(terminal_nodes$n, decreasing=TRUE), ]
#     
#     # Get top nodes
#     top_nodes <- head(terminal_nodes, min(4, nrow(terminal_nodes)))
#     
#     # 3. Distribution for top nodes
#     if (nrow(top_nodes) > 0 && !is.null(tree$y)) {
#       n_cats <- ncol(tree$y)
#       dist_matrix <- matrix(0, nrow=nrow(top_nodes), ncol=n_cats)
#       rownames(dist_matrix) <- paste0("Node ", top_nodes$node_id)
#       
#       # Try to get category names
#       if (!is.null(colnames(tree$y))) {
#         colnames(dist_matrix) <- colnames(tree$y)
#       } else {
#         colnames(dist_matrix) <- paste0("Cat", 1:n_cats)
#       }
#       
#       # Fill in distributions
#       for (i in 1:nrow(top_nodes)) {
#         node_id <- top_nodes$node_id[i]
#         obs_in_node <- which(tree$where == node_id)
#         
#         if (length(obs_in_node) > 0) {
#           node_counts <- colSums(tree$y[obs_in_node, , drop=FALSE])
#           dist_matrix[i,] <- node_counts / sum(node_counts)
#         }
#       }
#       
#       # Plot distributions as barplot
#       barplot(t(dist_matrix), beside=TRUE, 
#               main="Node Distributions",
#               col=rainbow(n_cats),
#               legend.text=colnames(dist_matrix),
#               args.legend=list(x="topright", cex=0.6))
#     }
#     
#     # 4. Variable importance (if available)
#     if (!is.null(tree$variable.importance)) {
#       var_imp <- tree$variable.importance
#       var_imp <- sort(var_imp, decreasing=TRUE)
#       top_vars <- head(var_imp, min(10, length(var_imp)))
#       
#       barplot(top_vars, 
#               main="Variable Importance",
#               las=2,
#               cex.names=0.7)
#     } else {
#       # Alternative: plot node sizes
#       barplot(top_nodes$n, 
#               names.arg=paste("Node", top_nodes$node_id),
#               main="Terminal Node Sizes",
#               las=2)
#     }
#   } else {
#     # Handle case with no terminal nodes
#     plot(1, type="n", axes=FALSE, xlab="", ylab="")
#     text(1, 1, "No terminal nodes found", cex=1.2)
#   }
#   
#   par(mfrow=c(1,1))
# }
# 
# # Function to try different plotting approaches
# plot_dm_tree <- function(tree) {
#   # Try to fix the tree for rpart.plot
#   fixed_tree <- fix_dm_tree_for_plotting(tree)
#   
#   cat("\n--- Simple base R tree plot ---\n")
#   visualize_dm_tree_simple(fixed_tree)
#   
#   cat("\n--- Publication-ready visualizations ---\n")
#   visualize_dm_publication(fixed_tree)
#   
#   # Only try rpart.plot if the tree seems compatible
#   tryCatch({
#     cat("\n--- Attempting rpart.plot visualization ---\n")
#     rpart.plot(fixed_tree, 
#                extra=106,  # Show node numbers and percentages
#                box.palette="auto", 
#                main="Dirichlet-Multinomial Tree")
#     cat("rpart.plot visualization succeeded\n")
#   }, error=function(e) {
#     cat("rpart.plot visualization failed:", e$message, "\n")
#     cat("Falling back to base R plot\n")
#   })
#   
#   return(fixed_tree)
# }
# 
# dm.tree <- fix_dm_tree_for_plotting(dm.tree)

# visualize_dm_tree_simple(dm.tree)
# visualize_dm_publication(dm.tree)
# plot_dm_tree(dm.tree)
