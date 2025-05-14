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


library(patchwork)  # For combining plots
library(dplyr)      # For data manipulation
high_risk_df <- data.frame(
  Category = factor(paste0("C", 1:K), levels = paste0("C", 1:K)),
  Probability = high.risk.probs,
  Lower = high.risk.lwr,
  Upper = high.risk.upr,
  Group = "High Risk Subgroup"
)

low_risk_df <- data.frame(
  Category = factor(paste0("C", 1:K), levels = paste0("C", 1:K)),
  Probability = low.risk.probs,
  Lower = low.risk.lwr,
  Upper = low.risk.upr,
  Group = "Low Risk Subgroup"
)

# Combine data for easy theming
combined_df <- rbind(high_risk_df, low_risk_df)

# Create High Risk Plot
p1 <- ggplot(high_risk_df, aes(x = Category, y = Probability)) +
  geom_col(fill = "#FF6666", color = "#990000", width = 0.7) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), 
                width = 0.2, color = "black", linewidth = 0.8) +
  labs(title = "High Risk Subgroup",
       x = "",
       y = "Probability") +
  scale_y_continuous(limits = c(0, max(combined_df$Upper) * 1.1),
                     expand = expansion(mult = c(0, 0.05))) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    panel.grid = element_blank(),
    plot.margin = margin(5, 10, 5, 10)
  )

# Create Low Risk Plot  
p2 <- ggplot(low_risk_df, aes(x = Category, y = Probability)) +
  geom_col(fill = "#377EB8", color = "#00008B", width = 0.7) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), 
                width = 0.2, color = "black", linewidth = 0.8) +
  labs(title = "Low Risk Subgroup",
       x = "",
       y = "Probability") +
  scale_y_continuous(limits = c(0, max(combined_df$Upper) * 1.1),
                     expand = expansion(mult = c(0, 0.05))) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    panel.grid = element_blank(),
    plot.margin = margin(5, 10, 5, 10)
  )

# Combine plots with patchwork
combined_plot <- p1 / p2 +
  plot_annotation(
    title = paste0("Birth Weight Category Bootstrap Probabilities"),
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      plot.caption = element_text(size = 10, color = "darkgray", hjust = 0.5),
      plot.margin = margin(20, 20, 20, 20)
    )
  )

# Save the plot
ggsave(
  filename = sprintf("%s/high_low_risk_pred_%d.pdf", boot.pwd, year),
  plot = combined_plot,
  width = 8.5, 
  height = 8,
  units = "in",
  dpi = 300
)

# Print the result table
cat("\nMean probabilities for all categories:\n")
result.table <- data.frame(
  cat = 1:K,
  high.risk.prob = high.risk.probs, 
  high.risk.lwr = high.risk.lwr,
  high.risk.upr = high.risk.upr,
  low.risk.prob = low.risk.probs,
  low.risk.lwr = low.risk.lwr,
  low.risk.upr = low.risk.upr
)
print(result.table)






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







