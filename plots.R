rm(list=ls())
library(rpart)
library(rpart.plot)
year <- 2021

# Determine type based on data path
use_without_2_5kg <- TRUE # Set to TRUE for type 2 (without 2.5kg), FALSE for type 1
type <- ifelse(use_without_2_5kg, 2, 1)

# Set appropriate data path based on type
if(type == 1) {
  data.pwd <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/birthweight_data/rebin/"
  plots.pwd <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/plots/type1"
} else {
  data.pwd <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/birthweight_data/rebin_without_2.5kg/"
  plots.pwd <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/plots/type2"
}

# data.pwd <- sprintf("/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/birthweight_data/rebin/")
load(sprintf("%s/birthweight_data_rebin_%d.RData", data.pwd, year))
load(sprintf("%s/informed_prior_%d.RData", data.pwd, 2020))
load(sprintf("%s/bootstrap_results_%d.RData", data.pwd, year))
load(sprintf("%s/dm_tree_%d.RData", data.pwd, year))

pdf(sprintf("%s/dm_tree_%d.pdf", plots.pwd, type), width = 10, height = 8)
plot(dm.tree, main="Dirichlet−Multinomial Decision Tree", uniform=TRUE, margin=0.1)
text(dm.tree, use.n=TRUE, cex=0.8, all=TRUE, pretty=0, xpd=NA, digits=3, font=3)
dev.off()


library(ggplot2)
library(plotly)
library(data.tree)
library(visNetwork)
library(rpart.plot)
library(DiagrammeR)
library(igraph)
rpart_to_dataframe <- function(model) {
  if (!inherits(model, "rpart")) {
    stop("Model must be an rpart object")
  }
  
  # Extract the frame information
  frame <- model$frame
  
  # Get splits information
  splits <- model$splits
  
  # Create a dataframe for nodes
  nodes_df <- data.frame(
    id = rownames(frame),
    label = row.names(frame),
    var = frame$var,
    n = frame$n,
    dev = frame$dev,
    yval = sapply(1:nrow(frame), function(i) {
      if (is.leaf(frame[i,])) {
        return(paste("Node", rownames(frame)[i], "n =", frame$n[i]))
      } else {
        return(paste(frame$var[i]))
      }
    }),
    terminal = is.leaf(frame),
    stringsAsFactors = FALSE
  )
  
  # Create a dataframe for edges
  # In rpart, if node i is not a leaf, then its children are 2i and 2i+1
  edges_df <- data.frame()
  for (i in nodes_df$id) {
    i_num <- as.numeric(i)
    left_child <- 2 * i_num
    right_child <- 2 * i_num + 1
    
    # Check if children exist in the nodes
    if (as.character(left_child) %in% nodes_df$id) {
      # Get the split value
      var_name <- nodes_df$var[nodes_df$id == i]
      
      # If not a leaf
      if (var_name != "<leaf>") {
        # Left child edge
        edges_df <- rbind(edges_df, data.frame(
          from = i,
          to = as.character(left_child),
          label = "yes",
          stringsAsFactors = FALSE
        ))
        
        # Right child edge
        edges_df <- rbind(edges_df, data.frame(
          from = i,
          to = as.character(right_child),
          label = "no",
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  return(list(nodes = nodes_df, edges = edges_df))
}
plot_interactive_tree <- function(model) {
  tree_data <- rpart_to_dataframe(model)
  
  # Create a list of nodes with customized appearance
  nodes <- data.frame(
    id = tree_data$nodes$id,
    label = sapply(1:nrow(tree_data$nodes), function(i) {
      if (tree_data$nodes$terminal[i]) {
        return(paste("Node", tree_data$nodes$id[i], "\nn =", tree_data$nodes$n[i]))
      } else {
        return(paste(tree_data$nodes$var[i]))
      }
    }),
    shape = ifelse(tree_data$nodes$terminal, "box", "ellipse"),
    color = ifelse(tree_data$nodes$terminal, 
                   list(background = "#E8F8F5", border = "#1ABC9C"), 
                   list(background = "#D4E6F1", border = "#3498DB")),
    font = list(size = 14),
    size = 25
  )
  
  # Create a list of edges with customized appearance
  edges <- data.frame(
    from = tree_data$edges$from,
    to = tree_data$edges$to,
    label = tree_data$edges$label,
    arrows = "to",
    color = "#85929E",
    width = 2,
    font = list(size = 12)
  )
  
  # Create the network visualization
  visNetwork(nodes, edges, width = "100%", height = "800px") %>%
    visLayout(hierarchical = list(
      enabled = TRUE,
      direction = "UD",
      sortMethod = "directed",
      levelSeparation = 100
    )) %>%
    visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
               nodesIdSelection = TRUE) %>%
    visInteraction(navigationButtons = TRUE)
}
plot_enhanced_rpart <- function(model) {
  # Set up a high-resolution PNG device for better quality
  # png("enhanced_tree_plot.png", width = 2000, height = 1600, res = 300)
  
  # Custom palette for better visualization
  custom_palette <- colorRampPalette(c("#D4E6F1", "#2E86C1", "#1A5276"))(20)
  
  rpart.plot(model, 
             type = 4,            # Type of plot (4 shows the splits)
             extra = 101,         # Display the number of observations at each node
             box.palette = custom_palette, # Custom colors
             shadow.col = "gray", # Add shadows for depth
             nn = TRUE,           # Display node numbers
             fallen.leaves = TRUE, # Align terminal nodes
             branch = 0.5,        # Set branch length
             branch.lty = 1,      # Branch line type
             branch.col = "gray30", # Branch color
             split.cex = 1.2,     # Size of split text
             split.box.col = "lightgray", # Split box color
             split.border.col = "darkgray", # Split box border
             split.round = 0.5,   # Round corners of split boxes
             cex = 0.8)           # Overall text size
  
  # dev.off()
}
plot_ggplot_tree <- function(model) {
  # Extract data from the model
  rpart_frame <- model$frame
  rpart_split <- model$splits
  
  # Convert rownames to numeric for easier manipulation
  node_ids <- as.numeric(rownames(rpart_frame))
  
  # Create a clean nodes dataframe (excluding leaf nodes)
  decision_nodes <- rpart_frame[rpart_frame$var != "<leaf>", ]
  decision_node_ids <- as.numeric(rownames(decision_nodes))
  
  # Extract variable names and create labels
  var_names <- as.character(decision_nodes$var)
  
  # Create a nodes dataframe with improved labels
  nodes <- data.frame(
    id = decision_node_ids,
    var = var_names,
    n = decision_nodes$n,
    stringsAsFactors = FALSE
  )
  
  # Create simplified tree structure
  tree_depth <- ceiling(log2(max(node_ids) + 1))
  max_nodes_at_depth <- 2^(tree_depth - 1)
  
  # Calculate x-coordinates based on node position in full binary tree
  nodes$depth <- floor(log2(nodes$id))
  
  # Get the count of nodes at each depth
  depth_counts <- table(nodes$depth)
  max_depth <- max(nodes$depth)
  
  # Create layout matrix for nodes
  x_coords <- list()
  for (d in 0:max_depth) {
    nodes_at_depth <- nodes[nodes$depth == d, ]
    if (nrow(nodes_at_depth) > 0) {
      n_nodes <- nrow(nodes_at_depth)
      x_coords[[d+1]] <- seq(0, 1, length.out = n_nodes + 2)[2:(n_nodes+1)]
    }
  }
  
  # Assign x and y coordinates
  nodes$y <- -nodes$depth  # Negative to put root at top
  
  # Assign x coordinates for each node based on its position at its depth
  for (d in 0:max_depth) {
    depth_nodes <- nodes[nodes$depth == d, ]
    if (nrow(depth_nodes) > 0) {
      # Sort nodes by ID to maintain binary tree structure
      depth_nodes <- depth_nodes[order(depth_nodes$id), ]
      # Assign coordinates
      coords <- x_coords[[d+1]]
      nodes$x[nodes$depth == d] <- coords
    }
  }
  
  # Create connections between nodes
  edges <- data.frame()
  for (i in nodes$id) {
    left_child <- 2 * i
    right_child <- 2 * i + 1
    
    # Check if children exist in the nodes
    if (left_child %in% nodes$id) {
      # Left child edge
      edges <- rbind(edges, data.frame(
        from = i,
        to = left_child,
        direction = "left",
        stringsAsFactors = FALSE
      ))
    }
    
    if (right_child %in% nodes$id) {
      # Right child edge
      edges <- rbind(edges, data.frame(
        from = i,
        to = right_child,
        direction = "right",
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Add coordinates to edges
  if (nrow(edges) > 0) {
    edges$x <- nodes$x[match(edges$from, nodes$id)]
    edges$y <- nodes$y[match(edges$from, nodes$id)]
    edges$xend <- nodes$x[match(edges$to, nodes$id)]
    edges$yend <- nodes$y[match(edges$to, nodes$id)]
  }
  
  # Create the ggplot
  p <- ggplot() +
    # Add nice grid for reference (optional)
    theme_minimal() +
    
    # Draw edges with different colors for left/right branches
    geom_segment(data = edges, 
                 aes(x = x, y = y, xend = xend, yend = yend,
                     color = direction, linetype = direction),
                 size = 1) +
    
    # Draw nodes as rectangles
    geom_rect(data = nodes,
              aes(xmin = x - 0.03, xmax = x + 0.03,
                  ymin = y - 0.15, ymax = y + 0.15),
              fill = "#D6EAF8", color = "#2E86C1", alpha = 0.8) +
    
    # Add variable names inside nodes
    geom_text(data = nodes,
              aes(x = x, y = y, label = var),
              size = 3.5, fontface = "bold") +
    
    # Add sample size under each node
    geom_text(data = nodes,
              aes(x = x, y = y + 0.25, 
                  label = paste("n =", n)),
              size = 2.8, color = "#566573") +
    
    # Add branch labels (yes/no or condition)
    geom_text(data = edges,
              aes(x = (x + xend) / 2, 
                  y = (y + yend) / 2,
                  label = ifelse(direction == "left", "yes", "no")),
              size = 2.5, color = "#566573", vjust = 1.5) +
    
    # Set colors for the branches
    scale_color_manual(values = c("left" = "#3498DB", "right" = "#E74C3C")) +
    scale_linetype_manual(values = c("left" = "solid", "right" = "solid")) +
    
    # Clean up the plot
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none"
    ) +
    
    # Set proper limits to avoid clipping
    xlim(-0.05, 1.05) +
    
    # Add title
    labs(title = "Dirichlet-Multinomial Decision Tree",
         subtitle = "Decision Nodes Only")
  
  return(p)
}

# df <- rpart_to_dataframe(dm.tree)
# plot_enhanced_rpart(dm.tree)
# plot_interactive_tree(dm.tree)
plot_ggplot_tree(dm.tree)



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
