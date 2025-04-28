rm(list=ls())
library(ggplot2)
library(dplyr)
library(forcats)
library(scales)
library(ggtext)  # For improved text rendering

year <- 2021

# Determine type based on data path
use_without_2_5kg <- TRUE  # Set to TRUE for type 2 (without 2.5kg), FALSE for type 1
type <- ifelse(use_without_2_5kg, 2, 1)

# Set appropriate data path based on type
if(type == 1) {
  data.pwd <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/birthweight_data/rebin/"
  model_label <- "Full Model"
} else {
  data.pwd <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/birthweight_data/rebin_without_2.5kg/"
  model_label <- "LBW-only Model"
}

# Set output directory for plots
plot.pwd <- sprintf("/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/plots/type%d", type)
dir.create(plot.pwd, showWarnings = FALSE, recursive = TRUE)

# Load data
load(sprintf("%s/birthweight_data_rebin_%d.RData", data.pwd, year))
load(sprintf("%s/informed_prior_%d.RData", data.pwd, 2020))
load(sprintf("%s/bootstrap_results_%d.RData", data.pwd, year))
load(sprintf("%s/dm_tree_%d.RData", data.pwd, year))

# Function to extract variable importance data from the tree
extract_variable_importance <- function(tree) {
  # Extract splits data as a data frame
  splits_data <- as.data.frame(tree$splits)
  
  # If the splits data is empty, return NULL
  if(is.null(splits_data) || nrow(splits_data) == 0) {
    cat("No splits data found in the tree object.\n")
    return(NULL)
  }
  
  # Add variable names as a column
  splits_data$variable <- rownames(splits_data)
  
  # Sort by improvement (descending)
  splits_data <- splits_data[order(-splits_data$improve), ]
  
  # Create variable importance summary
  var_importance <- aggregate(improve ~ variable, data = splits_data, sum)
  colnames(var_importance)[2] <- "total_improvement"
  var_importance <- var_importance[order(-var_importance$total_improvement), ]
  var_importance$rank <- 1:nrow(var_importance)
  
  return(var_importance)
}

# Create a publication-ready variable importance plot
create_publication_plot <- function(var_importance, model_type, output_dir) {
  # Create a mapping for clearer variable names - use abbreviated terms for space
  variable_mapping <- c(
    "mrace15" = "Mother's Race",
    "dmar" = "Marital Status",
    "cig_0" = "Cigarette Use",
    "precare5" = "Prenatal Care",
    "meduc" = "Mother's Education", 
    "sex" = "Child's Sex",
    "mager" = "Mother's Age"
  )
  
  # Apply mapping to variable names with improved formatting
  var_importance$display_name <- var_importance$variable
  for (i in 1:nrow(var_importance)) {
    var_base <- strsplit(as.character(var_importance$variable[i]), "\\.")[[1]][1]
    if (var_base %in% names(variable_mapping)) {
      if (grepl("\\.", var_importance$variable[i])) {
        level <- strsplit(as.character(var_importance$variable[i]), "\\.")[[1]][2]
        var_importance$display_name[i] <- sprintf("%s (%s)", variable_mapping[var_base], level)
      } else {
        var_importance$display_name[i] <- variable_mapping[var_base]
      }
    }
  }
  
  # Add percentages for context
  var_importance$percentage <- var_importance$total_improvement / sum(var_importance$total_improvement) * 100
  
  # Define groups based on natural breakpoints in the data
  # Use percentage-based cutoffs for better differentiation
  var_importance$group <- cut(var_importance$percentage,
                              breaks = c(-Inf, 1, 3, 10, Inf),
                              labels = c("< 1%", "1-3%", "3-10%", "> 10%"),
                              include.lowest = TRUE)
  
  # Create a visually appealing color palette (colorblind-friendly)
  # Use a professional blue/teal palette from colorbrewer
  color_palette <- c("#084594", "#2171b5", "#4292c6", "#6baed6")
  
  # Don't try to use specific fonts - use system default instead
  # This avoids the PostScript font database warnings
  
  # Calculate a reasonable figure height based on number of variables
  # Ensure adequate spacing between y-axis labels
  optimal_height <- max(8, min(nrow(var_importance) * 0.22, 14))
  
  # Create the enhanced bar plot with publication-quality styling
  p <- ggplot(var_importance, aes(x = reorder(display_name, total_improvement), 
                                  y = total_improvement,
                                  fill = group)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = color_palette, name = "Contribution") +
    coord_flip() +
    # Add rank labels for top variables
    geom_text(data = head(var_importance, 10),
              aes(label = sprintf("(#%d)", rank)),
              hjust = -0.2, size = 3.2, fontface = "bold") +
    # Format x-axis with commas for large numbers
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.1))) +
    labs(
      title = paste("Variable Importance in", model_type),
      subtitle = "Ranked by Improvement Score",
      caption = paste("Data from", year, "| Top variables labeled with rank | Colors indicate relative contribution"),
      x = NULL,
      y = "Improvement Score"
    ) +
    theme_minimal(base_size = 11) + # Removed base_family parameter
    theme(
      # Remove grid lines for cleaner appearance
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "gray90"),
      
      # Improve text elements
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "gray30"),
      plot.caption = element_text(size = 9, color = "gray30", hjust = 0),
      
      # Improve axis formatting
      axis.title.x = element_text(size = 11, margin = margin(t = 10)),
      axis.text.y = element_text(size = 9),
      
      # Position legend at top for better spacing
      legend.position = "top",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      
      # Add some breathing room
      plot.margin = margin(t = 20, r = 40, b = 20, l = 10)
    )
  
  # Save in high-quality formats for publication
  ggsave(file.path(output_dir, paste0("publication_variable_importance_", model_type, ".png")), 
         p, width = 8.5, height = optimal_height, dpi = 400)
  
  # Save as PDF for vector graphics
  ggsave(file.path(output_dir, paste0("publication_variable_importance_", model_type, ".pdf")), 
         p, width = 8.5, height = optimal_height)
  
  # Create a focused version showing only top 20 variables for better readability
  top20 <- head(var_importance, 20)
  
  p_top20 <- ggplot(top20, aes(x = reorder(display_name, total_improvement), 
                               y = total_improvement,
                               fill = group)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = color_palette, name = "Contribution") +
    coord_flip() +
    geom_text(aes(label = sprintf("#%d (%.1f%%)", rank, percentage)),
              hjust = -0.1, size = 3.5) +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.15))) +
    labs(
      title = paste("Top 20 Variables by Importance -", model_type),
      subtitle = "Ranked by Improvement Score with Percentage Contribution",
      caption = paste("Data from", year),
      x = NULL,
      y = "Improvement Score"
    ) +
    theme_minimal(base_size = 11) + # Removed base_family parameter
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "gray90"),
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "gray30"),
      plot.caption = element_text(size = 9, color = "gray30", hjust = 0),
      axis.title.x = element_text(size = 11, margin = margin(t = 10)),
      axis.text.y = element_text(size = 10, face = "bold"),
      legend.position = "top",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      plot.margin = margin(t = 20, r = 60, b = 20, l = 10)
    )
  
  # Save focused version
  ggsave(file.path(output_dir, paste0("publication_top20_importance_", model_type, ".png")), 
         p_top20, width = 8.5, height = 8, dpi = 400)
  
  ggsave(file.path(output_dir, paste0("publication_top20_importance_", model_type, ".pdf")), 
         p_top20, width = 8.5, height = 8)
  
  cat(sprintf("Created publication-ready variable importance plots for %s in %s\n", 
              model_type, output_dir))
  
  return(list(full = p, top20 = p_top20))
}

# Extract variable importance from the tree
var_importance <- extract_variable_importance(dm.tree)

# Create publication-ready plots
plots <- create_publication_plot(var_importance, model_label, plot.pwd)

cat("\nPublication-ready variable importance plots created successfully!\n")