# Modified Tree Split Visualization Script
# Creates visualization showing top variables used in tree splits

# Load required libraries
library(ggplot2)
library(dplyr)
library(forcats)
library(scales)
library(ggtext)

# Set parameters
year <- 2021
max_variables <- 20  # Limit visualization to top N variables

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

# Load tree data
load(sprintf("%s/dm_tree_%d.RData", data.pwd, year))

# Function to extract and process variable importance data from the tree
extract_tree_split_contribution <- function(tree, max_vars = Inf) {
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
  
  # Create variable importance summary by aggregating improvement by variable
  var_contribution <- aggregate(improve ~ variable, data = splits_data, sum)
  colnames(var_contribution)[2] <- "total_improvement"
  var_contribution <- var_contribution[order(-var_contribution$total_improvement), ]
  var_contribution$rank <- 1:nrow(var_contribution)
  
  # Calculate percentage contribution
  var_contribution$percentage <- var_contribution$total_improvement / sum(var_contribution$total_improvement) * 100
  
  # Limit to top N variables if specified
  if(nrow(var_contribution) > max_vars) {
    var_contribution <- head(var_contribution, max_vars)
  }
  
  return(var_contribution)
}

# Function to create publication-ready visualization of tree split contributions
create_tree_split_visualization <- function(var_contribution, model_type, output_dir) {
  # Create a mapping for clearer variable names
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
  var_contribution$display_name <- var_contribution$variable
  for (i in 1:nrow(var_contribution)) {
    var_base <- strsplit(as.character(var_contribution$variable[i]), "\\.")[[1]][1]
    if (var_base %in% names(variable_mapping)) {
      if (grepl("\\.", var_contribution$variable[i])) {
        level <- strsplit(as.character(var_contribution$variable[i]), "\\.")[[1]][2]
        var_contribution$display_name[i] <- sprintf("%s (%s)", variable_mapping[var_base], level)
      } else {
        var_contribution$display_name[i] <- variable_mapping[var_base]
      }
    }
  }
  
  # Define groups based on natural breakpoints in the data
  var_contribution$group <- cut(var_contribution$percentage,
                                breaks = c(-Inf, 1, 3, 10, Inf),
                                labels = c("< 1%", "1-3%", "3-10%", "> 10%"),
                                include.lowest = TRUE)
  
  # Create a visually appealing color palette (colorblind-friendly)
  color_palette <- c("#084594", "#2171b5", "#4292c6", "#6baed6")
  
  # Calculate optimal figure height based on number of variables
  optimal_height <- max(8, min(nrow(var_contribution) * 0.25, 14))
  
  # Create the enhanced bar plot with publication-quality styling
  p <- ggplot(var_contribution, aes(x = reorder(display_name, total_improvement), 
                                    y = total_improvement,
                                    fill = group)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = color_palette, name = "Contribution") +
    coord_flip() +
    # Add rank labels for all variables
    geom_text(aes(label = sprintf("(#%d)", rank)),
              hjust = -0.2, size = 3.2, fontface = "bold") +
    # Format x-axis with commas for large numbers
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.1))) +
    
      labs(
        title = paste("Key Split Variables in", model_type, "Decision Tree"),
        subtitle = "Ranked by Total Split Contribution",
        caption = paste("Data from", year, 
                        "| Values represent summed improvement in tree fit when variable is used for splitting"),  # Change "Numbers indicate variable" to "Numbers indicate variable rank"
        x = NULL,
        y = "Split Contribution Value"
    ) +
    theme_minimal(base_size = 11) +
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
  filename_base <- paste0("tree_split_contribution_top", nrow(var_contribution), "_", model_type)
  
  # Save as PNG for screen viewing
  ggsave(file.path(output_dir, paste0(filename_base, ".png")), 
         p, width = 8.5, height = optimal_height, dpi = 400)
  
  # Save as PDF for vector graphics (better for publication)
  ggsave(file.path(output_dir, paste0(filename_base, ".pdf")), 
         p, width = 8.5, height = optimal_height)
  
  cat(sprintf("Created tree split contribution visualization for %s with top %d variables in %s\n", 
              model_type, nrow(var_contribution), output_dir))
  
  return(p)
}

# Main execution
cat("Extracting split contribution data from decision tree...\n")
var_contribution <- extract_tree_split_contribution(dm.tree, max_variables)

if(!is.null(var_contribution)) {
  cat(sprintf("Found %d variables, limiting to top %d...\n", 
              nrow(var_contribution), min(nrow(var_contribution), max_variables)))
  
  # Create visualization
  plot <- create_tree_split_visualization(var_contribution, model_label, plot.pwd)
  
  cat("\nVisualization created successfully!\n")
  cat("Number of variables included:", nrow(var_contribution), "\n")
  cat("Top 5 variables by contribution:\n")
  
  # Display top 5 variables with proper column selection
  if(nrow(var_contribution) >= 5) {
    top5 <- head(var_contribution, 5)
    print(top5[, c("rank", "variable", "total_improvement", "percentage")])
  } else {
    print(var_contribution[, c("rank", "variable", "total_improvement", "percentage")])
  }
} else {
  cat("ERROR: Could not extract variable contribution data from the tree.\n")
}