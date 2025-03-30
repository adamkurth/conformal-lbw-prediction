rm(list=ls())
library(xtable)
library(knitr)
year <- 2021

# Determine type based on data path
use_without_2_5kg <- FALSE  # Set to TRUE for type 2 (without 2.5kg), FALSE for type 1
type <- ifelse(use_without_2_5kg, 2, 1)

# Set appropriate data path based on type
if(type == 1) {
  data.pwd <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/birthweight_data/rebin/"
} else {
  data.pwd <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/birthweight_data/rebin_without_2.5kg/"
}

# Set output directory based on type
table.pwd <- sprintf("/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/table/results%d", type)
dir.create(table.pwd, showWarnings = FALSE, recursive = TRUE)

# Load data
load(sprintf("%s/birthweight_data_rebin_%d.RData", data.pwd, year))
load(sprintf("%s/informed_prior_%d.RData", data.pwd, 2020))
load(sprintf("%s/bootstrap_results_%d.RData", data.pwd, year))
load(sprintf("%s/dm_tree_%d.RData", data.pwd, year))

# Function to create formatted variable importance tables
create_importance_table <- function(tree, type, min_rows_for_two_columns = 20) {
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
  
  # Create the appropriate table - one or two columns based on the number of rows
  if(nrow(var_importance) >= min_rows_for_two_columns) {
    # Create two-column table
    cat("% Two-column variable importance table\n",
        file = sprintf("%s/dm_tree_var_imp_type%d.tex", table.pwd, type))
    
    cat("\\begingroup\n\\begin{table}[htbp]\n\\centering\n\\setlength{\\tabcolsep}{0.5em}\n\\renewcommand{\\arraystretch}{0.9}\n\\footnotesize\n",
        file = sprintf("%s/dm_tree_var_imp_type%d.tex", table.pwd, type),
        append = TRUE)
    
    cat(sprintf("\\caption{DM Tree Split Variables Importance Summary Ranked by Improvement (Type %d)}\n\\label{tab:var_imp_summary_type%d}\n", type, type),
        file = sprintf("%s/dm_tree_var_imp_type%d.tex", table.pwd, type),
        append = TRUE)
    
    # Calculate split point for two columns
    half_rows <- ceiling(nrow(var_importance) / 2)
    
    cat("\\begin{tabular}{lcr|lcr}\n\\hline\nVariable & Improvement & Rank & Variable & Improvement & Rank \\\\ \\hline\n",
        file = sprintf("%s/dm_tree_var_imp_type%d.tex", table.pwd, type),
        append = TRUE)
    
    # Write each row of the two-column table
    for(i in 1:half_rows) {
      second_idx <- i + half_rows
      if(second_idx <= nrow(var_importance)) {
        cat(sprintf("%s & %.2f & %d & %s & %.2f & %d \\\\\n", 
                    var_importance$variable[i], 
                    var_importance$total_improvement[i],
                    var_importance$rank[i],
                    var_importance$variable[second_idx], 
                    var_importance$total_improvement[second_idx],
                    var_importance$rank[second_idx]),
            file = sprintf("%s/dm_tree_var_imp_type%d.tex", table.pwd, type),
            append = TRUE)
      } else {
        cat(sprintf("%s & %.2f & %d & & & \\\\\n", 
                    var_importance$variable[i], 
                    var_importance$total_improvement[i],
                    var_importance$rank[i]),
            file = sprintf("%s/dm_tree_var_imp_type%d.tex", table.pwd, type),
            append = TRUE)
      }
    }
    
    cat("\\hline\n\\end{tabular}\n\\end{table}\n\\endgroup\n",
        file = sprintf("%s/dm_tree_var_imp_type%d.tex", table.pwd, type),
        append = TRUE)
    
    cat(sprintf("Created two-column importance table for Type %d with %d variables\n", 
                type, nrow(var_importance)))
  } else {
    # Create single-column table
    cat("% Single-column variable importance table\n",
        file = sprintf("%s/dm_tree_var_imp_type%d.tex", table.pwd, type))
    
    cat("\\begingroup\n\\begin{table}[htbp]\n\\centering\n\\setlength{\\tabcolsep}{0.5em}\n\\renewcommand{\\arraystretch}{0.9}\n\\footnotesize\n",
        file = sprintf("%s/dm_tree_var_imp_type%d.tex", table.pwd, type),
        append = TRUE)
    
    cat(sprintf("\\caption{DM Tree Split Variables Importance Summary Ranked by Improvement (Type %d)}\n\\label{tab:var_imp_summary_type%d}\n", type, type),
        file = sprintf("%s/dm_tree_var_imp_type%d.tex", table.pwd, type),
        append = TRUE)
    
    cat("\\begin{tabular}{lcr}\n\\hline\nVariable & Improvement & Rank \\\\ \\hline\n",
        file = sprintf("%s/dm_tree_var_imp_type%d.tex", table.pwd, type),
        append = TRUE)
    
    # Write each row of the single-column table
    for(i in 1:nrow(var_importance)) {
      cat(sprintf("%s & %.2f & %d \\\\\n", 
                  var_importance$variable[i], 
                  var_importance$total_improvement[i],
                  var_importance$rank[i]),
          file = sprintf("%s/dm_tree_var_imp_type%d.tex", table.pwd, type),
          append = TRUE)
    }
    
    cat("\\hline\n\\end{tabular}\n\\end{table}\n\\endgroup\n",
        file = sprintf("%s/dm_tree_var_imp_type%d.tex", table.pwd, type),
        append = TRUE)
    
    cat(sprintf("Created single-column importance table for Type %d with %d variables\n", 
                type, nrow(var_importance)))
  }
  
  # Also create a multipage longtable version for reference
  cat("\\begingroup\\small\\setlength{\\tabcolsep}{0.5em}\\renewcommand{\\arraystretch}{0.9}\n", 
      file = sprintf("%s/dm_tree_var_imp_type%d_multipage.tex", table.pwd, type))
  
  cat(sprintf("\\begin{longtable}{lcr}\n\\caption{DM Tree Split Variables Importance Summary Ranked by Improvement (Type %d)}\\\\\n\\label{tab:var_imp_summary_type%d_multipage}\\\\\n", type, type),
      file = sprintf("%s/dm_tree_var_imp_type%d_multipage.tex", table.pwd, type),
      append = TRUE)
  
  cat("\\hline\nVariable & Improvement & Rank \\\\ \\hline\n\\endhead\n",
      file = sprintf("%s/dm_tree_var_imp_type%d_multipage.tex", table.pwd, type),
      append = TRUE)
  
  # Write each row
  for(i in 1:nrow(var_importance)) {
    cat(sprintf("%s & %.2f & %d \\\\\n", 
                var_importance$variable[i], 
                var_importance$total_improvement[i],
                var_importance$rank[i]),
        file = sprintf("%s/dm_tree_var_imp_type%d_multipage.tex", table.pwd, type),
        append = TRUE)
  }
  
  cat("\\hline\n\\end{longtable}\n\\endgroup\n",
      file = sprintf("%s/dm_tree_var_imp_type%d_multipage.tex", table.pwd, type),
      append = TRUE)
  
  return(var_importance)
}

# Function to create just the simplified splits table
create_simplified_splits_table <- function(tree, type) {
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
  
  # Create simplified splits data (removing Categories and Adjustment columns)
  simplified_splits <- data.frame(
    Count = splits_data$count,
    Improvement = splits_data$improve,
    Index = splits_data$index,
    Variable = splits_data$variable,
    stringsAsFactors = FALSE
  )
  
  # Create multipage longtable version for simplified splits table
  cat("\\begingroup\\small\\setlength{\\tabcolsep}{0.6em}\\renewcommand{\\arraystretch}{0.9}\n", 
      file = sprintf("%s/dm_tree_splits_type%d_multipage.tex", table.pwd, type))
  
  cat(sprintf("\\begin{longtable}{rrrp{3cm}}\n\\caption{DM Tree Split Variables Ranked by Improvement (Type %d)}\\\\\n\\label{tab:splits_type%d_multipage}\\\\\n", type, type),
      file = sprintf("%s/dm_tree_splits_type%d_multipage.tex", table.pwd, type),
      append = TRUE)
  
  cat("\\hline\nCount & Improvement & Index & Variable \\\\ \\hline\n\\endhead\n",
      file = sprintf("%s/dm_tree_splits_type%d_multipage.tex", table.pwd, type),
      append = TRUE)
  
  # Write each row
  for(i in 1:nrow(simplified_splits)) {
    cat(sprintf("%d & %.2f & %d & %s \\\\\n", 
                simplified_splits$Count[i], 
                simplified_splits$Improvement[i],
                simplified_splits$Index[i],
                simplified_splits$Variable[i]),
        file = sprintf("%s/dm_tree_splits_type%d_multipage.tex", table.pwd, type),
        append = TRUE)
  }
  
  cat("\\hline\n\\end{longtable}\n\\endgroup\n",
      file = sprintf("%s/dm_tree_splits_type%d_multipage.tex", table.pwd, type),
      append = TRUE)
  
  cat(sprintf("Created simplified splits table for Type %d with %d rows\n", 
              type, nrow(simplified_splits)))
  
  return(simplified_splits)
}

# Create the tables for the current type
var_importance <- create_importance_table(dm.tree, type)
splits_table <- create_simplified_splits_table(dm.tree, type)

# Print summary
cat("\nCreated tables for Type", type, "in", table.pwd, "\n")
cat("\nVariable importance summary (top rows):\n")
print(head(var_importance, 5))

# Now switch to the other type and generate tables for it too
use_without_2_5kg <- !use_without_2_5kg
type <- ifelse(use_without_2_5kg, 2, 1)

# Set appropriate data path based on type
if(type == 1) {
  data.pwd <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/birthweight_data/rebin/"
} else {
  data.pwd <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/birthweight_data/rebin_without_2.5kg/"
}

# Set output directory based on type
table.pwd <- sprintf("/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/table/results%d", type)
dir.create(table.pwd, showWarnings = FALSE, recursive = TRUE)

# Load data for the other type
load(sprintf("%s/birthweight_data_rebin_%d.RData", data.pwd, year))
load(sprintf("%s/informed_prior_%d.RData", data.pwd, 2020))
load(sprintf("%s/bootstrap_results_%d.RData", data.pwd, year))
load(sprintf("%s/dm_tree_%d.RData", data.pwd, year))

# Create tables for the other type too
var_importance_other <- create_importance_table(dm.tree, type)
splits_table_other <- create_simplified_splits_table(dm.tree, type)

# Print summary for the other type
cat("\nCreated tables for Type", type, "in", table.pwd, "\n")
cat("\nVariable importance summary (top rows):\n")
print(head(var_importance_other, 5))

cat("\nAll tables have been created successfully!\n")