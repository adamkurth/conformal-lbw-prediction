# Load necessary packages
library(xtable)
library(dplyr)

# Determine type based on data path
use_without_2_5kg <- FALSE  # Set to TRUE for type 2 (without 2.5kg), FALSE for type 1
type <- ifelse(use_without_2_5kg, 2, 1)

# Set appropriate data path based on type
if(type == 1) {
  results.pwd <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/results/results1"
  save.path <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/results/results1/varimp"
} else {
  save.path <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/results/results2/varimp"
  results.pwd <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/results/results2"
}
dir.create(save.path, showWarnings = FALSE)

# Define B (number of bootstrap samples)
B <- 10000

# Data extraction
pred.results <- bootstrap.tree.results$pred.results
var.freq.df <- bootstrap.tree.results$var.freq.df
top.var.df <- bootstrap.tree.results$top.var.df
n.vars.summary <- bootstrap.tree.results$n.vars.summary

# Helper function for creating publication-quality LaTeX tables
create_thesis_table <- function(data, caption, label, digits = NULL, 
                                include.rownames = FALSE, include.colnames = TRUE,
                                font_size = "normalsize") {
  
  # Create alignment vector with correct length (ncol + 1)
  ncols <- ncol(data)
  align_vector <- c("l", rep("c", ncols))
  
  # Create properly sized digits vector if specified
  if (!is.null(digits)) {
    if (length(digits) == 1) {
      # If a single value, repeat it for all columns
      digits_vector <- rep(digits, ncols + 1)
    } else if (length(digits) != ncols + 1) {
      # If length doesn't match, warn and use default
      warning("digits vector length should be ncol(data) + 1. Using default.")
      digits_vector <- NULL
    } else {
      # Use as provided
      digits_vector <- digits
    }
  } else {
    digits_vector <- NULL
  }
  
  # Create xtable object with appropriate digits
  xtable_obj <- xtable(data, caption = caption, label = label, 
                       align = align_vector, digits = digits_vector)
  
  # Convert to character string with enhanced formatting
  table_content <- capture.output(
    print(xtable_obj, 
          include.rownames = include.rownames,
          include.colnames = include.colnames,
          booktabs = TRUE,
          floating = TRUE,
          tabular.environment = "tabular",
          hline.after = NULL,
          sanitize.text.function = function(x) x,
          type = "latex")
  )
  
  # Join lines into a single string
  table_text <- paste(table_content, collapse = "\n")
  
  # Format with booktabs - using fixed=TRUE to avoid regex interpretation
  table_text <- gsub("\\hline", "\\toprule", table_text, fixed = TRUE)
  
  # Add midrule after headers
  if(include.colnames) {
    # We need a literal pattern for the line ending after column headers
    lines <- strsplit(table_text, "\n")[[1]]
    header_end_idx <- which(grepl("\\\\\\\\", lines))[1]
    
    if(!is.na(header_end_idx)) {
      lines[header_end_idx] <- paste0(lines[header_end_idx], "\n\\midrule")
      table_text <- paste(lines, collapse = "\n")
    }
  }
  
  # Replace last hline with bottomrule - using string manipulation to avoid regex issues
  table_text <- sub("\\hline\n\\end", "\\bottomrule\n\\end", table_text, fixed = TRUE)
  
  # Add font size command if specified (without using regex)
  if (font_size != "normalsize") {
    # Insert size command after begin{table}
    lines <- strsplit(table_text, "\n")[[1]]
    for (i in 1:length(lines)) {
      if (grepl("begin{table}", lines[i], fixed = TRUE)) {
        lines[i+1] <- paste0("\\", font_size, "\n", lines[i+1])
        break
      }
    }
    table_text <- paste(lines, collapse = "\n")
  }
  
  return(table_text)
}

# Function to add footnote to table
add_footnote <- function(table_text, footnote) {
  # Find where the tabular environment ends
  end_pos <- regexpr("\\end{tabular}", table_text, fixed = TRUE)
  
  if(end_pos > 0) {
    # Split the string
    first_part <- substr(table_text, 1, end_pos + 12) # 12 is length of "\\end{tabular}"
    last_part <- substr(table_text, end_pos + 13, nchar(table_text))
    
    # Insert footnote
    result <- paste0(first_part, "\n", footnote, last_part)
    return(result)
  } else {
    warning("Could not find end of tabular environment")
    return(table_text)
  }
}

# Function to write LaTeX tables to files
write_latex_table <- function(content, filename) {
  file_path <- file.path(save.path, filename)
  
  tryCatch({
    # Make sure the output directory exists
    dir.create(dirname(file_path), showWarnings = FALSE, recursive = TRUE)
    
    # Write to file
    writeLines(content, file_path)
    cat(sprintf("Saved LaTeX table to: %s\n", file_path))
  }, error = function(e) {
    cat("Error writing to file:", e$message, "\n")
    cat("Attempting alternative method...\n")
    
    # Alternative approach with explicit connection
    con <- file(file_path, "w")
    writeLines(content, con)
    close(con)
    cat(sprintf("Saved LaTeX table to: %s\n", file_path))
  })
}

# 1. Variable Frequency Table
# ---------------------------
var_freq_table <- var.freq.df %>%
  arrange(desc(frequency))

var_freq_caption <- "Variable Usage Frequency in Bootstrap Trees"
var_freq_label <- "tab:var_freq"

var_freq_latex <- create_thesis_table(
  var_freq_table, 
  var_freq_caption, 
  var_freq_label,
  digits = 4,
  font_size = "normalsize"
)

# Add footnote
var_freq_footnote <- "\\caption*{\\textit{Note:} Frequency indicates the proportion of bootstrap trees that include each variable.}"
var_freq_latex <- add_footnote(var_freq_latex, var_freq_footnote)

write_latex_table(var_freq_latex, "variable_frequency_table.tex")

# 2. Top Split Variable Table
# ---------------------------
top_var_table <- top.var.df %>%
  arrange(desc(frequency))

top_var_caption <- "Frequency of Variables as Top Split in Bootstrap Trees"
top_var_label <- "tab:top_var"

top_var_latex <- create_thesis_table(
  top_var_table, 
  top_var_caption, 
  top_var_label,
  digits = 4,
  font_size = "normalsize"
)

# Add footnote
top_var_footnote <- "\\caption*{\\textit{Note:} Frequency indicates proportion of bootstrap trees using the variable as the first (root) split.}"
top_var_latex <- add_footnote(top_var_latex, top_var_footnote)

write_latex_table(top_var_latex, "top_split_table.tex")

# 3. Variable Count Summary
# -------------------------
var_count_caption <- "Summary of Number of Variables Used in Bootstrap Trees"
var_count_label <- "tab:var_count"

# Create data frame with just one row for summary statistics
var_count_table <- data.frame(
  Mean = n.vars.summary$mean,
  Median = n.vars.summary$median,
  Min = n.vars.summary$min,
  Max = n.vars.summary$max,
  SD = n.vars.summary$sd
)

# Figure out the correct digits vector length for this table
var_count_digits <- c(0, 4, 0, 0, 0, 6)

var_count_latex <- create_thesis_table(
  var_count_table, 
  var_count_caption, 
  var_count_label,
  include.rownames = FALSE,
  digits = var_count_digits,
  font_size = "normalsize"
)

# Add footnote
var_count_footnote <- sprintf("\\caption*{\\textit{Note:} Statistics describing the number of unique variables used across %d bootstrap trees.}", B)
var_count_latex <- add_footnote(var_count_latex, var_count_footnote)

write_latex_table(var_count_latex, "variable_count_table.tex")

# 4. Prediction Results Table (improved for thesis)
# ---------------------------------------------------
# Select subset of columns to display (removing the predictor columns)
pred_subset <- pred.results %>%
  select(node, bagged.pred, oob.pred, se, lower.ci, upper.ci, width) %>%
  head(10)

pred_caption <- "Prediction Results from Bootstrap Trees"
pred_label <- "tab:pred_results"

# Better column names for publication
colnames(pred_subset) <- c("Node", "Bagged Pred", "OOB Pred", "SE", 
                           "Lower CI", "Upper CI", "Width")

# Count the number of columns in pred_subset and create appropriate digits vector
pred_cols <- ncol(pred_subset)
# Create appropriate digits vector for each column
pred_digits <- c(0, 0, 2, 2, 3, 2, 2, 3)

pred_latex <- create_thesis_table(
  pred_subset, 
  pred_caption, 
  pred_label,
  digits = pred_digits,
  font_size = "footnotesize"  # Smaller font to fit all columns
)

# Add more descriptive footnote
pred_footnote <- "\\caption*{\\textit{Note:} First 10 rows of prediction results shown. Bagged Pred = bagged predictions, OOB Pred = out-of-bag predictions, SE = standard error, Lower/Upper CI = 95\\% confidence interval bounds, Width = confidence interval width.}"
pred_latex <- add_footnote(pred_latex, pred_footnote)

write_latex_table(pred_latex, "prediction_results_table.tex")

cat("Successfully created publication-quality LaTeX tables in", save.path, "\n")



# using var.across.depths
combined.content <- c(
  "% Table 1: Variable Usage Across Different Tree Depths",
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Variable Importance and Usage Analysis from Bootstrap Trees}",
  "\\label{tab:var-across-depths}",
  "\\begin{tabular}{lrrrr}",
  "\\hline",
  "Variable & Depth 2 & Depth 3 & Depth 4 & Depth 5 & Total Usage \\",
  "\\hline"
)

for ( i in 1:nrow(var.across.depths)){

    var <- var.across.depths$variable[i]
    # want integer format
    depth2 <- sprintf("%d", var.across.depths$depth_2[i])
    depth3 <- sprintf("%d", var.across.depths$depth_3[i])
    depth4 <- sprintf("%d", var.across.depths$depth_4[i])
    depth5 <- sprintf("%d", var.across.depths$depth_5[i])
    total.usage <- sprintf("%d", var.across.depths$total.usage[i])

    combined.content <- c(combined.content,
                         paste0(var, " & ", depth2, " & ", depth3, " & ", depth4, " & ", depth5, " & ", total.usage, " \\ "))
}

combined.content <- c(
  combined.content,
  "\\hline",
  "\\end{tabular}",
  paste0("\\caption*{\\textit{Note:}}"),
  "\\end{table}")

write.latex.table(combined.content, "var_across_depths.tex")
# 
# 
# #-----------------------------------------------------------
# # 1) Combined Metrics Table
# #-----------------------------------------------------------
# combined_content <- c(
#   "% Table 1: Combined metrics table",
#   "\\begin{table}[htbp]",
#   "\\centering",
#   "\\caption{Variable Importance and Usage Analysis from Bootstrap Trees}",
#   "\\label{tab:combined_metrics}",
#   "\\begin{tabular}{lrrrr}",
#   "\\hline",
#   "Variable & Frequency & Top Split & Mean Deviance & Importance \\\\",
#   "& (in Tree) & Frequency & Reduction & Score \\\\ ",
#   "\\hline"
# )
# 
# # Add rows with improved formatting
# for(i in 1:nrow(var.metrics)) {
#   var <- var.metrics$variable[i]
#   freq <- sprintf("%.4f", var.metrics$frequency[i])
# 
#   # Find top split frequency if available
#   top_freq <- "0.0000"
#   for(j in 1:nrow(top.var.df)) {
#     if(top.var.df$variable[j] == var) {
#       top_freq <- sprintf("%.4f", top.var.df$frequency[j])
#       break
#     }
#   }
# 
#   # Format deviance reduction with commas as thousands separator and 1 decimal place
#   # Use mean.dev.reduction or mean_dev_reduction depending on your column name
#   if("mean.dev.reduction" %in% colnames(var.metrics)) {
#     dev_red <- formatC(var.metrics$mean.dev.reduction[i], format = "f", digits = 1, big.mark = ",")
#   } else if("mean_dev_reduction" %in% colnames(var.metrics)) {
#     dev_red <- formatC(var.metrics$mean_dev_reduction[i], format = "f", digits = 1, big.mark = ",")
#   } else {
#     stop("Could not find deviance reduction column in var.metrics")
#   }
# 
#   # Format importance score with 4 decimal places and a + sign for positive values
#   # Use importance.score or importance_score depending on your column name
#   if("importance.score" %in% colnames(var.metrics)) {
#     imp_score <- sprintf("%.4f", var.metrics$importance.score[i])
#     if(var.metrics$importance.score[i] > 0) {
#       imp_score <- paste0("+", imp_score)
#     }
#   } else if("importance_score" %in% colnames(var.metrics)) {
#     imp_score <- sprintf("%.4f", var.metrics$importance_score[i])
#     if(var.metrics$importance_score[i] > 0) {
#       imp_score <- paste0("+", imp_score)
#     }
#   } else {
#     stop("Could not find importance score column in var.metrics")
#   }
# 
#   # Add the formatted row
#   combined_content <- c(combined_content,
#                         paste0(var, " & ", freq, " & ", top_freq, " & ", dev_red, " & ", imp_score, " \\\\ "))
# }
# 
# combined_content <- c(
#   combined_content,
#   "\\hline",
#   "\\end{tabular}",
#   paste0("\\caption*{\\textit{Note:} Analysis based on ", B, " bootstrap samples. 'Frequency (in Tree)' shows how often a variable appears anywhere in the tree. 'Top Split Frequency' shows how often it's chosen as the first split. 'Importance Score' is the normalized mean deviance reduction.}"),
#   "\\end{table}"
# )
# 
# # Save to file
# write_latex_table(combined_content, "variable_metrics_table.tex")
# 
# #-----------------------------------------------------------
# # 2) Top Split Table
# #-----------------------------------------------------------
# top_split_content <- c(
#   "% Table 2: Top split frequencies",
#   "\\begin{table}[htbp]",
#   "\\centering",
#   "\\caption{Frequency of Variables as Top Split in Bootstrap Trees}",
#   "\\label{tab:top_split}",
#   "\\begin{tabular}{lr}",
#   "\\hline",
#   "Variable & Frequency \\\\ ",
#   "\\hline"
# )
# 
# # Add rows with consistent formatting
# for(i in 1:nrow(top.var.df)) {
#   var <- top.var.df$variable[i]
#   freq <- sprintf("%.4f", top.var.df$frequency[i])
#   top_split_content <- c(top_split_content, paste0(var, " & ", freq, " \\\\ "))
# }
# 
# top_split_content <- c(
#   top_split_content,
#   "\\hline",
#   "\\end{tabular}",
#   "\\caption*{\\textit{Note:} Frequency indicates proportion of bootstrap trees using the variable as the first (root) split.}",
#   "\\end{table}"
# )
# 
# # Save to file
# write_latex_table(top_split_content, "top_split_table.tex")
# 
# #-----------------------------------------------------------
# # 3) Variable Count Summary
# #-----------------------------------------------------------
# var_count_content <- c(
#   "% Table 3: Variable count summary",
#   "\\begin{table}[htbp]",
#   "\\centering",
#   "\\caption{Summary of Number of Variables Used in Bootstrap Trees}",
#   "\\label{tab:var_count}",
#   "\\begin{tabular}{ccccc}",
#   "\\hline",
#   "Mean & Median & Min & Max & SD \\\\ ",
#   "\\hline"
# )
# 
# # Format numbers consistently, handle potential column name differences
# # Some R data frames use '.' and others use '_' in column names
# if("mean" %in% colnames(n.vars.summary)) {
#   mean_val <- sprintf("%.4f", n.vars.summary$mean)
# } else {
#   mean_val <- sprintf("%.4f", n.vars.summary$mean[1])
# }
# 
# if("median" %in% colnames(n.vars.summary)) {
#   median_val <- sprintf("%d", n.vars.summary$median)
# } else {
#   median_val <- sprintf("%d", n.vars.summary$median[1])
# }
# 
# if("min" %in% colnames(n.vars.summary)) {
#   min_val <- sprintf("%d", n.vars.summary$min)
# } else {
#   min_val <- sprintf("%d", n.vars.summary$min[1])
# }
# 
# if("max" %in% colnames(n.vars.summary)) {
#   max_val <- sprintf("%d", n.vars.summary$max)
# } else {
#   max_val <- sprintf("%d", n.vars.summary$max[1])
# }
# 
# if("sd" %in% colnames(n.vars.summary)) {
#   sd_val <- sprintf("%.6f", n.vars.summary$sd)
# } else {
#   sd_val <- sprintf("%.6f", n.vars.summary$sd[1])
# }
# 
# var_count_content <- c(
#   var_count_content,
#   paste0(mean_val, " & ", median_val, " & ", min_val, " & ", max_val, " & ", sd_val, " \\\\ "),
#   "\\hline",
#   "\\end{tabular}",
#   paste0("\\caption*{\\textit{Note:} Statistics describing the number of unique variables used across ", B, " bootstrap trees.}"),
#   "\\end{table}"
# )
# 
# # Save to file
# write_latex_table(var_count_content, "variable_count_table.tex")
# 
# cat("Successfully created all three LaTeX tables in", save.path, "\n")