## Read data and sort by bagged probability
#-----------------------------------------------------------
year <- 2021

# Determine type based on data path
use_without_2_5kg <- FALSE  # Set to TRUE for type 2 (without 2.5kg), FALSE for type 1
type <- ifelse(use_without_2_5kg, 2, 1)

# Set appropriate data path based on type
if(type == 1) {
  # load WITH 2.5kg (below)
  data.pwd <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/birthweight_data/rebin/"
  results.pwd <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/results/results1"
  table.pwd <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/table/results1"
  load(file = sprintf("%s/bootstrap_tree_results_%d.RData", results.pwd, year))
} else {
  # load WITHOUT 2.5kg (below)
  data.pwd <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/birthweight_data/rebin_without_2.5kg/"
  results.pwd <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/results/results2"
  table.pwd <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/table/results2"
  load(file = sprintf("%s/bootstrap_tree_results_%d.RData", results.pwd, year))
}

# Convert factor variables to numeric if needed
if(is.factor(pred.results$sex)) {
  pred.results$sex <- as.numeric(as.character(pred.results$sex))
  pred.results$dmar <- as.numeric(as.character(pred.results$dmar))
  pred.results$mrace15 <- as.numeric(as.character(pred.results$mrace15))
  pred.results$mager <- as.numeric(as.character(pred.results$mager))
  pred.results$meduc <- as.numeric(as.character(pred.results$meduc))
  pred.results$precare5 <- as.numeric(as.character(pred.results$precare5))
  pred.results$cig_0 <- as.numeric(as.character(pred.results$cig_0))
}

preds.sorted <- pred.results[order(pred.results$bagged.pred), ]
lowest.5.pred <- head(preds.sorted, 5)
highest.5.pred <- tail(preds.sorted, 5)
highest.5.pred <- highest.5.pred[order(-highest.5.pred$bagged.pred), ]

## Function to find common predictors
find.common.preds <- function(data) {
  preds <- c("sex", "dmar", "mrace15", "mager", "meduc", "precare5", "cig_0")
  common <- c()
  
  for (pred in preds) {
    values <- unique(data[[pred]])
    if (length(values) == 1) {
      common <- c(common, pred)
    }
  }
  
  return(common)
}

## Find common predictors for lowest and highest groups
lowest.5.common <- find.common.preds(lowest.5.pred)
highest.5.common <- find.common.preds(highest.5.pred)

## Create the common predictor value strings
lowest.vals <- c(); highest.vals <- c()
for (pred in lowest.5.common) {
  lowest.vals <- c(lowest.vals, paste0(pred, " = ", lowest.5.pred[1, pred]))
}
for (pred in highest.5.common) {
  highest.vals <- c(highest.vals, paste0(pred, " = ", highest.5.pred[1, pred]))
}

## Create more compact publication-ready table with fixed LaTeX
## Create more compact publication-ready table with CI information
create.compact.table <- function() {
  # Convert binary values to readable descriptions
  predictor.labels <- c(
    "sex" = "Sex",
    "dmar" = "Marital Status",
    "mrace15" = "Race (Black)",
    "mager" = "Age > 33",
    "meduc" = "High School Ed.",
    "precare5" = "Full Prenatal",
    "cig_0" = "Smoker"
  )
  
  binary.meanings <- list(
    "sex" = c("Female", "Male"),
    "dmar" = c("Not Married", "Married"),
    "mrace15" = c("Not Black", "Black"),
    "mager" = c("≤ 33", "> 33"),
    "meduc" = c("No", "Yes"),
    "precare5" = c("No", "Yes"),
    "cig_0" = c("No", "Yes")
  )
  
  # Calculate average CI width for both groups
  mean.low.ci.width <- mean(lowest.5.pred$width)
  mean.high.ci.width <- mean(highest.5.pred$width)
  
  # Calculate average lower and upper bounds
  mean.low.lower <- mean(lowest.5.pred$lower.ci)
  mean.low.upper <- mean(lowest.5.pred$upper.ci)
  mean.high.lower <- mean(highest.5.pred$lower.ci)
  mean.high.upper <- mean(highest.5.pred$upper.ci)
  
  # Create rows for the table
  table.rows <- c()
  
  # Add header row for groups
  predictor.row <- c("Predictor", "Lowest 5", "Highest 5")
  table.rows <- c(table.rows, paste(predictor.row, collapse = " & "), "\\\\")
  table.rows <- c(table.rows, "\\midrule")
  
  # Add mean bagged probability
  mean.prob.row <- c("Mean Bagged Probability", 
                     sprintf("%.4f", mean(lowest.5.pred$bagged.pred)),
                     sprintf("%.4f", mean(highest.5.pred$bagged.pred)))
  table.rows <- c(table.rows, paste(mean.prob.row, collapse = " & "), "\\\\")
  
  # Add mean CI bounds
  mean.ci.row <- c("Mean 95\\% CI", 
                   sprintf("[%.4f, %.4f]", mean.low.lower, mean.low.upper),
                   sprintf("[%.4f, %.4f]", mean.high.lower, mean.high.upper))
  table.rows <- c(table.rows, paste(mean.ci.row, collapse = " & "), "\\\\")
  
  # Add mean CI width
  mean.width.row <- c("Mean CI Width", 
                      sprintf("%.4f", mean.low.ci.width),
                      sprintf("%.4f", mean.high.ci.width))
  table.rows <- c(table.rows, paste(mean.width.row, collapse = " & "), "\\\\")
  
  table.rows <- c(table.rows, "\\midrule")
  
  # Add rows for each predictor
  all.predictors <- c("sex", "dmar", "mrace15", "mager", "meduc", "precare5", "cig_0")
  
  for (pred in all.predictors) {
    # Get values for lowest and highest groups
    if (pred %in% lowest.5.common) {
      # Value is consistent across all classes in this group
      val <- lowest.5.pred[1, pred]
      lowest.val <- binary.meanings[[pred]][val + 1]
    } else {
      # Mixed values in this group
      lowest.val <- "\\textemdash"
    }
    
    if (pred %in% highest.5.common) {
      # Value is consistent across all classes in this group
      val <- highest.5.pred[1, pred]
      highest.val <- binary.meanings[[pred]][val + 1]
    } else {
      # Mixed values in this group
      highest.val <- "\\textemdash"
    }
    
    # Format row with bold for common predictors
    pred.label <- predictor.labels[pred]
    if (pred %in% lowest.5.common && pred %in% highest.5.common) {
      pred.label <- paste0("\\textbf{", pred.label, "}")
      lowest.val <- paste0("\\textbf{", lowest.val, "}")
      highest.val <- paste0("\\textbf{", highest.val, "}")
    }
    
    pred.row <- c(pred.label, lowest.val, highest.val)
    table.rows <- c(table.rows, paste(pred.row, collapse = " & "), "\\\\")
  }
  
  # Create LaTeX table with fixed structure
  latex.table <- c(
    "% Compact table comparing lowest and highest bagged probability classes",
    "% Add these packages to your LaTeX preamble:",
    "% \\usepackage{booktabs}",
    "% \\usepackage{siunitx}",
    "% \\usepackage{caption}",
    "% \\usepackage{threeparttablex}",
    "",
    "\\begin{table}[htbp]",
    "\\centering",
    "\\caption{Comparison of Common Predictors Between Lowest and Highest Bagged Probability Classes}",
    "\\label{tab:bagged_probability_comparison}",
    "\\begin{threeparttable}",
    "\\begin{tabular}{lcc}",
    "\\toprule",
    paste(table.rows, collapse = "\n"),
    "\\bottomrule",
    "\\end{tabular}",
    "\\begin{tablenotes}[flushleft]",
    "\\small",
    "\\item \\textit{Note}: Bold values indicate predictors that are consistent across all classes within the group.",
    "\\item ``\\textemdash'' indicates mixed values within the group.",
    "\\item CI = Confidence Interval for bagged probability estimates.",
    "\\end{tablenotes}",
    "\\end{threeparttable}",
    "\\end{table}"
  )
  
  return(latex.table)
}


## Generate the compact table and save to the proper directory
compact.table <- create.compact.table()

# Save to the appropriate directory
compact.table.file <- sprintf("%s/compact_bagged_probability_table_type%d.tex", table.pwd, type)
writeLines(compact.table, compact.table.file)

# Summary report to console
cat("Compact table saved to:", compact.table.file, "\n\n")
cat("Common predictors in lowest 5 classes:", paste(lowest.5.common, collapse=", "), "\n")
cat("Common predictors in highest 5 classes:", paste(highest.5.common, collapse=", "), "\n\n")

cat("Lowest 5 class statistics:\n")
print(lowest.5.pred[, c("node", "bagged.pred", "sex", "dmar", "mrace15", "mager", "meduc", "precare5", "cig_0")])

cat("\nHighest 5 class statistics:\n")
print(highest.5.pred[, c("node", "bagged.pred", "sex", "dmar", "mrace15", "mager", "meduc", "precare5", "cig_0")])