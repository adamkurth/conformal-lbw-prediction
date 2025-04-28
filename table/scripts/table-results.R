rm(list=ls())
library(rpart)
library(rpart.plot)
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

# Load data
load(sprintf("%s/birthweight_data_rebin_%d.RData", data.pwd, year))
load(sprintf("%s/informed_prior_%d.RData", data.pwd, 2020))
load(sprintf("%s/bootstrap_results_%d.RData", data.pwd, year))
load(sprintf("%s/dm_tree_%d.RData", data.pwd, year))

clean.table <- function(tree){
  tree.frame <- tree$frame
  split.info <- attr(tree, "splits")
  node.nums <- as.numeric(rownames(tree.frame))
  is.term <- tree.frame$var == "<leaf>"
  
  results.table <- data.frame(
    Node = node.nums,
    Terminal = ifelse(is.term, "Yes", "No"),
    Split = ifelse(is.term, NA, as.character(tree.frame$var)),
    n = tree.frame$n,
    Deviance = tree.frame$dev,
    stringsAsFactors = FALSE
  )
  
  results.table <- results.table[order(results.table$Node),]
  results.table$SplitCondition <- NA
  
  # Simpler approach to get split conditions - focus on direct split details
  if(any(!is.term)) {
    for (i in which(!is.term)){
      var.name <- as.character(tree.frame$var[i])
      
      # Create split condition based on variable type
      if(var.name %in% c("mrace15", "cig_0", "dmar", "mager", "sex", "precare5", "meduc")) {
        # For binary/categorical variables, use "= 0" or "= 1" format
        results.table$SplitCondition[i] <- paste0(var.name, " = 0")
      } else {
        # For numeric variables, handle with actual split value if available
        if(!is.null(split.info) && nrow(split.info) > 0) {
          # Try to extract the split value if possible
          if(as.character(i) %in% rownames(split.info)) {
            split.val <- split.info[as.character(i), "index"]
            if(!is.na(split.val)) {
              results.table$SplitCondition[i] <- paste0(var.name, " < ", round(split.val, 3))
            } else {
              results.table$SplitCondition[i] <- var.name
            }
          } else {
            results.table$SplitCondition[i] <- var.name
          }
        } else {
          results.table$SplitCondition[i] <- var.name
        }
      }
    }
  }
  
  return(results.table)
}

# Create the table
tree.table <- clean.table(dm.tree)

# Output the formatted table
library(knitr)
kable(tree.table, caption = sprintf("DM Tree Results using RPART (Type %d)", type), 
      format = "markdown", row.names = FALSE)

# Save as CSV
save.path.csv <- sprintf("%s/dm_tree_results_type%d.csv", table.pwd, type)
write.csv(tree.table, save.path.csv, row.names = FALSE)

# For LaTeX output
if(requireNamespace("xtable", quietly = TRUE)) {
  library(xtable)
  xtable.output <- xtable(tree.table, 
                          caption = sprintf("DM Tree Results using RPART (Type %d)", type),
                          label = sprintf("tab:rpart_results_type%d", type))
  
  # Create the save path for LaTeX file
  save.path.tex <- sprintf("%s/dm_tree_results_type%d.tex", table.pwd, type)
  
  # OPTION 1: Smaller font size, compact table
  print(xtable.output, 
        file = save.path.tex, 
        include.rownames = FALSE,
        size = "footnotesize",  # Reduces font size
        tabular.environment = "tabular",
        floating = TRUE,
        table.placement = "htbp",
        sanitize.text.function = function(x) x,
        add.to.row = list(pos = list(-1),
                          command = "\\setlength{\\tabcolsep}{3pt}\n")  # Reduces column spacing
  )
  
  # OPTION 2: Split across multiple pages using longtable
  save.path.tex.long <- sprintf("%s/dm_tree_results_type%d_multipage.tex", table.pwd, type)
  print(xtable.output,
        file = save.path.tex.long,
        include.rownames = FALSE,
        tabular.environment = "longtable",  # For multi-page tables
        floating = FALSE,
        size = "small",
        sanitize.text.function = function(x) x,
        add.to.row = list(pos = list(-1),
                          command = paste0(sprintf("\\caption{DM Tree Results using RPART (Type %d)}", type),
                                           sprintf("\\label{tab:rpart_results_type%d}", type),
                                           "\\\\\n"))
  )
  
  # OPTION 3: Landscape orientation
  save.path.tex.landscape <- sprintf("%s/dm_tree_results_type%d_landscape.tex", table.pwd, type)
  cat("\\begin{landscape}\n", file = save.path.tex.landscape)
  print(xtable.output,
        file = save.path.tex.landscape,
        include.rownames = FALSE,
        append = TRUE,
        size = "normalsize",  # Normal size since we have more space in landscape
        tabular.environment = "tabular",
        floating = TRUE,
        table.placement = "htbp",
        sanitize.text.function = function(x) x
  )
  cat("\\end{landscape}\n", file = save.path.tex.landscape, append = TRUE)
}