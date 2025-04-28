# table-boot-freq.R

## This script generates a LaTeX table comparing variable importance across two models 
## over all bootstrap samples: 1) full model, 2) LBW-only model.

# Load required libraries
library(dplyr)

# Define paths for both models
year <- 2021
results1.pwd <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/results/results1"
results2.pwd <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/results/results2"
table.pwd <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/table"

# Load data for Type 1 (Full Model)
load(sprintf("%s/bootstrap_tree_results_%d.RData", results1.pwd, year))
var_freq_type1 <- bootstrap.tree.results$var.freq.df
top_var_type1 <- bootstrap.tree.results$top.var.df

# Load data for Type 2 (LBW-Only Model)
load(sprintf("%s/bootstrap_tree_results_%d.RData", results2.pwd, year))
var_freq_type2 <- bootstrap.tree.results$var.freq.df
top_var_type2 <- bootstrap.tree.results$top.var.df

# Process the top split variable data
# Only keep top 3 variables from each model
top_vars_combined <- c(
  top_var_type1$variable, 
  top_var_type2$variable
) %>% 
  unique()

# Initialize dataframe for top variables
top_vars_df <- data.frame(
  variable = top_vars_combined,
  full_model = 0,
  lbw_only_model = 0
)

# Fill in the values
for (var in top_vars_combined) {
  # Full model
  idx <- which(top_var_type1$variable == var)
  if (length(idx) > 0) {
    top_vars_df$full_model[top_vars_df$variable == var] <- top_var_type1$frequency[idx]
  }
  
  # LBW-only model
  idx <- which(top_var_type2$variable == var)
  if (length(idx) > 0) {
    top_vars_df$lbw_only_model[top_vars_df$variable == var] <- top_var_type2$frequency[idx]
  }
}

# Sort by full_model frequency and take top 3
top_vars_df <- top_vars_df %>%
  arrange(desc(full_model)) %>%
  head(3)

# Process variable frequency data
# Keep all variables of interest
freq_vars_of_interest <- c("sex", "dmar", "mrace15", "mager", "precare5", "cig_0", "meduc")

freq_vars_df <- data.frame(
  variable = freq_vars_of_interest,
  full_model = 0,
  lbw_only_model = 0
)

# Fill in the values
for (var in freq_vars_of_interest) {
  # Full model
  idx <- which(var_freq_type1$variable == var)
  if (length(idx) > 0) {
    freq_vars_df$full_model[freq_vars_df$variable == var] <- var_freq_type1$frequency[idx]
  }
  
  # LBW-only model
  idx <- which(var_freq_type2$variable == var)
  if (length(idx) > 0) {
    freq_vars_df$lbw_only_model[freq_vars_df$variable == var] <- var_freq_type2$frequency[idx]
  }
}

# Generate LaTeX table
latex_table <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Comparison of Variable Importance in Full Model and LBW-Only Model}",
  "\\label{tab:var_comparison}",
  "\\begin{tabular}{lcc}",
  "\\toprule",
  "& \\multicolumn{1}{c}{Full Model} & \\multicolumn{1}{c}{LBW-Only Model} \\\\",
  "\\cmidrule(lr){2-2} \\cmidrule(lr){3-3}",
  "\\multicolumn{3}{l}{\\textit{Initial Split Variable}} \\\\",
  "\\midrule"
)

# Add top variables rows
for (i in 1:nrow(top_vars_df)) {
  latex_table <- c(
    latex_table,
    sprintf("%s & %.4f & %.4f \\\\", 
            top_vars_df$variable[i], 
            top_vars_df$full_model[i], 
            top_vars_df$lbw_only_model[i])
  )
}

latex_table <- c(
  latex_table,
  "\\midrule",
  "\\multicolumn{3}{l}{\\textit{Variable Frequency}} \\\\",
  "\\midrule"
)

# Add frequency rows
for (i in 1:nrow(freq_vars_df)) {
  latex_table <- c(
    latex_table,
    sprintf("%s & %.4f & %.4f \\\\", 
            freq_vars_df$variable[i], 
            freq_vars_df$full_model[i], 
            freq_vars_df$lbw_only_model[i])
  )
}

latex_table <- c(
  latex_table,
  "\\bottomrule",
  "\\multicolumn{3}{p{0.8\\textwidth}}{\\textit{Note:} The table shows variable frequency (proportion of trees containing each variable) and initial split variable (normalized measure of predictive contribution) for both models.} \\\\",
  "\\end{tabular}",
  "\\end{table}"
)

# Write the LaTeX table to a file
output_file <- file.path(table.pwd, "var_importance_comparison.tex")
writeLines(latex_table, output_file)

cat(sprintf("LaTeX table saved to: %s\n", output_file))
cat("Table generation complete!\n")
