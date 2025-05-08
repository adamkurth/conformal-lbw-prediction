# quantile_table.R - Generate the exact LaTeX table format as requested
# Paths to the data
data_pwd_type1 <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/birthweight_data/rebin"
data_pwd_type2 <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/birthweight_data/rebin_without_2.5kg"
output_dir <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/results"
year <- 2020

# Load the data for both types
load(sprintf("%s/quantile_cutpoints_%d.RData", data_pwd_type1, year))
cut_points_type1 <- cut.points

load(sprintf("%s/quantile_cutpoints_%d.RData", data_pwd_type2, year))
cut_points_type2 <- cut.points

load(sprintf("%s/informed_prior_%d.RData", data_pwd_type1, year))
alphavec_type1 <- alphavec

load(sprintf("%s/informed_prior_%d.RData", data_pwd_type2, year))
alphavec_type2 <- alphavec

# Check lengths
cat("Length of alphavec_type1:", length(alphavec_type1), "\n")
cat("Length of alphavec_type2:", length(alphavec_type2), "\n")

# Normalize priors to percentages
normalized_prior_type1 <- alphavec_type1/sum(alphavec_type1) * 100
normalized_prior_type2 <- alphavec_type2/sum(alphavec_type2) * 100

# Instead of creating a data frame, we'll directly create the LaTeX table
# with the exact format requested

# Create the exact table string
latex_table <- c(
  "% insert table of quantile cutpoints",
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Birth-weight quantile cut points and Dirichlet priors}",
  "\\label{tab:birthweight_quantiles}",
  "% siunitx is needed only for the S columns",
  "\\begin{tabular}{@{}l c S[table-format=2.2] c S[table-format=2.2]@{}}",
  "\\toprule",
  "& \\multicolumn{2}{c}{\\textbf{Type 1: LBW + Normal}} &",
  "  \\multicolumn{2}{c}{\\textbf{Type 2: LBW only}} \\\\",
  "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}",
  "\\textbf{Quantile} &",
  "  \\textbf{Range (g)} & {\\textbf{Prior (\\%)}} &",
  "  \\textbf{Range (g)} & {\\textbf{Prior (\\%)}} \\\\",
  "\\midrule"
)

# Add the quantile rows
for (i in 1:10) {
  # Format the ranges with double dashes
  range_type1 <- paste0(round(cut_points_type1[i]), "--", round(cut_points_type1[i+1]))
  range_type2 <- paste0(round(cut_points_type2[i]), "--", round(cut_points_type2[i+1]))
  
  # Format the prior percentages (already normalized to percentages above)
  prior_type1 <- sprintf("%.2f", normalized_prior_type1[i])
  prior_type2 <- sprintf("%.2f", normalized_prior_type2[i])
  
  # Add the row
  row <- paste0(
    "Q", i, "     & ", 
    range_type1, " & ", 
    prior_type1, " & ", 
    range_type2, " & ", 
    prior_type2, " \\\\"
  )
  
  latex_table <- c(latex_table, row)
}

# Add the Normal row (specific format with \textgreater and --)
normal_prior <- sprintf("%.2f", normalized_prior_type1[11])
normal_row <- paste0(
  "Normal & \\textgreater{}2500 & ", 
  normal_prior, " & -- & -- \\\\"
)
latex_table <- c(latex_table, normal_row)

# Add the bottom rule
latex_table <- c(latex_table, "\\bottomrule", "\\end{tabular}", "\\end{table}")

# Combine into a single string
latex_table_str <- paste(latex_table, collapse = "\n")

# Print to console
cat("Generated LaTeX Table:\n\n")
cat(latex_table_str)
cat("\n\n")

# Save to file
output_file <- file.path(output_dir, "birthweight_quantiles_table.tex")
writeLines(latex_table_str, output_file)

cat("Table saved to:", output_file, "\n")