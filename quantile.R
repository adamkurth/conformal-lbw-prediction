rm(list = ls())

# ----- Step 1: Load the Quantile Cutpoints for 2020 -----
data.rebin.pwd <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/birthweight_data/rebin/"
load(sprintf("%s/quantile_cutpoints_%d.RData", data.rebin.pwd, 2020))
# 'cut_points' now holds the quantile boundaries (e.g., 227, 1170, 1644, ..., 2500)

# ----- Step 2: Load the 2020 Natality Data -----
data.pwd <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/birthweight_data"
file.path <- sprintf("%s/natality%dus-original.csv", data.pwd, 2020)
natalitydata <- read.csv(file.path)

# Use only birthweights <= 2500 g (as in your quantile computation)
dat_lte2.5 <- natalitydata$dbwt[natalitydata$dbwt <= 2500]

# ----- Step 3: Bin the Birthweights Using the Cutpoints -----
# Each birthweight is assigned to a quantile interval defined by cut_points
birth_bins <- cut(
  dat_lte2.5,
  breaks = cut_points,
  include.lowest = TRUE,
  right = TRUE,
  labels = FALSE
)

# ----- Step 4: Count and Normalize the Bin Frequencies -----
# Count the number of observations in each interval (there are length(cut_points)-1 intervals)
bin_counts <- table(factor(birth_bins, levels = 1:(length(cut_points) - 1)))

# Create the informed prior by normalizing these counts to sum to 1
alphavec <- as.numeric(bin_counts) / sum(bin_counts)

# Optionally, scale the prior by a pseudo-count weight (alpha0)
alpha0 <- 1  # Change this if you need a stronger/weaker prior effect
alphavec <- alpha0 * alphavec

# Print the informed prior vector
print(alphavec)

# ----- Step 5: Save the Informed Prior for Future Use -----
save(alphavec, file = sprintf("%s/informed_prior_%d.RData", data.rebin.pwd, 2020))