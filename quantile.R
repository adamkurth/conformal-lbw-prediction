rm(list = ls())

# ----- Step 1: Load the Quantile Cutpoints for 2020 -----
data.pwd <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/birthweight_data"
# data.rebin.pwd <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/birthweight_data/rebin"
data.rebin.pwd <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/birthweight_data/rebin_without_2.5kg"

# ----- Step 2: Load the 2020 Natality Data -----
file.path <- sprintf("%s/natality%dus-original.csv", data.pwd, 2020)
natalitydata <- read.csv(file.path)

# Use only birthweights <= 2500 g (as in your quantile computation)
dat_lte2.5 <- natalitydata$dbwt[natalitydata$dbwt <= 2500]

# ----- Step 3: calculate cutpoints -----
# num.quantiles <- 11
num.quantiles <- 10

# We add the min and max boundaries to ensure complete coverage
cut.points <- c(
  min(dat_lte2.5),
  quantile(dat_lte2.5, probs = seq(1/num.quantiles, (num.quantiles-1)/num.quantiles, 1/num.quantiles)),
  max(dat_lte2.5)
)
save(cut.points, file = sprintf("%s/quantile_cutpoints_%d.RData", data.rebin.pwd, 2020))

# ----- Step 4: Bin the Birthweights Using the Cutpoints -----
birth.bins <- cut(
  dat_lte2.5,
  breaks = cut.points,
  include.lowest = TRUE,
  right = TRUE,
  labels = FALSE
)

# ----- Step 4: Count and Normalize the Bin Frequencies -----
# Count the number of observations in each interval (there are length(cut.points)-1 intervals)
bin.counts <- table(factor(birth.bins, levels = 1:(length(cut.points) - 1)))

# Create the informed prior by normalizing these counts to sum to 1
alphavec <- as.numeric(bin.counts) / sum(bin.counts)

# Optionally, scale the prior by a pseudo-count weight (alpha0)
alpha0 <- 1  # Change this if you need a stronger/weaker prior effect
alphavec <- alpha0 * alphavec

# Print the informed prior vector
print(alphavec)

# ----- Step 5: Save the Informed Prior for Future Use -----
save(alphavec, file = sprintf("%s/informed_prior_%d.RData", data.rebin.pwd, 2020))

cat("\nQuantile Cutpoints:\n")
print(cut.points)

cat("\nBin Counts:\n")
print(bin.counts)

cat("\nProportions per bin:\n")
print(as.numeric(bin.counts) / sum(bin.counts))
