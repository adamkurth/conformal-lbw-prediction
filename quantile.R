rm(list = ls())
# ----- Step 1: Set up paths and determine type -----
use.without.2.5kg <- TRUE  # Set to TRUE for type 2 (without 2.5kg), FALSE for type 1
type <- ifelse(use.without.2.5kg, 2, 1)

data.pwd <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/birthweight_data"
if(type == 1) {
  # load WITH 2.5kg
    data.rebin.pwd <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/birthweight_data/rebin"
    plots.pwd <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/plots/type1"
} else {
    # load WITHOUT 2.5kg
    data.rebin.pwd <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/birthweight_data/rebin_without_2.5kg"
    plots.pwd <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/plots/type2"
}

# ----- Step 2: Load the 2020 Natality Data -----
file.path <- sprintf("%s/natality%dus-original.csv", data.pwd, 2020)
natalitydata <- read.csv(file.path)

# ----- Step 3: Calculate LBW and normal birth weight proportions -----
# Get all birth weights
all.weights <- natalitydata$dbwt
# Count low birth weights (≤ 2500g) and normal birth weights (> 2500g)
count.lbw <- sum(all.weights <= 2500)
count.normal <- sum(all.weights > 2500)
total.count <- length(all.weights)

# Check the proportions
prop.lbw <- count.lbw / total.count
prop.normal <- count.normal / total.count

cat("Proportion of LBW:", prop.lbw, "\n")
cat("Proportion of normal birth weights:", prop.normal, "\n")

# ----- Step 4: Calculate quantiles for LBW data -----
dat.lte2.5 <- all.weights[all.weights <= 2500]
num.quantiles <- 10 # DO NOT CHANGE

# We add the min and max boundaries to ensure complete coverage
cut.points <- c(
  min(dat.lte2.5),
  quantile(dat.lte2.5, probs = seq(1/num.quantiles, (num.quantiles-1)/num.quantiles, 1/num.quantiles)),
  max(dat.lte2.5)
)

# ----- Step 5: Bin the LBW data using the cutpoints -----
birth.bins <- cut(
  dat.lte2.5,
  breaks = cut.points,
  include.lowest = TRUE,
  right = TRUE,
  labels = FALSE
)

# ----- Step 6: Count and normalize to create informed prior -----
# Count observations in each LBW interval
bin.counts <- table(factor(birth.bins, levels = 1:(length(cut.points) - 1)))

# Create informed prior based on actual data distribution
if(type == 1) {
      # Type 1: Include normal birth weights as 11th category
      # calculate the proportion within LBW categories
  
      lbw.props <- as.numeric(bin.counts) / sum(bin.counts)
      
      # Scale the LBW proportions to account for the normal birth weight category
      # ensures that when we add the normal birth weight proportion,all proportions will sum to 1
      alphavec.lbw <- lbw.props * prop.lbw
      
      # Add the normal birth weight proportion as the last category
      alphavec <- c(alphavec.lbw, prop.normal)
      
      # category labels for plotting
      category.labels <- c(paste0("Q", 1:10), "Normal")
      
} else if(type == 2) {
  
      # Type 2: Exclude normal birth weights, use only 10 categories
      # Just create uniform probabilities for the LBW-only case
      alphavec <- rep(1/num.quantiles, num.quantiles)
      
      # Still calculate the cut points for binning future data
      cut.points <- c(
        min(dat.lte2.5),
        quantile(dat.lte2.5, probs = seq(1/num.quantiles, (num.quantiles-1)/num.quantiles, 1/num.quantiles)),
        max(dat.lte2.5)
      )
      
      # Bin the data (only needed for checking, not for setting alphavec)
      birth.bins <- cut(
        dat.lte2.5,
        breaks = cut.points,
        include.lowest = TRUE,
        right = TRUE,
        labels = FALSE
      )
      
      # Count observations (only for reporting)
      bin.counts <- table(factor(birth.bins, levels = 1:(length(cut.points) - 1)))
      
      # Create category labels for plotting
      category.labels <- paste0("Q", 1:10)
}

alpha0 <- 1
alphavec <- alpha0 * alphavec

cat("\nInformed Prior Vector (alphavec):\n")
print(alphavec)

# ----- Step 7: Save the data -----
save(cut.points, file = sprintf("%s/quantile_cutpoints_%d.RData", data.rebin.pwd, 2020))
save(alphavec, file = sprintf("%s/informed_prior_%d.RData", data.rebin.pwd, 2020))

# ----- Step 8: Print summary statistics -----
cat("\nQuantile Cutpoints (g):\n")
print(cut.points)

cat("\nBin Counts (LBW only):\n")
print(bin.counts)

cat("\nProportions per bin (distribution):\n")
print(alphavec / sum(alphavec))

# Create a bar plot to visualize the prior
pdf(sprintf("%s/alphavec_plot_%d.pdf", plots.pwd, 2020), width = 10, height = 8)
if(type == 1){
  main.label <- "Informed Dirichlet Prior (alphavec) with Normal Birth Weight"
} else {
  main.label <- "Informed Dirichlet Prior (alphavec) LBW-only"
}
barplot(alphavec, names.arg=category.labels, 
        col = "skyblue",
        border = "blue",
        main=main.label,
        xlab="Birth Weight Category", 
        ylab="Prior Probabability",
        ylim = c(0, max(alphavec) * 1.2))
dev.off()

# barplot(alphavec,
#         names.arg = paste("Q", 1:length(alphavec)),
#         col = "skyblue",
#         border = "blue",
#         main = "Informed Dirichlet Prior (alphavec)",
#         xlab = "Birth Weight Categories",
#         ylab = "Prior Probability",
#         ylim = c(0, max(alphavec) * 1.2))
# dev.off()
