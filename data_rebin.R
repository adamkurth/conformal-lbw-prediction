rm(list=ls())
library(rnn)
year <- 2021
data.pwd <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/birthweight_data"
file.path <- sprintf("%s/natality%dus-original.csv", data.pwd, year)

# rebin (1) vs. without with 2.5kg (2)
# data.rebin.pwd <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/birthweight_data/rebin/"
data.rebin.pwd <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/birthweight_data/rebin_without_2.5kg/"

natalitydata <- read.csv(file.path)
names(natalitydata)

# BOY: sex
# MARRIED: dmar
# BLACK: mrace15 (=2 means Black Only)
# OVER33: mager
# HIGH SCHOOL: meduc
# FULL PRENATAL: precare5
# SMOKER: cig_0
# BIRTH WEIGHT: dbwt (response variable)

# Select and preprocess features
dat <- natalitydata[,c('sex', 'dmar', 'mrace15', 'mager', 'meduc', 'precare5', 'cig_0', 'dbwt')]
dat$sex <- as.factor(ifelse(dat$sex=="M", 1, 0))
dat$dmar <- as.factor(ifelse(dat$dmar==1, 1, 0))
dat$mrace15 <- as.factor(ifelse(dat$mrace15==2, 1, 0))
dat$mager <- as.factor(ifelse(dat$mager > 33, 1, 0))
dat$meduc <- as.factor(ifelse(dat$meduc==3, 1, 0))
dat$precare5 <- as.factor(ifelse(dat$precare5==1, 1, 0))
dat$cig_0 <- as.factor(ifelse(dat$cig_0 > 0, 1, 0))

# get rid of factor: turn into numeric 
# dat$sex <- ifelse(dat$sex=="M", 1, 0)
# dat$dmar <- ifelse(dat$dmar==1, 1, 0)
# dat$mrace15 <- ifelse(dat$mrace15==2, 1, 0)
# dat$mager <- ifelse(dat$mager > 33, 1, 0)
# dat$meduc <- ifelse(dat$meduc==3, 1, 0)
# dat$precare5 <- ifelse(dat$precare5==1, 1, 0)
# dat$cig_0 <- ifelse(dat$cig_0 > 0, 1, 0)

p <- ncol(dat) - 1
num_node <- 2^p

# preprocessing
dat <- na.omit(dat)
X <- dat[-ncol(dat)] # remove the response variable
y <- log(dat$dbwt) # log-transform the birth weight

xnode <- apply(X[, 1:p], 1, function(a) bin2int(t(as.matrix(as.numeric(a), 1, p))) + 1)

# generate all possible combinations
feature.names <- colnames(X)
X.all <- expand.grid(replicate(n=7, list(0:1)))
colnames(X.all) <- feature.names

# convert to factors w/ correct levels
X.all[] <- lapply(X.all, function(col) factor(col, levels = c(0, 1)))

binary.matrix.all <- as.matrix(
  sapply(X.all, function(col) as.numeric(as.character(col)))
)
xnode.all <- apply(binary.matrix.all, 1, function(row) bin2int(matrix(row, nrow=1)) + 1)
X.all$xnode <- xnode.all

#aggregate statistics
agg.dbwt <- aggregate(dat$dbwt, list(xnode = xnode), mean)
agg.log.dbwt <- aggregate(y, list(xnode = xnode), mean)
count.observations <- aggregate(list(count = rep(1, nrow(dat))), list(xnode = xnode), sum)

# Merge with all combinations
final.data <- merge(X.all, agg.dbwt, by = "xnode", all.x = TRUE)
final.data <- merge(final.data, agg.log.dbwt, by = "xnode", all.x = TRUE)
final.data <- merge(final.data, count.observations, by = "xnode", all.x = TRUE)

# Handle missing values
#final.data$x[is.na(final.data$x)] <- 0
final.data$x.x[is.na(final.data$x.x)] <- 0
final.data$count[is.na(final.data$count)] <- 0

# Create final matrices
X.matrix <- final.data[, feature.names]
Y.matrix <- data.frame(
  # actual_weight = final.data$x.y, 
  # log_weight = final.data$x.x,
  count = final.data$count
  )

# ---- ONLY CHANGE: use 10% quantiles instead of 100g increments ----
# (A) Restrict to birthweights <= 2500g, then compute 10% quantiles
dat_lte2.5 <- dat$dbwt[dat$dbwt <= 2500]
probs <- seq(0, 1, by = 0.1)
cut_points <- quantile(dat_lte2.5, probs = probs)  
# (B) Create counts for each quantile-based bin from 0 up to 2.5 kg
counts <- list()
for(i in 1:(length(cut_points) - 1)) {
  counts_table <- tapply(
    (dat$dbwt > cut_points[i] & dat$dbwt <= cut_points[i+1]), 
    xnode, 
    sum
  )
  name <- paste0(
    "counts_",
    formatC(probs[i],   format="f", digits=2), "_",
    formatC(probs[i+1], format="f", digits=2),
    "_quantile_",
    round(cut_points[i]   / 1000, 2), "_",
    round(cut_points[i+1] / 1000, 2),
    "kg"
  )
  counts[[name]] <- counts_table
}

# (C) Create final interval for > 2.5 kg
count_above_2.5kg <- tapply(dat$dbwt > 2500, xnode, sum)


# toggle this to include or exclude > 2.5kg
# counts[["counts_above_2.5kg"]] <- count_above_2.5kg



# Combine all counts
counts.df <- do.call(cbind, counts)


# Save data
write.csv(final.data, sprintf("%s/final_data_%d.csv", data.rebin.pwd, year), row.names = FALSE)
write.csv(X.matrix, sprintf("%s/X_matrix_%d.csv", data.rebin.pwd, year), row.names = FALSE)
write.csv(Y.matrix, sprintf("%s/Y_matrix_%d.csv", data.rebin.pwd, year), row.names = FALSE)
write.csv(counts.df, sprintf("%s/counts_df_%d.csv", data.rebin.pwd, year), row.names = FALSE)

save(
  X.matrix, 
  Y.matrix, 
  counts.df, 
  final.data, 
  file = sprintf("%s/birthweight_data_rebin_%d.RData", data.rebin.pwd, year)
)
