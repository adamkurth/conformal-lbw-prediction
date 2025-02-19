rm(list=ls())
library(rnn)
year <- 2021
data.pwd <- "/Users/adamkurth/Documents/RStudio/ms-thesis-kurth/birthweight_data/"
data.rebin.pwd <- "/Users/adamkurth/Documents/RStudio/ms-thesis-kurth/birthweight_data/rebin/"
file.path <- sprintf("%s/natality%dus-original.csv", data.pwd, year)
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
Y.matrix <- data.frame(actual_weight = final.data$x.y, 
                       log_weight = final.data$x.x,
                       count = final.data$count)

# Keep original cutoff points and counts code
# Change the cut_points sequence to 100g (0.10kg) increments
cut_points <- seq(0, 2500, by = 100)  # From 0g to 2500g in 100g steps

# Update the counts generation loop
counts <- list()
for(i in 1:(length(cut_points)-1)){
  counts_table <- tapply((dat$dbwt > cut_points[i] & dat$dbwt <= cut_points[i+1]), xnode, sum)
  name <- paste("count_btw_", 
                cut_points[i]/1000, "_",  # Convert grams to kg
                cut_points[i+1]/1000, 
                "kg", sep = "")
  counts[[name]] <- counts_table
}

# Keep the final "above 2.5kg" category unchanged
count_above_2.5kg <- tapply(dat$dbwt > 2500, xnode, sum)
counts.df <- do.call(cbind, counts)
counts.df <- cbind(counts.df, count_above_2.5kg)

# Save data
write.csv(final.data, sprintf("%s/final_data_%d.csv", data.rebin.pwd, year), row.names = FALSE)
write.csv(X.matrix, sprintf("%s/X_matrix_%d.csv", data.pwd, year), row.names = FALSE)
write.csv(Y.matrix, sprintf("%s/Y_matrix_%d.csv", data.pwd, year), row.names = FALSE)
write.csv(counts.df, sprintf("%s/counts_df_%d.csv", data.pwd, year), row.names = FALSE)

save(X.matrix, Y.matrix, counts.df, final.data, file = sprintf("%s/birthweight_data_%d.RData", data.pwd, year))
