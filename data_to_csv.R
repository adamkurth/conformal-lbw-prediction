library(rnn)

# 2021 Data

year <- 2021
data.pwd <- "/Users/adamkurth/Documents/RStudio/ms-thesis-kurth/birthweight_data"
#file_path <- sprintf("/Users/adamkurth/Documents/RStudio/ms-thesis-kurth/birthweight_data/natality%dus.csv", year)
file_path <- sprintf("%s/natality%dus.csv", data.pwd, year)
natalitydata <- read.csv(file_path)
names(natalitydata)


# BOY: sex
# MARRIED: dmar
# BLACK: mrace15 (=2 means Black Only)
# OVER33: mager
# HIGH SCHOOL: meduc
# FULL PRENATAL: precare5
# SMOKER: cig_0
# BIRTH WEIGHT: dbwt

dat <- natalitydata[,c('sex', 'dmar', 'mrace15', 'mager', 'meduc', 'precare5', 'cig_0', 'dbwt')]
# dat <- natalitydata[,c('mrace15', 'mager','cig_0', 'dbwt')]
dat$sex <- as.factor(ifelse(dat$sex=="M", 1, 0))
dat$dmar <- as.factor(ifelse(dat$dmar==1, 1, 0))
dat$mrace15 <- as.factor(ifelse(dat$mrace15==2, 1, 0))
dat$mager <- as.factor(ifelse(dat$mager > 33, 1, 0))
dat$meduc <- as.factor(ifelse(dat$meduc==3, 1, 0))
dat$precare5 <- as.factor(ifelse(dat$precare5==1, 1, 0))
dat$cig_0 <- as.factor(ifelse(dat$cig_0 > 0, 1, 0))

p <- ncol(dat) - 1
num_node <- 2^p


dat <- na.omit(dat)
X <- dat[-ncol(dat)] # remove the response variable
y <- log(dat$dbwt) # log-transform the birth weight

xnode <- apply(X[, 1:p], 1, function(a) bin2int(t(as.matrix(as.numeric(a), 1, p))) + 1)

# Sufficient statistics

n_in_node <- table(xnode)

var_y <- tapply(y, xnode, var)

sum_y <- tapply(y, xnode, sum)
sum_y_squared <- tapply(y^2, xnode, sum)

sum_y2_squared <- tapply(y^4, xnode, sum)

sum_y3 <- tapply(y^3, xnode, sum)
sum_y3_squared <- tapply(y^6, xnode, sum)


# cut_points <-c(seq(0, 800, by = 200), seq(900, 1500, by = 100))
cut_points <- seq(0, 2500, by = 250)

counts <- list()

for(i in 1:(length(cut_points)-1)){
  counts_table <- tapply((dat$dbwt > cut_points[i] & dat$dbwt <= cut_points[i+1]), xnode, sum)
  name <- paste("count_btw_", cut_points[i]/1000,"_",cut_points[i+1]/1000, "kg",sep = "")
  counts[[name]] <- counts_table
}

counts_df <- do.call(cbind, counts)

# count_btw_1.5_2.5kg <- tapply((dat$dbwt > 1400 & dat$dbwt < 2500), xnode, sum)
count_above_2.5kg <- tapply(dat$dbwt>2500, xnode, sum)


# Number of data points
n <- sum(n_in_node)

# Combine the results into a data frame
results <- data.frame(
  var_y = var_y,
  sum_y = sum_y,
  sum_y_squared = sum_y_squared,
  sum_y2_squared = sum_y2_squared,
  sum_y3 = sum_y3,
  sum_y3_squared = sum_y3_squared,
  n_in_node = as.vector(n_in_node)
  # count_above_2.5kg = count_above_2.5kg
  # count_btw_1.5_2.5kg = count_btw_1.5_2.5kg
  
)

# Write the data frame to a CSV file
file_path_out <- sprintf("%s/results%dus.csv", data.pwd, year)
write.csv(cbind(results, counts_df, count_above_2.5kg = count_above_2.5kg), file_path_out, row.names = FALSE)
