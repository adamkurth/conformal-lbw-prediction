library(rnn)

year <- 2021
data.pwd <- "/Users/adamkurth/Documents/RStudio/ms-thesis-kurth/birthweight_data"
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
# BIRTH WEIGHT: dbwt (response variable)

# Select and preprocess features
dat <- natalitydata[,c('sex', 'dmar', 'mrace15', 'mager', 'meduc', 'precare5', 'cig_0', 'dbwt')]

# dat <- natalitydata[,c('mrace15', 'mager','cig_0', 'dbwt')]
dat$sex <- as.factor(ifelse(dat$sex=="M", 1, 0))
dat$dmar <- as.factor(ifelse(dat$dmar==1, 1, 0))
dat$mrace15 <- as.factor(ifelse(dat$mrace15==2, 1, 0))
dat$mager <- as.factor(ifelse(dat$mager > 33, 1, 0))
dat$meduc <- as.factor(ifelse(dat$meduc==3, 1, 0))
dat$precare5 <- as.factor(ifelse(dat$precare5==1, 1, 0))
dat$cig_0 <- as.factor(ifelse(dat$cig_0 > 0, 1, 0))

# preprocessing
dat <- na.omit(dat)
X.raw <- dat[-ncol(dat)] # remove the response variable (birth weight)
y <- log(dat$dbwt) # log-transform the birth weight


# generate all possible combinations
feature.names <- c('sex', 'dmar', 'mrace15', 'mager', 'meduc', 'precare5', 'cig_0')
X.all <- expand.grid(replicate(n=7, list(0:1))) # all possible combinations of 0 and 1
colnames(X.all) <- feature.names

# convert to factors w/ correct levels
X.all[] <- lapply(X.all, function(col), factor(col, levels = c(0, 1)))

binary.matrix.all <- as.matrix(
  sapply(X.all, function(col) as.numeric(as.character(col))) # convert to numeric type matrix (0 and 1)
)
# 7 binary features, xnode converts combinations into unique integer identifier 
# xnode: integer representation of binary matrix 
xnode.all <- apply(binary.matrix.all, 1, function(row) bin2int(matrix(data=row, nrow=1)) + 1) # +1 to start from 1
X.all$xnode <- xnode.all

# p <- ncol(dat) - 1
# num_node <- 2^p
# xnode <- apply(X[, 1:p], 1, function(a) bin2int(t(as.matrix(as.numeric(a), 1, p))) + 1)

# Merge aggregated data with all combinations
agg.dbwt <- aggregate(dat$dbwt, list(xnode=xnode), mean)
agg.log.dbwt <- aggregate(y, list(xnode=xnode), mean)
count.observations <- aggregate(list(count=rep(1,nrow(dat))), list(xnode=xnode), sum)

# handle missing values 
final.data$x[is.na(final.data$x)] <- 0          # actual birthweight mean
final_data$x.x[is.na(final_data$x.x)] <- 0      # Log-transformed mean
final_data$count[is.na(final_data$count)] <- 0  # Observation count

# create final matrices 
X <- final.data[,feature.names]
Y.actual <- final.data$x
Y.log <- final.data$x.x
observation.counts <- final.data$count

valid.combinations <- observation.counts > 0 # only consider combinations with at least one observation
X <- X[valid.combinations,]
Y.actual <- Y.actual[valid.combinations]
Y.log <- Y.log[valid.combinations]

write.csv(X, sprintf("%s/X_matrix_%d.csv", data.pwd, year), row.names = FALSE)
write.csv(data.frame(actual_weight = Y_actual, log_weight = Y_log), 
          sprintf("%s/Y_values_%d.csv", data.pwd, year), row.names = FALSE)
