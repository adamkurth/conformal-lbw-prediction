rm(list=ls())
library(rpart)
library(rpart.plot)

#-----------------------------------------------------------
# A. Load Data & Construct Matrices
#-----------------------------------------------------------
year <- 2021

# Determine type based on data path
use_without_2_5kg <- FALSE  # Set to TRUE for type 2 (without 2.5kg), FALSE for type 1
type <- ifelse(use_without_2_5kg, 2, 1)

# Set appropriate data path based on type
if(type == 1) {
  # load WITH 2.5kg (below)
  data.pwd <- sprintf("/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/birthweight_data/rebin/")
  results.pwd <- sprintf("/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/results/results1")
  boot.pwd <- sprintf("/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/results/results1/plots")
  load(sprintf("%s/informed_prior_%d.RData", data.pwd, 2020))
  load(sprintf("%s/birthweight_data_rebin_%d.RData", data.pwd, year))
} else {
  # load WITHOUT 2.5kg (below)
  data.pwd <- sprintf("/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/birthweight_data/rebin_without_2.5kg/")
  results.pwd <- sprintf("/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/results/results2")
  boot.pwd <- sprintf("/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/results/results2/plots")
  load(sprintf("%s/informed_prior_%d.RData", data.pwd, 2020))
  load(sprintf("%s/birthweight_data_rebin_%d.RData", data.pwd, year))
}

# load('birthweight_data_2021.Rdata')
# load(sprintf("%s/birthweight_data_rebin_%d.RData", data.pwd, year))

print(alphavec)
length(alphavec)

Y.df <- counts.df
response.cols <- colnames(Y.df)

# prior
# informed prior based on 2020 data
# alphavec <- load(sprintf("%s/informed_prior_%d.RData", data.pwd, 2020))

# alphavec <- 1*rep(1/ncol(Y.df),ncol(Y.df)) # uniform
# alphavec <- 1*colSums(Y.df)/sum(Y.df) # informed

#-----------------------------------------------------------
# B. Simplified Dirichlet-Multinomial Functions
#-----------------------------------------------------------
log.dm.likelihood <- function(counts, alpha = 1) {
 # cat("** counts: ", counts, "\n")
  # 'counts': an integer vector (n_1, ..., n_K)
  # 'alpha':  scalar Dirichlet hyperparameter (assuming alpha_k = alpha for all k)
  #
  # Returns the log-likelihood of the Dirichlet-Multinomial model:
  #   log p(x | alpha) = log Gamma(alpha_0) - log Gamma(N+alpha_0)
  #                     + sum_k [ log Gamma(n_k + alpha) - log Gamma(alpha) ]
  # where alpha_0 = K * alpha,  N = sum(counts), and K = length(counts).
  
  N <- sum(counts)
  K <- length(counts)
  #alpha <- alpha*rep(1,K)
  
  # cat("N: ", N, "\n")
  # cat("K: ", K, "\n")

  # sum of alpha over all categories
  alpha_0 <- sum(alpha)

  # Term 1: log Gamma(alpha_0) + log Gamma(N+1) - log Gamma(N + alpha_0)
  term1 <- lgamma(alpha_0) + 0*lgamma(N + 1) - lgamma(N + alpha_0)
  
  # Term 2: sum over k of [log Gamma(n_k + alpha) - log Gamma(alpha) - log Gamma(n_k + 1)]
  term2 <- sum(lgamma(counts + alpha) - lgamma(alpha) - 0*lgamma(counts + 1))
  
  # Total log-likelihood
  ll <- term1 + term2

  # if numeric issues occur
  if (is.na(ll) || is.infinite(ll)) {
    ll <- -Inf
  }

  return(ll)
}

#-----------------------------------------------------------
# C. Define Custom rpart Method (DM with log-likelihood deviance)
#-----------------------------------------------------------
# rpart requires a list with the following named functions:
#    init, eval, split
#
# We'll define alpha=1 for each category (simple uniform prior).
# The 'eval' function must return a $deviance scalar.
# Here, deviance = - log_dm_likelihood( colSums(Y) ) -- i.e. negative LL.
# 
#    - We ignore weights (wt)
#    - We store in "label" something descriptive (i.e. colsum)
#    - We use the "summary" function to print the deviance in a readable way.

#----------------- 1) init() ------------------#
myinit <- function(y, offset, parms=NULL, wt=NULL) {
  # y: matrix of response counts 
  # offset: not used
  # parms: list containing alpha parameter
  # wt: not used
  # cat("** myinit() called. Y dimension:", dim(y), "\n")

  list(
    method="dm",
    y = y,
    parms = list(alpha=1),
    numresp = ncol(y),
    numy = ncol(y),
    summary = function(yval, dev, wt, ylevel, digits) {
      cat("DM Model\n")
      cat("Deviance: ", format(dev,digits=digits), "\n")
      cat("Counts: ", paste(format(yval, digits=digits), collapse=", "), "\n")
    },
    print = function(yval, dev, wt, ylevel, digits) {
      cat("Deviance:", format(dev, digits = digits), "\n")
    },
    text = function(yval, dev, wt, ylevel, digits, n, use.n){
      # yval is the colsum from myeval(); a vector of category counts.
      # keep labels short.
      total.counts <- sum(yval)
      lbl <- paste0(format(which.max(total.counts), digits=4))
      if (use.n) lbl <- paste(lbl, "\nn", n)
      return(lbl)
    }
  )
}

#----------------- 2) eval() ------------------#
myeval <- function(y, wt=NULL, parms=1) {
  # y: response matrix subset for current node
  # wt: weights (not used here)
  if (nrow(y)==1) {counts = y}else{
    counts <- colSums(y)}
  # dev <- dm.deviance(counts,alpha=1) # return -log.dm.likelihood
  dev <- -log.dm.likelihood(counts, alpha=alphavec) # counts is a vector!
  # 'label' can be something to print for the node. We'll just store colSums  
  return(list(label=counts, deviance=dev))
}

mysplit <- function(y, wt, x, parms, continuous = FALSE) {
  
  if (nrow(y)==1) {counts = y}else{
    counts <- colSums(y)}
  
  
  parent.dev <- -log.dm.likelihood(counts, alpha=alphavec)
  
  # Suppose 'x' is a factor or a binary 0/1 variable
  ux <- sort(unique(x)) # unique values of x
  goodness <- numeric(length(ux) - 1)# store improvement in deviance
  direction <- ux
  
  for(i in 1:(length(ux) - 1)) {
    split.val <- ux[i]
    left.idx <- x == split.val
    right.idx <- !left.idx
    
    if (sum(left.idx)>1){
         left.counts  <- colSums(y[left.idx,  , drop=FALSE])}else{left.counts <- y[left.idx, , drop = FALSE]}
         
         
         if (sum(right.idx)>1){
           right.counts  <- colSums(y[right.idx,  , drop=FALSE])}else{right.counts <- y[right.idx, , drop = FALSE]}
         
    
    
    child.dev <- -log.dm.likelihood(left.counts,alpha=alphavec) - log.dm.likelihood(right.counts,alpha=alphavec)

    goodness[i] <- parent.dev - child.dev
  }
  
  # 'goodness' is a vector of length = #possible splits. We have only 1 here.
  # 'direction' is a vector of length = #possible splits. Typically -1 means "x < cut" go left
  #   For a categorical, you can store an integer code. We'll just do -1.
  return(list(goodness = pmax(goodness, 0), direction = direction))

}

dm.method <- list(init=myinit, eval=myeval, split=mysplit, method="dm")

#-----------------------------------------------------------
# D. Fit rpart Tree Using Our Custom DM Method
#-----------------------------------------------------------
cat("\n--- Fitting the DM-based rpart tree ---\n")

# We set various control parameters to reduce complexity:
#   - usesurrogate=0, maxsurrogate=0 => no surrogate splits
#   - maxcompete=0 => do not evaluate competing splits
#   - xval=0 => no cross-validation
#   - cp=0 => allow splits with minimal improvement
#   - minsplit=5 => each node must have at least 5 obs
#   - maxdepth=5 => limit tree depth to 5

dm.control <- rpart.control(minsplit=2, cp=0, maxdepth=8, xval=0, usesurrogate = 0)
dm.tree <- rpart(Y.df~.,
  data = data.frame(X.matrix),
  method = dm.method,
  control = dm.control
)

#-----------------------------------------------------------
# E. Save and Inspect Results
#-----------------------------------------------------------
save(dm.tree, file=sprintf("%s/dm_tree_rebin_%d.RData", data.pwd, year))
print(dm.tree)
printcp(dm.tree)
plot(dm.tree,main="DM-based rpart Tree", uniform=TRUE, compress=TRUE)
text(dm.tree, use.n=TRUE, cex=1,font=3)

#-----------------------------------------------------------
# F. Multinomial Parametric Bootstrap Sampling
#-----------------------------------------------------------
# - calculate the multinomial probabilities for each cell using observed counts with Dirichlet prior
# - perform parametric bootstrap sampling, generate new counts based on multinomial dist/prob
# - fit new DM tree model on bootstrap samples
# - collect predictions/ci

B <- 10000  # Number of bootstrap samples
n.rows <- nrow(X.matrix)
lbw.cols <- 1:10  # indices of LBW categories (1:10, if 11th is above_2.5kg)
alpha.sig <- 0.05  # Significance level for confidence intervals

# Initialize matrices for predictions
Yhat.oob <- matrix(NA, nrow = n.rows, ncol = B)  # Out-of-bag predictions
Yhat.all <- matrix(0, nrow = n.rows, ncol = B)   # Full dataset predictions

# For tracking progress
cat("Starting multinomial parametric bootstrap with", B, "samples...\n")
pb <- txtProgressBar(min = 0, max = B, style = 3)

# First, calculate the multinomial probabilities for each node/cell
# These will be used to generate bootstrap samples
multinomial.probs <- matrix(0, nrow = n.rows, ncol = ncol(Y.df))
for(i in 1:n.rows) {
    # get counts for this cell 
    cell.counts <- as.numeric(Y.df[i,])
    # calculate prob w/ Dirichlet smoothing
    cell.probs <- (cell.counts + alphavec) / (sum(cell.counts) + sum(alphavec))
    multinomial.probs[i, ] <- cell.probs
}
  
  
# Loop over the bootstrap samples
for (b in 1:B) {
    # For parametric bootstrap, we generate new counts based on the multinomial model
    bootstrap.counts <- matrix(0, nrow = n.rows, ncol = ncol(Y.df))
    
    for (i in 1:n.rows) {
        # get original total count for this cell
        # Total counts in the original data for this cell 
        # i.e. sum of all categories in that row
        total.count <- sum(Y.df[i, ])
        if(total.count > 0){
          # generate new counts from multinomial dist using est. probs
          new.counts <- rmultinom(1, size=total.count, prob=multinomial.probs[i, ])
          bootstrap.counts[i,] <- new.counts
      }
  }
  # convert to dataframe for rpart
  bootstrap.counts.df <- as.data.frame(bootstrap.counts)
  colnames(bootstrap.counts.df) <- colnames(Y.df)
  
  # sample idx for bootstrap evaluation (mimic oob evaluation)
  idx <- sample(seq_len(n.rows), size = n.rows, replace = TRUE) # sample with replacement
  oob.idx <- setdiff(seq_len(n.rows), unique(idx)) # identify out-of-bag indices
  
  #fit model to parametric bootstrap sample
  dm.tree.b <- rpart(
    Y.df[idx, ] ~ .,
    data = data.frame(X.matrix[idx, , drop =FALSE]),
    method = dm.method,
    control = dm.control
  )
  
  # predictions on oob samples
  if(length(oob.idx > 0)) {
    
    preds.mat.oob <- predict(dm.tree.b, 
                             newdata=data.frame(X.matrix[oob.idx, , drop=FALSE]), 
                             type="matrix")
    # calculate LBW prob (sum of first 10 categories)
    # normalize by rowsums to get valid prob
    lbw.prob.oob <- rowSums(preds.mat.oob[, lbw.cols]) / rowSums(preds.mat.oob) # Nx1 matrix of LBW probabilities
    Yhat.oob[oob.idx, b] <- lbw.prob.oob
  }
  
  preds.mat.all <- predict(dm.tree.b, 
                       newdata=data.frame(X.matrix), 
                       type="matrix") # make predictions on full dataset
  
  # calculate LBW prob
  lbw.prob.all <- rowSums(preds.mat.all[, lbw.cols]) / rowSums(preds.mat.all) # Nx1 matrix of LBW probabilities
  Yhat.all[, b] <- lbw.prob.all # store the predictions for the full dataset
  
  setTxtProgressBar(pb, b) # Update the progress bar
  
}
close(pb)
cat("\n--- Multinomial bootstrap sampling completed ---\n")

#-----------------------------------------------------------
# G. Compute Predictions and Confidence Uncertainty Intervals
#-----------------------------------------------------------
# bagged predictions (across all bootstrap models)
bagged.preds <- rowMeans(Yhat.all)

# oob predictions (for evaluations)
oob.preds <- apply(Yhat.oob, 1, function(x) mean(x, na.rm = TRUE)) # mean of each row (for each observation)

# calculate std err 
pred.se <- apply(Yhat.all, 1, function(x) sd(x))

# confidence intervals 
lower.ci <- apply(Yhat.all, 1, function(x) quantile(x, alpha.sig / 2))
upper.ci <- apply(Yhat.all, 1, function(x) quantile(x, 1 - alpha.sig / 2))
interval.width <- upper.ci - lower.ci

#-----------------------------------------------------------
# H. Create Prediction Results Dataset
#-----------------------------------------------------------
pred.results <- data.frame(
  node = 1:n.rows,
  bagged.pred = bagged.preds,
  oob.pred = oob.preds,
  se = pred.se,
  lower.ci = lower.ci,
  upper.ci = upper.ci,
  width = interval.width
)

pred.results <- cbind(pred.results, X.matrix)

#-------------------------------------------------------------
# I. Visualize Bootstrap Results
#-------------------------------------------------------------
library(ggplot2)
library(tidyr)
library(dplyr)

# 0. Internal consistency check
coverage.check <- matrix(FALSE,nrow=n.rows,ncol=B)
for(b in 1:B) {
    # for each bootstrap sample, check if predictions are within CI bounds
    coverage.check[, b] <- (Yhat.all[, b] >= lower.ci) & (Yhat.all[, b] <= upper.ci)
}
empirical.coverage <- rowMeans(coverage.check)
pred.results$empirical.coverage <- empirical.coverage

# visualize coverage check
ggplot(pred.results, aes(x = empirical.coverage)) +
  geom_histogram(bins = 30, fill = "darkblue", color = "white") +
  geom_vline(xintercept = 0.95, color = "red", linetype = "dashed") +
  theme_minimal() +
  labs(title = "Empirical Coverage Rates of 95% Confidence Intervals",
       x = "Coverage Rate", 
       y = "Count",
       caption = "Red line shows nominal 95% coverage") +
  xlim(0.90, 1.0)
ggsave(sprintf("%s/empirical_coverage_%d.png", boot.pwd, year), width = 8, height = 6)

# 1. calibration check (bagged vs. oob) 
calibration.data <- data.frame(
  bagged = pred.results$bagged.pred,
  oob = pred.results$oob.pred
)
calibration.data <- na.omit(calibration.data)
ggplot(calibration.data, aes(x = bagged, y = oob)) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  geom_smooth(method = "loess", color = "blue") +
  theme_minimal() +
  labs(title = "Calibration Plot: Bagged vs. OOB Predictions",
       x = "Bagged Predictions", 
       y = "Out-of-Bag Predictions",
       caption = "Points close to diagonal indicate good calibration") +
  coord_equal()
ggsave(sprintf("%s/calibration_plot_%d.png", boot.pwd, year), width = 8, height = 8)


# 1. prediction intervals
ggplot(pred.results, aes(x = node, y = bagged.pred)) +
  geom_point() +
  geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), alpha = 0.2, fill = "blue") +
  geom_errorbar(aes(ymin = lower.ci, ymax = upper.ci), width = 0.2, alpha = 0.5) +
  theme_minimal() +
  labs(title = "Bootstrap Prediction Intervals for Low Birth Weight Probability",
       x = "Node Index", 
       y = "Predicted LBW Probability",
       caption = sprintf("Based on %d bootstrap samples", B)) +
  theme(plot.title = element_text(size = 14, face = "bold"))
ggsave(sprintf("%s/pred_interval_%d.png", boot.pwd, year), height = 10, width = 12)


# 2. uncertainty vs. prediction values
ggplot(pred.results, aes(x = bagged.pred, y = width)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess") +
  theme_minimal() +
  labs(title = "Uncertainty vs. Predicted Probability",
       x = "Predicted LBW Probability", 
       y = "Width of 95% Prediction Interval")
ggsave(sprintf("%s/uncertainty_vs_prediction_%d.png", boot.pwd, year), width = 8, height = 6)


# 3. dist of predictions
ggplot(pred.results, aes(x=bagged.pred)) + 
  geom_histogram(bins = 50, fill = "steelblue", color="white") + 
  theme_minimal() +
  labs(title = "Distribution of Bagged Predictions for LBW Probability", 
       x = "Predicted LBW Probability", 
       y = "Count") +
  theme(plot.title = element_text(size = 14, face = "bold"))
ggsave(sprintf("%s/pred_dist_%d.png", boot.pwd, year), height = 8, width = 6)


# 4. coefficient stability 
nodes.to.plot <- c(1,10,50,100) 
nodes.to.plot <- nodes.to.plot[nodes.to.plot <= n.rows]
node.dists <- data.frame(matrix(NA,nrow=length(nodes.to.plot) * B, ncol=3))
colnames(node.dists) <- c("node", "sample", "prediction")
# extract prediction distribution
row.idx <- 1
for (i in 1:length(nodes.to.plot)) {
  node <- nodes.to.plot[i]
  for (b in 1:B) {
    node.dists[row.idx, ] <- c(node, b, Yhat.all[node, b])
    row.idx <- row.idx + 1
  }
}

ggplot(node.dists, aes(x = prediction, fill = factor(node))) +
  geom_density(alpha = 0.5) +
  theme_minimal() +
  labs(title = "Bootstrap Distribution of Predictions for Selected Nodes",
       x = "Predicted LBW Probability", 
       y = "Density",
       fill = "Node") +
  scale_fill_brewer(palette = "Set1")
ggsave(sprintf("%s/bootstrap_dist_nodes_%d.png", boot.pwd, year), height=10,width=15)


#-----------------------------------------------------------
# K. Save Results
#-----------------------------------------------------------
# Save primary prediction results
save(pred.results, file = sprintf("%s/bootstrap_pred_results_%d.RData", results.pwd, year))

# Save raw bootstrap matrices for future analysis
save(Yhat.all, file = sprintf("%s/Yhat_all_%d.RData", results.pwd, year))
save(Yhat.oob, file = sprintf("%s/Yhat_oob_%d.RData", results.pwd, year))

# Save aggregated predictions
save(bagged.preds, oob.preds, lower.ci, upper.ci, 
     file = sprintf("%s/bagged_oob_preds_%d.RData", results.pwd, year))

# Print summary information
cat("\nBootstrap Analysis Summary:\n")
cat("------------------------\n")
cat("Number of bootstrap samples:", B, "\n")
cat("Average prediction interval width:", round(mean(interval.width), 4), "\n")
cat("Min/Max predictions:", round(min(bagged.preds), 4), "/", round(max(bagged.preds), 4), "\n")
if (exists("error.metrics")) {
  cat("OOB RMSE:", round(oob.rmse, 4), "\n")
}
cat("Results saved to:", results.pwd, "\n")

#-----------------------------------------------------------

# use the alphavec parameter to explain approaches and results of sensitivity of multi. coefficient adjustment
# see which rebinning techniques work the best.
# use marginal quantiles to split/rebin data 10%, 20% etc. until 2.5kg. 

# explain the different splitting results

# dirichlet uniform: alphavec <- 1*rep(1/ncol(Y.df),ncol(Y.df))
# vs. 
# marginal: alphavec <- 1*colSums(Y.df)/sum(Y.df)

# dont need to explain: quantile idea instead of incremental splits. split into 10 or 20 etc. 

# bootstrap, sampling with replacement. drawing data from multinomial model from probabilities of that cell. 

# dm.tree["splits"]
# dm.tree["variable.importance"]
# dm.tree["cptable"]
# dm.tree["frame"]
# printcp(dm.tree)

# rpart.rules(dm.tree)