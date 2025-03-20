rm(list=ls())
library(rpart)
library(rpart.plot)

#-----------------------------------------------------------
# A. Load Data & Construct Matrices
#-----------------------------------------------------------
year <- 2021

# load WITH 2.5kg 
data.pwd <- sprintf("/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/birthweight_data/rebin/")
load(sprintf("%s/birthweight_data_rebin_%d.RData", data.pwd, year))

# load WITHOUT 2.5kg
# data.pwd <- sprintf("/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/birthweight_data/rebin_without_2.5kg/")
# load(sprintf("%s/birthweight_data_rebin_%d.RData", data.pwd, year))

# load('birthweight_data_2021.Rdata')

Y.df <- counts.df
response.cols <- colnames(Y.df)

alphavec <- 1*rep(1/ncol(Y.df),ncol(Y.df)) # uniform
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
      paste("Deviance:", format(dev, digits = digits))
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
  
  # We only handle the case !continuous, i.e. factor or ordered factor
 
    # levs <- levels(x)         # <-- name them "levs"
    # nlevs <- length(levs)
    # if (nlevs < 2) return(NULL)  # can't split with <2 levels
    # 
    # 
    # i <- 1
    # 
    #   left.idx  <- (x == levs[i])
    #   right.idx <- !left.idx
    #   
    #   if (length(left.idx)>1){
    #   left.counts  <- colSums(y[left.idx,  , drop=FALSE])}else{left.counts <- y[left.idx, , drop = FALSE]}
    #   
    #   
    #   if (length(right.idx)>1){
    #     right.counts  <- colSums(y[right.idx,  , drop=FALSE])}else{right.counts <- y[right.idx, , drop = FALSE]}
    #   
    #   child.dev <- -log.dm.likelihood(left.counts) - log.dm.likelihood(right.counts)
    #   improvement <- parent.dev - child.dev
    #   
    #   goodness  <- max(improvement,0)
    #   direction <- 1  # code "2" is typical for factor splits
    # 
    # 
    # 
    # return(list(
    #   goodness  = goodness,
    #   direction = direction
    # ))
  
  
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

#save(dm.tree, file=sprintf("%s/dm_tree_rebin_%d.RData", data.pwd, year))
print(dm.tree)
printcp(dm.tree)

# # first base R plot
# if(nrow(dm.tree$frame) > 1) {
#   plot(dm.tree, uniform=TRUE, margin=0.1)
#   text(dm.tree, use.n=TRUE, cex=0.5)
# } else {
#   cat("No splits found - check data or model\n")
# }

# rpart.plot
plot(dm.tree)
#-----------------------------------------------------------

# use the alphavec parameter to explain approaches and results of sensitivity of multi. coefficient adjustment
# see which rebinning techniques work the best.
# use marginal quantiles to split/rebin data 10%, 20% etc. until 2.5kg. 

# explain the different splitting results

# dirichlet uniform: alphavec <- 1*rep(1/ncol(Y.df),ncol(Y.df))
# vs. 
# marginal: alphavec <- 1*colSums(Y.df)/sum(Y.df)

# dont need to explain: 
# quantile idea instead of incremental splits. split into 10 or 20 etc. 

# bootstrap, sampling with replacement. drawing data from multinomial model from probabilities of that cell. 