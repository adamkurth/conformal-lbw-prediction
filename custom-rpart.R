rm(list=ls())
library(rpart)
library(rpart.plot)
#-----------------------------------------------------------
# A. Load Data & Construct Matrices
#-----------------------------------------------------------
year <- 2021
data.pwd <- sprintf("/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/birthweight_data/rebin/")
load(sprintf("%s/birthweight_data_rebin_%d.RData", data.pwd, year))
Y.counts <- as.matrix(counts.df)
Y.df <- as.data.frame(Y.matrix)
colnames(Y.df) <- paste0("Y",seq(1:ncol(Y.df))) # name columns Y1, Y2, ...
combined.data <- data.frame(X.matrix, Y.df)
response.cols <- colnames(Y.df)
cat("Dimensions of X.matrix:", dim(X.matrix), "\n")
cat("Dimensions of Y.matrix:", dim(Y.matrix), "\n")

# Y = birthweight (actual, log, count)
# X = 7 binary features (from original large.csv of binary values)

#-----------------------------------------------------------
# B. Simplified Dirichlet-Multinomial Functions
#-----------------------------------------------------------
# 'counts' is a numeric vector of counts for each category at a node. 
# 'alpha'  is a numeric vector of Dirichlet parameters, typically alpha_j > 0.
log.dm.likelihood <- function(counts, alpha=1) {
  # counts: vector of counts for each category 
  N <- sum(counts)
  k <- length(counts)
  
  # term1 = ln Γ(Σ alpha_j) - ln Γ(N + Σ alpha_j)
  term1 <- lgamma(k*alpha) - lgamma(N + k*alpha)
  
  # term2 = sum over j of [ ln Γ(count_j + alpha_j) - ln Γ(alpha_j) ]
  term2 <- sum(lgamma(counts + alpha) - lgamma(alpha))
  
  ll <- term1 + term2
  
  if(is.na(ll) || is.infinite(ll)) return(0)  # Handle edge cases
  ll
}

dm.deviance <- function(counts.matrix, alpha=1){
  # Pass counts.matrix directly (already summed per category)
      # cat("input: counts.matrix", counts.matrix, "\n")
  log.dm.likelihood(counts.matrix, alpha=alpha)
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
  # 
  # cat("** myinit() called. Y dimension:", dim(y), "\n")

  list(
    y = y,
    parms = list(alpha =1),
    numresp = ncol(y),
    numy = ncol(y),
    summary = function(yval, dev, wt, ylevel, digits) {
      paste("Deviance:", format(dev, digits = digits))
    },
    text = function(yval, dev, wt, ylevel, digits, n, use.n){
      # yval is the colsum from myeval(); a vector of category counts.
      # keep labels short.
      total.counts <- sum(yval)
      # dev.str <- format(dev, digits=2)
      # lbl <- paste0("Total=", total.counts, ", Dev=", dev.str)
      lbl <- paste0(format(total.counts, digits=4))
      if (use.n) lbl <- paste(lbl, "\nn", n)
      return(lbl)
    }
  )
}

#----------------- 2) eval() ------------------#
myeval <- function(y, wt=NULL, parms=1) {
  # y: response matrix subset for current node
  # wt: weights (not used here)
  
  # notes
  # - we are using a uniform prior (alpha=1) for each category
  # - we are ignoring weights (wt)
  # - we are storing the colSums of y as the "label" for the node
  # - we are using the negative log-likelihood as the deviance
  
  # integrated likelihood of the prior (out of denominator)
  # mysplit = continuous = FALSE
  # NOT USING ANY WEIGHTS in eval/split

  # CALCULATING GOODNESS: 
  # RATIO SCALE: (V(counts_left) + V(counts_left)) / V(counts)
  # LOG SCALE: (V(counts_left) + V(counts_left)) - V(counts)
  
  # Write split function in terms of one X
  counts <- colSums(y)
  dev <- -log.dm.likelihood(counts, alpha=1) # counts is a vector!
  # 'label' can be something to print for the node. We'll just store colSums  
  return(list(label=counts, deviance=dev))
}

#----------------- 3) split() -----------------#
mysplit <- function(y, wt, x, parms, continuous=FALSE) {
  # y: response matrix
  # wt: weights (not used here)
  # x: splitting variable
  # parms: list containing alpha parameter
  # continuous: whether x is continuous
  
  # Using the category counts, we compute the valid deviance improvement for splits.
  
# 
#   This function evaluates all possible splits on 'x'(predictor)
#   and returns 'goodness' + 'direction' for each possible way of splitting 'x'.
# 
#   If x is binary or categorical with levels, rpart calls mysplit with continuous=FALSE
#     and a single pass checking the categories. For each subset left vs right, we measure
#     improvement in deviance.  #
#   The key formula for improvement is:
#     improvement = parent_dev - (left_dev + right_dev).
#   Because rpart tries to *reduce* deviance

  parent.dev <- -log.dm.likelihood(colSums(y), alpha=1)
    # debug
    # cat("input: y",y, "\n")
    # cat("output: parent.dev", parent.dev, "\n")
    
  if(!continuous) {
    # Suppose 'x' is a factor or a binary 0/1 variable
    ux <- sort(unique(x)) # unique values of x
    goodness <- numeric(length(ux) - 1)# store improvement in deviance
    direction <-ux
    
    for(i in 1:(length(ux) - 1)) {
      split.val <- ux[i]
      left.idx <- x == split.val
      left.counts <- colSums(y[left.idx, , drop=FALSE])
      right.counts <- colSums(y[!left.idx, , drop=FALSE])
      child.dev <- -log.dm.likelihood(left.counts) - log.dm.likelihood(right.counts)
      goodness[i] <- parent.dev - child.dev
    }
    

    # 'goodness' is a vector of length = #possible splits. We have only 1 here.
    # 'direction' is a vector of length = #possible splits. Typically -1 means "x < cut" go left
    #   For a categorical, you can store an integer code. We'll just do -1.
    return(list(goodness = pmax(goodness, 0), direction = direction))
  }
}

dm.method <- list(init=myinit, eval=myeval, split=mysplit)

#-----------------------------------------------------------
# D. Fit rpart Tree Using Our Custom DM Method
#-----------------------------------------------------------
X <- X.matrix
Y <- Y.matrix

# my.formula <- Y$actual_weight ~ X$sex + X$dmar + X$mrace15 + X$mager + X$meduc + X$precare5 + X$cig_0
my.formula <- as.formula(
  paste("cbind(", paste(response.cols, collapse = ", "), ") ~ .", sep = "")
)
print(my.formula)

cat("\n--- Fitting the DM-based rpart tree ---\n")

# We set various control parameters to reduce complexity:
#   - usesurrogate=0, maxsurrogate=0 => no surrogate splits
#   - maxcompete=0 => do not evaluate competing splits
#   - xval=0 => no cross-validation
#   - cp=0 => allow splits with minimal improvement
#   - minsplit=5 => each node must have at least 5 obs
#   - maxdepth=5 => limit tree depth to 5

dm.control <- rpart.control(minsplit=0, cp=0, maxdepth=5, xval=0)
dm.tree <- rpart(
  formula = as.formula(paste("cbind(", paste(response.cols, collapse=", "), ") ~ .")),
  data = combined.data,
  method = dm.method,
  control = dm.control
)

#-----------------------------------------------------------
# E. Save and Inspect Results
#-----------------------------------------------------------

save(dm.tree, file=sprintf("%s/dm_tree_rebin_%d.RData", data.pwd, year))
print(dm.tree)
printcp(dm.tree)

# first base R plot
if(nrow(dm.tree$frame) > 1) {
  plot(dm.tree, uniform=TRUE, margin=0.1)
  text(dm.tree, use.n=TRUE, cex=0.5)
} else {
  cat("No splits found - check data or model\n")
}

# rpart.plot
rpart.plot(
  dm.tree,
  # extra = 1,            # or 0, or 101—experiment
  # under = FALSE,        # put node “n=” under the box
  # faclen = 1,           # don’t truncate factor names
  # varlen = 1,           # don’t truncate variable names
  # cex = 1.0,           # shrink text
  # compress = TRUE,      # try to compact the tree horizontally
  # fallen.leaves = FALSE,# place leaves at the bottom
  # tweak = 0.55,         # adjust the spacing between nodes
  # clip.facs = FALSE,     # clip factor levels
  # # shadow.col = "gray",  # color of shadow text
  # nn = FALSE,            # display node numbers
)
#-----------------------------------------------------------


# before rebin: 
# make bar charts /sum(each) throw into pdf
# use x_1, ... x_7 predictors to decode the classes using the legend.

# extract out counts from each bin/node.
# re geneate the counts data for 0.1 kg to create more bins
# write up progress 

# Writing: 
# Page about dataset: US natality data
# Dirichlet Multinomial Criterion 

# More advanced version using statistical uncertainty


