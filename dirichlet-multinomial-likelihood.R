rm(list=ls())
#-----------------------------------------------------------
# A. Load Data & Construct Count Matrix
#-----------------------------------------------------------
library(rpart)
year <- 2021
data.pwd <- "/Users/adamkurth/Documents/RStudio/ms-thesis-kurth/birthweight_data/"
file.path <- sprintf("%s/results%dus.csv", data.pwd, year)
results <- read.csv(file.path)

# Because n_in_node can be very large (on the order of 400k),
# we floor-divide by 100 to keep these values moderate and
# reduce the chance of numeric overflow in the Dirichlet-multinomial.
results$n_in_node <- floor(results$n_in_node / 100)

# count matrix: each row = one "node" or group, each column = a birthweight bin
# Here, "dm.counts" will be an N x K matrix where N=128 (in your case), and K=5 categories.
dm.counts <- as.matrix(results[, c("count_btw_0_0.25kg", 
                                   "count_btw_0.25_0.5kg", 
                                   "count_btw_0.5_0.75kg", 
                                   "count_btw_0.75_1kg", 
                                   "count_above_2.5kg")])

cat("Dimension of dm.counts:", dim(dm.counts), "\n")
cat("First row of dm.counts:", dm.counts[1, ], "\n")

#-----------------------------------------------------------
# B. Define Dirichlet-Multinomial (DM) Log-Likelihood
#-----------------------------------------------------------
log.dirmult.likelihood.row <- function(counts, alpha) {
  # 'counts': integer vector (length K) of observed counts in each category
  # 'alpha': numeric vector (length K) of Dirichlet parameters (priors)
  
  N <- sum(counts) # total number of observations in this row/group
  
  # term1 = ln Γ(Σ alpha_j) - ln Γ(N + Σ alpha_j)
  #   => This is the difference in gamma functions for the sum of alpha vs. the sum plus total N.
  term1 <- lgamma(sum(alpha)) - lgamma(N + sum(alpha))
  
  # term2 = sum over j of [ ln Γ(count_j + alpha_j) - ln Γ(alpha_j) ]
  #   => This is the multinomial-like part, summing individual gamma log-likelihood components.
  term2 <- sum(lgamma(counts + alpha) - lgamma(alpha))
  
  ll <- term1 + term2
  
  # If extreme values cause numeric issues, clamp it to 0 and print a warning
  if (is.na(ll) || is.infinite(ll)) {
    cat("WARNING: DM loglike is NA/Inf. Counts=", counts, " => forcing 0.\n")
    ll <- 0
  }
  ll
}

#-----------------------------------------------------------
# C. Define Custom DM objective function for rpart
#-----------------------------------------------------------
# We define alpha=1 for each category (uniform Dirichlet prior, simple case).
# The objective parts:
#   1) init:  define how 'rpart' initializes the model with the multi-column response.
#   2) eval:  compute the deviance at each node as -2 * (sum of DM log-likelihood).
#   3) split: define how to evaluate potential splits on a predictor.

# How does it split?
# deviance = -2*(sum of DM log-likelihood) sum observations in a node.
# proposed split, partitions those observations in the node to left/right subsets 
# then compute left/right deviances and compare to parent node (called improvement from the split).

K <- ncol(dm.counts)            # number of categories in your multi-column response
alpha.global <- rep(0.05, K)   # uniform Dirichlet prior, alpha=1 for each category

# Lower alpha => DM prior penalizes the "spread-out" categories less, emphasizing the differences among rows. 
#   splits can become more valuable, resulting in deeper trees.
# Higher alpha => DM prior is more "uniform", so within-node variation is not as strongly rewarded or penalize. Often few/no splits occr.
#   splits can become less valuable, resulting in shallower trees.
#----------------- 1) init() ------------------#
myinit <- function(y, offset, parms=NULL, wt=NULL) {
  cat("In myinit, dim(y) =", dim(y), "\n")
  
  if (is.null(parms)) {
    parms <- alpha.global
  }
  
  # We define both `summary` and `text` so rpart knows how to label nodes.
  list(
    y       = y,           
    parms   = parms,       
    numresp = ncol(y),     
    numy    = ncol(y),     
    
    # This is displayed in printcp or summary.rpart for each node
    summary = function(yval, dev, wt, ylevel, digits) {
      paste("Node deviance =", format(dev, digits=digits))
    },
    
    # This function is used by text.rpart to label nodes in the tree plot
    text = function(yval, dev, wt, ylevel, digits, n, use.n) {
      # yval is the node's "label"
      # dev   is the node's deviance
      # n     is the number of observations in each node (or child)
      # use.n is a logical whether to show the #observations
      if (use.n) {
        # Show both deviance and the total n in the node
        paste0("dev=", round(dev, digits), "\nn=", sum(n))
      } else {
        # Just show deviance
        paste("dev=", round(dev, digits))
      }
    }
  )
}

#----------------- 2) eval() ------------------#
myeval <- function(y, wt, parms) {
  # 'y' is the subset of rows in this node,
  # 'wt' is the vector of weights (if used),
  # 'parms' is where alpha is stored.
  
  alpha <- parms
  node.loglike <- 0
  
  # We sum the Dirichlet-multinomial log-likelihood for each row in the node.
  for (i in seq_len(nrow(y))) {
    ll_i <- log.dirmult.likelihood.row(y[i,], alpha)
    node.loglike <- node.loglike + ll_i
  }
  
  # The node deviance is -2 * the total log-likelihood across these rows.
  dev <- -2 * node.loglike
  
  # If dev is infinite or NaN, we clamp it to 0 and warn
  if (is.na(dev) || is.infinite(dev)) {
    cat("WARNING: dev is NA/Inf => forcing 0.\n")
    dev <- 0
  }
  
  # The 'label' for the node is typically the proportion in each category
  total.counts <- colSums(y)
  denom <- sum(total.counts)
  if (denom <= 0) {
    cat("WARNING: sum of counts=0 => using uniform label.\n")
    label <- rep(1/ncol(y), ncol(y))  # fallback uniform if no counts
  } else {
    label <- total.counts / denom
  }
  
  # Return the node label (estimate) and deviance
  list(label=label, deviance=dev, weight=sum(wt))
}

#----------------- 3) split() -----------------#
##
## IMPORTANT: rpart calls this as split(y, wt, x, parms, continuous), not (x,y,wt,...)
##
mysplit <- function(y, wt, x, parms, continuous) {
  #
  # rpart calls this function to evaluate all possible split points on 'x'.
  # 'y', 'wt' are the rows/weights in the current node,
  # 'parms' is alpha, 'continuous' indicates numeric vs factor.
  #
  
  # Attempt to infer predictor name from parent.frame (hacky; optional for debug):
  varname <- NA
  pf <- parent.frame()
  if ("xn" %in% names(pf)) {
    varname <- pf$xn
  }
  
  cat("\n--- mysplit() called on predictor:", varname, "\n")
  cat("continuous =", continuous, "\n")
  cat("range(x) =", range(x), " #unique=", length(unique(x)), "\n")
  
  # If x is a factor, we skip. The code only handles continuous splits:
  if (!continuous) {
    cat("Factor => returning NULL.\n")
    return(NULL)
  }
  
  # Evaluate the deviance of the parent node
  parent.eval <- myeval(y, wt, parms)
  parent.dev  <- parent.eval$deviance
  if (is.na(parent.dev) || is.infinite(parent.dev)) {
    parent.dev <- 0
  }
  
  # Sort x to check potential splits in ascending order
  ord <- order(x)
  x.sorted  <- x[ord]
  y.sorted  <- y[ord, , drop=FALSE]
  wt.sorted <- wt[ord]
  
  # Identify the distinct values of x
  x.unique <- unique(x.sorted)
  n.splits <- length(x.unique) - 1
  
  # If there's only one unique value, no splits are possible:
  if (n.splits < 1) {
    cat("No valid split for this predictor.\n")
    return(NULL)
  }
  
  goodness  <- numeric(n.splits)
  direction <- numeric(n.splits)
  
  # For each possible split, cut at the midpoint between consecutive unique x values.
  for (s in seq_len(n.splits)) {
    cut.val <- (x.unique[s] + x.unique[s+1]) / 2
    idx.left <- (x.sorted < cut.val)
    
    # If split partitions all rows left or all right, no improvement:
    if (!any(idx.left) || all(idx.left)) {
      goodness[s]  <- 0
      direction[s] <- 0
      next
    }
    
    # Evaluate the child nodes' deviances
    left.dev  <- myeval(y.sorted[idx.left, , drop=FALSE],  wt.sorted[idx.left],  parms)$deviance
    right.dev <- myeval(y.sorted[!idx.left, , drop=FALSE], wt.sorted[!idx.left], parms)$deviance
    
    # improvement = parent's dev - (left.dev + right.dev)
    imp <- parent.dev - (left.dev + right.dev)
    if (is.na(imp) || is.infinite(imp)) {
      imp <- 0
    }
    
    cat(sprintf(
      "Split %2d cut=%.5g => parent.dev=%.3f, left.dev=%.3f, right.dev=%.3f, improvement=%.3f\n",
      s, cut.val, parent.dev, left.dev, right.dev, imp
    ))
    
    goodness[s]  <- imp
    # direction < 0 conventionally means "x < cut goes left"
    direction[s] <- -1
  }
  
  # Final clamp of any NA/Inf in goodness
  idx.bad <- which(is.na(goodness) | is.infinite(goodness))
  if (length(idx.bad) > 0) {
    goodness[idx.bad] <- 0
  }
  
  # Return vectors that tell rpart the "quality" of each split
  list(goodness=goodness, direction=direction)
}

# Combine our custom functions into a list for rpart
dm.method <- list(
  init  = myinit, 
  eval  = myeval, 
  split = mysplit
)

#-----------------------------------------------------------
# D. Build the Single Tree using Custom DM Objective
#-----------------------------------------------------------
# 
# We specify a formula with a 5-column response ~ some numeric predictors.
# The "method=dm.method" tells rpart to use our custom deviance code.
#
# We have scaled 'n_in_node' to avoid huge gamma values, but we omit it here
# for demonstration. You can re-add it once you confirm stability.

my.formula <- cbind(
  count_btw_0_0.25kg,
  count_btw_0.25_0.5kg,
  count_btw_0.5_0.75kg,
  count_btw_0.75_1kg,
  count_above_2.5kg
) ~ var_y + sum_y + sum_y_squared + sum_y2_squared + sum_y3 + sum_y3_squared
# If you want to re-add the scaled n_in_node, append + n_in_node

# want to model 4 birthweight buns as DM response, excluding count_above_2.5kg, 
# from the 4 bins, and use count_above_2.5kg as the predictor. 
# i.e.) predict distribution among the 4 lower-birthweight categories based on count_above_2.5kg

cat("\n--- Now fitting the rpart tree with our DM method ---\n")

# We set various control parameters to reduce complexity:
#   - usesurrogate=0, maxsurrogate=0 => no surrogate splits
#   - maxcompete=0 => do not evaluate competing splits
#   - xval=0 => no cross-validation
#   - cp=0 => allow splits with minimal improvement
#   - minsplit=5 => each node must have at least 5 obs
#   - maxdepth=5 => limit tree depth to 5

control = rpart.control(
  minsplit=2, 
  cp=0, 
  maxdepth=10
)

single.tree.dm <- rpart(
  formula  = my.formula,
  data     = results,
  method   = dm.method,
  model    = TRUE,
  na.action= na.omit,
  control  = control
)
#-----------------------------------------------------------
#  E. Inspect / Plot the Resulting Tree
#-----------------------------------------------------------
cat("\n--- Finished building a single DM-based tree ---\n")
print(single.tree.dm)
printcp(single.tree.dm)

# Only plot if there's more than 1 node
if (nrow(single.tree.dm$frame) > 1) {
  plot(single.tree.dm, uniform=TRUE, compress=TRUE, margin=0.1)
  text(single.tree.dm, use.n=TRUE, cex=0.8)
} else {
  cat("No splits were found (just a root node). No tree to plot.\n")
}

#-----------------------------------------------------------
#  F. Predictions
#-----------------------------------------------------------  


#-----------------------------------------------------------
#  Notes
#-----------------------------------------------------------  
# default 'method' options for rpart include:
# - 'anova'  => for continuous response, uses SSE-based deviance
# - 'class'  => for classification, uses Gini/Entropy
# - 'poisson'=> for count data, uses Poisson log-likelihood
#
# With our custom DM method, the deviance = -2 * (sum of DM log-likelihood).
# Splitting tries to reduce that deviance, i.e., find partitions where the DM
# likelihood is improved (higher log-likelihood => lower deviance).