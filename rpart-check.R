#########################################
# rpart-check.R
# ---------------------------------------
# Demonstration of a custom rpart method
# that uses a Dirichlet-Multinomial deviance
# for a *factor* predictor only.
#########################################

rm(list=ls())

# 1) Check we have a normal (recent) rpart
library(rpart)
library(rpart.plot)
cat("rpart package version:", as.character(packageVersion("rpart")), "\n")

##
## 1. Generate small synthetic data
##
set.seed(123)
set.seed(123)
n <- 20
df <- data.frame(
  group = factor(sample(c("A","B","C"), size=n, replace=TRUE))
)
Y <- matrix(rpois(n*3, lambda=5), nrow=n, ncol=3)
colnames(Y) <- c("count_cat1", "count_cat2", "count_cat3")
mydata <- data.frame(df, Y)


##
## 2. Dirichlet-Multinomial log-likelihood
##
log.dm.likelihood <- function(counts, alpha = 1) {
  N <- sum(counts)
  K <- length(counts)
  alpha_0 <- K * alpha
  
  term1 <- lgamma(alpha_0) + lgamma(N+1) - lgamma(N+alpha_0)
  term2 <- sum(lgamma(counts + alpha) - lgamma(alpha) - lgamma(counts + 1))
  
  ll <- term1 + term2
  if (is.na(ll) || is.infinite(ll)) ll <- -Inf
  ll
}

##
## 3. Custom rpart method
##
dm.init <- function(y, offset, parms=NULL, wt=NULL) {
  list(
    method = "dm",
    y = y,
    parms = list(alpha=1),
    numresp = ncol(y),
    numy = ncol(y),
    summary = function(yval, dev, wt, ylevel, digits) {
      paste("Deviance:", format(dev, digits=digits))
    },
    text = function(yval, dev, wt, ylevel, digits, n, use.n){
      total.counts <- sum(yval)
      lbl <- paste0(format(total.counts, digits=4))
      if (use.n) lbl <- paste(lbl, "\nn", n)
      lbl
    }
  )
}

dm.eval <- function(y, wt=NULL, parms=1) {
  counts <- colSums(y)
  dev <- -log.dm.likelihood(counts, alpha=1)
  list(label=counts, deviance=dev)
}

dm.split <- function(y, wt, x, parms, continuous) {
  # Only handle factor splitting
  if (continuous) return(NULL)
  
  # Convert x to a factor if it is numeric (using original levels A/B/C)
  if (is.numeric(x)) {
    x <- factor(
      x, 
      levels = 1:3, 
      labels = c("A", "B", "C")  # Hardcode levels from your synthetic data
    )
  }
  
  # Drop unused levels IN THIS NODE
  x_current <- droplevels(x)
  levs <- levels(x_current)
  nlevs <- length(levs)
  
  # No splits possible if <2 levels
  if (nlevs < 2) {
    return(list(
      goodness = numeric(0),
      direction = integer(0),
      split = matrix(numeric(0), nrow=0, ncol=0)
    ))
  }
  
  parent.dev <- -log.dm.likelihood(colSums(y), alpha=1)
  
  # Initialize candidate splits
  cand.good <- numeric(0)
  cand.dir  <- integer(0)
  cand.mat  <- list()
  
  # Loop over levels PRESENT in this node
  for (i in seq_len(nlevs)) {
    left.idx  <- (x_current == levs[i])
    right.idx <- !left.idx
    
    # Skip splits where a child has zero observations
    if (sum(left.idx) == 0 || sum(right.idx) == 0) next
    
    left.counts  <- colSums(y[left.idx, , drop=FALSE])
    right.counts <- colSums(y[right.idx, , drop=FALSE])
    
    child.dev <- -log.dm.likelihood(left.counts) - log.dm.likelihood(right.counts)
    improvement <- parent.dev - child.dev
    
    if (improvement > 0) {
      cand.good <- c(cand.good, improvement)
      cand.dir  <- c(cand.dir, 2)
      
      # Map split to ORIGINAL factor levels (A/B/C)
      subsetVec <- rep(3, 3)  # 3 = right (original 3 levels)
      subsetVec[which(levels(x) == levs[i])] <- 1  # 1 = left
      cand.mat[[length(cand.mat) + 1]] <- subsetVec
    }
  }
  
  # Return empty structure if no valid splits
  if (length(cand.good) < 1) {
    return(list(
      goodness = numeric(0),
      direction = integer(0),
      split = matrix(numeric(0), nrow=0, ncol=0)
    ))
  }
  
  # Build the split matrix
  split.mat <- do.call(cbind, cand.mat)
  
  list(
    goodness  = cand.good,
    direction = cand.dir,
    split     = split.mat
  )
}
dm.method <- list(init=dm.init, eval=dm.eval, split=dm.split, method="dm")

##
## 4. Fit the tree
##
my.formula <- cbind(count_cat1, count_cat2, count_cat3) ~ group
dm.control <- rpart.control(
  minsplit=2,
  cp=0,
  maxdepth=5,
  xval=0,
  usesurrogate=0,
  maxsurrogate=0,
  maxcompete=0
)

cat("\n--- Fitting the DM-based rpart tree on synthetic data ---\n")
dm.tree <- rpart(
  formula = cbind(count_cat1, count_cat2, count_cat3) ~ group,
  data    = mydata,
  method  = dm.method,
  control = dm.control
)

cat("\n## Resulting tree:\n")
print(dm.tree)
cat("\n## CP Table:\n")
printcp(dm.tree)

if (nrow(dm.tree$frame) > 1) {
  cat("\nWe got at least one split. Plotting...\n")
  rpart.plot(dm.tree)
} else {
  cat("\nNo splits found. Possibly no >0 improvement.\n")
}