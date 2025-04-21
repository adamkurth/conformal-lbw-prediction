rm(list=ls())
library(rpart)
# library(rpart.plot)

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

#----------------- 3) spit() ------------------#
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
dm.tree <- rpart(
  formula = Y.df~.,
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

plot(dm.tree, main="DM Tree",uniform=TRUE)
text(dm.tree,use.n=TRUE,cex=1,font=3)

#-----------------------------------------------------------
# F. Multinomial Parametric Bootstrap Sampling
#-----------------------------------------------------------
# - Calculate Dirichlet-smoothed multinomial probabilities for each cell/row
# - For each of B iterations:
#   (1) Generate new counts (parametric bootstrap) from these probabilities
#   (2) Sample rows in proportion to row totals for fitting the model (stratified bootstrap)
#   (3) Collect predictions, including OOB

B <- 10000  # Number of bootstrap samples
n.rows <- nrow(X.matrix)
lbw.cols <- 1:10  # Indices of LBW categories (1:10, if 11th is above_2.5kg)
alpha.sig <- 0.05 # Significance level for confidence intervals

# Initialize matrices for predictions
Yhat.oob <- matrix(NA, nrow = n.rows, ncol = B)  # Out-of-bag predictions
Yhat.all <- matrix(0,  nrow = n.rows, ncol = B)  # Full-dataset predictions

# initialize list for tree structure information
tree.structure <- list()
top.vars.used <- character(B)
n.vars.used <- integer(B)

cat("Starting multinomial parametric + category-proportional bootstrap with", B, "samples...\n")
pb <- txtProgressBar(min = 0, max = B, style = 3)

#-----------------------------------------------------------
# 1) Dirichlet-smoothed Multinomial Probabilities per Row
#-----------------------------------------------------------
# For each row (unique predictor combination):
# - Get the observed counts across categories
# - Apply Dirichlet smoothing using the alphavec prior
# - Convert to probabilities for multinomial sampling
multinomial.probs <- matrix(0, nrow = n.rows, ncol = ncol(Y.df))

for(i in 1:n.rows) {
  
    # Get counts for this cell/row
    cell.counts <- as.numeric(Y.df[i,])
    
    # Calculate probability with Dirichlet smoothing
    # (counts + alpha) / (sum(counts) + sum(alpha))
    cell.probs <- (cell.counts + alphavec) / (sum(cell.counts) + sum(alphavec))
    
    multinomial.probs[i, ] <- cell.probs
}

#-----------------------------------------------------------
# 2) Prepare Weights for Row-Bootstrap
#-----------------------------------------------------------
# Each row in counts.df represents a unique predictor combination
# We want to sample rows in proportion to their occurrence frequency
row.totals <- rowSums(counts.df)
row.probs <- row.totals / sum(row.totals)  # Probability of selecting each row

# Sample size for row-level bootstrap
sample.size <- n.rows


#-----------------------------------------------------------
# 3) Main Loop over B Bootstrap Samples
#-----------------------------------------------------------
for (b in seq_len(B)) {
  
    #---------------------------------------
    # (A) Parametric bootstrap of responses
    #---------------------------------------
    # Generate new counts from the multinomial distribution
    # for each row, using the cell's original total count
    bootstrap.counts <- matrix(0, nrow = n.rows, ncol = ncol(Y.df))
      
      for (i in seq_len(n.rows)) {
          # get original total count for this predictor combination
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

    #---------------------------------------
    # (B) Row bootstrap in proportion to row totals
    #---------------------------------------
    # Generate a vector of row indices with replacement
    # where rows are sampled with probability proportional to row.probs
    row.indicies <- sample(seq_len(n.rows),
                           size = sample.size,
                           replace = TRUE,
                           prob = row.probs)
    
    # Out-of-bag indices (rows not selected in this bootstrap sample)
    oob.idx <- setdiff(seq_len(n.rows), unique(row.indicies))
  
  #---------------------------------------
  # (C) Fit rpart to row- and param-boot data
  #---------------------------------------
  # The "dm.method" expects a matrix of counts (the response)
  # combined with predictor columns. thus:   bootstrap.counts.df[idx, ] ~ .
  # and X.matrix[idx, ] for the predictor data.
  
  model.data <- cbind(
    bootstrap.counts.df[row.indicies, ], # response
    X.matrix[row.indicies, ]             # predictors
  )
  
  # create formula for rpart
  response.cols <- colnames(bootstrap.counts.df) # extract count colnames
  response.cbind <- paste(response.cols,collapse = ",") # counts_..., counts_..., 
  formula <- paste0("cbind(", response.cbind, ") ~ .") # cbind(counts_..., counts_...) ~ .
  r.formula <- as.formula(formula) # convert to formula object
  
  dm.tree.b <- rpart(
    formula = r.formula,
    data = model.data,
    method  = dm.method,
    control = dm.control
  )
  
  
  #---------------------------------------
  # (D) Extract tree structure information
  #---------------------------------------
  # Get variables used in this tree
  vars.used <- unique(dm.tree.b$frame$var[dm.tree.b$frame$var != "<leaf>"])
  n.vars.used[b] <- length(vars.used)
  
  # store top variables (first split)
  if (length(vars.used) > 0 ){
      
      top.vars.used[b] <- as.character(dm.tree.b$frame$var[1])
  } else {
      top.vars.used[b] <- "none" # no splits just root
      
  }
  
  # store tree structure for later 
  tree.structure[[b]] <- list(
    variables = vars.used,
    frame = dm.tree.b$frame,
    splits = dm.tree.b$splits
  )
  
  #---------------------------------------
  # (E) Fit reduced model (smaller tree)
  #---------------------------------------
  # Create a more restricted tree (smaller)
  dm.control.small <- rpart.control(minsplit = 5, cp = 0.01, maxdepth = 3, xval = 0, usesurrogate = 0)
  
  dm.tree.small.b <- rpart(
    formula = r.formula,
    data = model.data,
    method = dm.method,
    control = dm.control.small
  )
  
  # compare smaller tree structure to full tree
  vars.uses.small <- unique(dm.tree.small.b$frame$var[dm.tree.small.b$frame$var != "<leaf>"])
  
  #---------------------------------------
  # (F) Predictions: OOB and full dataset
  #---------------------------------------
  # Out-of-bag predictions (if any OOB observations)
  
  if(length(oob.idx > 0)) {
      
      preds.oob <- predict(
          dm.tree.b, 
          newdata=data.frame(X.matrix[oob.idx, , drop=FALSE]), 
          type="matrix"
      )
      
      # calculate LBW prob (sum of first 10 categories)
      # normalize by rowsums to get valid prob
      lbw.prob.oob <- rowSums(preds.oob[, lbw.cols]) / rowSums(preds.oob)
      Yhat.oob[oob.idx, b] <- lbw.prob.oob
      
  }
  
  
  # Full dataset predictions
  preds.all <- predict(
    dm.tree.b, 
    newdata = data.frame(X.matrix), 
    type = "matrix"
  )
  
  # Calculate LBW probability for all rows
  lbw.prob.all <- rowSums(preds.all[, lbw.cols]) / rowSums(preds.all)
  Yhat.all[, b] <- lbw.prob.all
  
  setTxtProgressBar(pb, b) # Update the progress bar
  
}
close(pb)
cat("\n--- Multinomial + category-proportional bootstrap sampling completed ---\n")


#-----------------------------------------------------------
# G. Decision Tree Depth Analysis for cig_0 Predictor
#-----------------------------------------------------------
cat("\n--- Analyzing cig_0 predictor across different tree depths ---\n")
depths <- c(2, 3, 4, 5)
trees <- list()
tree.sum <- list()
var.mats <- list()

cat("\n--- Comparing Decision Trees at Different Depths ---\n")

for (depth in depths) {
    
    cat(sprintf("\nFitting tree with max depth = %d...\n", depth))
    
    # control parameters with specific depths
    depth.control <- rpart.control(
      minsplit = 2,     # min obs in node for split
      cp = 0,           # complexity param
      maxdepth = depth, # max depth 
      xval = 0,         # no cross-validation
      usesurrogate = 0  # no surrogate splits
    )
    
    depth.tree <- rpart(
      formula = r.formula,
      data = model.data,
      method = dm.method,
      control = depth.control
    )
  
    # store tree 
    trees[[as.character(depth)]] <- depth.tree
    
    # extract info about variables used
    vars.used <- unique(depth.tree$frame$var[depth.tree$frame$var != "<leaf>"])
    
    tree.sum[[as.character(depth)]] <- list(
      variables = vars.used,
      n.variables = length(vars.used),
      n.terminal.nodes = sum(depth.tree$frame$var == "<leaf>"),
      first.split = ifelse(length(vars.used) > 0, 
                          as.character(depth.tree$frame$var[1]), 
                          "none")
    )
    
    # check if "cig_0" is used in this tree
    cig.used <- "cig_0" %in% vars.used
    
    if (cig.used){
        
        # find where cig0 used
        cig.nodes <- which(depth.tree$frame$var == "cig_0")
        node.depths <- floor(log2(cig.nodes))
        cat("cig_0 appears at tree level(s):", paste(unique(node.depths), collapse = ", "), "\n")
        
    } else {

        cat("cig_0 not used in this tree.\n")
    }
    
      # create presence/absence matrix for variables used
      # visualize which variables used at different depths
      var.names <- colnames(X.matrix)
      var.mat <- matrix(0, nrow = length(var.names), ncol = depth)
      rownames(var.mat) <- var.names  
      colnames(var.mat) <- paste0("level_", 1:depth)
      
      # for each node in tree 
      for (i in 1:nrow(depth.tree$frame)) {
        
          var.name <- depth.tree$frame$var[i]
          
          if (var.name != "<leaf>"){
              # calculate the level of this node (approximately)
              node.level <- min(floor(log2(i)) + 1, depth)
              var.mat[var.name, node.level] <- 1
          }
      }
      
      var.mats[[as.character(depth)]] <- var.mat
      
      # Print a summary of tree structure
      cat("\nTree Summary (Depth =", depth, "):\n")
      cat("Number of terminal nodes:", sum(depth.tree$frame$var == "<leaf>"), "\n")
      cat("Number of variables used:", length(vars.used), "\n")
      cat("First split on variable:", ifelse(length(vars.used) > 0, as.character(depth.tree$frame$var[1]), "none"), "\n")
  
       # print top 3 splits
      if (length(vars.used) >= 3) {
          top.splits <- as.character(depth.tree$frame$var[1:min(3, length(vars.used) + 1)])
          top.splits <- top.splits[top.splits != "<leaf>"]
          cat("Top 3 splits: ", paste(top.splits, collapse = ", "), "\n")
      }
      
      # Predict LBW probabilities with this tree
      preds <- predict(depth.tree, newdata = data.frame(X.matrix), type = "matrix")
      lbw.prob <- rowSums(preds[, lbw.cols]) / rowSums(preds)
      
      # Store prediction statistics directly in tree summary (not in a nested pred.stats list)
      tree.sum[[as.character(depth)]]$mean.lbw.prob = mean(lbw.prob)
      tree.sum[[as.character(depth)]]$median.lbw.prob = median(lbw.prob)
      tree.sum[[as.character(depth)]]$min.lbw.prob = min(lbw.prob)
      tree.sum[[as.character(depth)]]$max.lbw.prob = max(lbw.prob)
      tree.sum[[as.character(depth)]]$sd.lbw.prob = sd(lbw.prob)
      # These are already stored earlier, but we'll store them again for consistency
      tree.sum[[as.character(depth)]]$n.terminal.nodes = sum(depth.tree$frame$var == "<leaf>")
      tree.sum[[as.character(depth)]]$n.variables = length(vars.used)
      
      cat("LBW Probability - Mean:", round(mean(lbw.prob), 4), 
          "Min:", round(min(lbw.prob), 4), 
          "Max:", round(max(lbw.prob), 4), "\n")
}


#-----------------------------------------------------------
# Comparative Analysis of Trees
#-----------------------------------------------------------
cat("\n--- Comparative Analysis of Trees at Different Depths ---\n")

# Compare variable usage across depths
var.across.depths <- data.frame(
  variable = colnames(X.matrix)
)

for(depth in depths) {
  
  depth.str <- as.character(depth)
  
  var.across.depths[[paste0("depth_", depth)]] <- 
    ifelse(var.across.depths$variable %in% tree.sum[[depth.str]]$variables, 1, 0)
  
}

# Calculate total usage across all depths
var.across.depths$total.usage <- rowSums(var.across.depths[, paste0("depth_", depths)])

# Sort by total usage
var.across.depths <- var.across.depths[order(-var.across.depths$total.usage), ]

cat("\nVariable Usage Across Different Tree Depths:\n")
print(var.across.depths)

cat("\nFocus on cig_0 Predictor:\n")
cig.usage <- var.across.depths[var.across.depths$variable == "cig_0", ]
print(cig.usage)

pred.stats.compare <- data.frame(
  depth = depths,
  mean.lbw.prob = sapply(as.character(depths), function(d) tree.sum[[d]]$mean.lbw.prob),
  min.lbw.prob = sapply(as.character(depths), function(d) tree.sum[[d]]$min.lbw.prob),
  max.lbw.prob = sapply(as.character(depths), function(d) tree.sum[[d]]$max.lbw.prob),
  sd.lbw.prob = sapply(as.character(depths), function(d) tree.sum[[d]]$sd.lbw.prob),
  n.terminal.nodes = sapply(as.character(depths), function(d) tree.sum[[d]]$n.terminal.nodes),
  n.variables = sapply(as.character(depths), function(d) tree.sum[[d]]$n.variables)
)


cat("\nPrediction and Tree Structure Comparison:\n")
print(pred.stats.compare)

# check consistency
first.splits <- sapply(as.character(depths), function(d) tree.sum[[d]]$first.split)
cat("\nFirst Split Variable at Each Depth:", paste(first.splits, collapse = ", "), "\n")

if (length(unique(first.splits)) == 1){
  
    cat("First split is consistent across all depths\n")
} else {
    
    cat("First split varies across different depths\n")
}

save(trees, tree.sum, var.mats, var.across.depths, pred.stats.compare,
     file = sprintf("%s/tree_depth_comparison_%d.RData", results.pwd, year))


plot.data <- data.frame(
  depth = rep(depths, 3),
  type = c(rep("min", length(depths)), rep("mean", length(depths)), rep("max", length(depths))),
  value = c(
    pred.stats.compare$min.lbw.prob,
    pred.stats.compare$mean.lbw.prob,
    pred.stats.compare$max.lbw.prob
  )
)
p.stats <- ggplot(plot.data, aes(x = factor(depth), y = value, group = type, color = type)) +
  geom_line() +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "LBW Probability Statistics by Tree Depth",
       x = "Maximum Tree Depth",
       y = "LBW Probability",
       color = "Statistic") +
  scale_color_brewer(palette = "Set1")
ggsave(sprintf("%s/tree_depth_comparison_%d.png", boot.pwd, year), p.stats, width = 8, height = 6)


p.complexity <- ggplot(pred.stats.compare, aes(x = factor(depth))) +
  geom_line(aes(y = n.terminal.nodes, group = 1), color = "blue") +
  geom_point(aes(y = n.terminal.nodes), color = "blue", size = 3) +
  geom_line(aes(y = n.variables * 5, group = 1), color = "red") +  # Scale up for visibility
  geom_point(aes(y = n.variables * 5), color = "red", size = 3) +
  scale_y_continuous(
    name = "Number of Terminal Nodes",
    sec.axis = sec_axis(~ . / 5, name = "Number of Variables Used")
  ) +
  theme_minimal() +
  labs(title = "Tree Complexity by Maximum Depth",
       x = "Maximum Tree Depth") +
  theme(
    axis.title.y.left = element_text(color = "blue"),
    axis.title.y.right = element_text(color = "red")
  )

# Save the complexity plot
ggsave(sprintf("%s/tree_complexity_comparison_%d.png", boot.pwd, year), p.complexity, width = 8, height = 6)


#-----------------------------------------------------------
# H) Analyze Bootstrap Tree Structure Results
#-----------------------------------------------------------
# Summarize variable usage across bootstrap samples
var.names <- colnames(X.matrix) 
var.usage <- matrix(0, nrow = B, ncol=length(var.names))
colnames(var.usage) <- var.names

for ( b in 1:B ){
  
  if(length(tree.structure[[b]]$variables) > 0 ){
    
    for (var in tree.structure[[b]]$variables){
      
      if(var %in% var.names) {
        var.usage[b,var] <- 1
      }
    }
  }
}

# calclate frequency table of variables used
var.freq <- colSums(var.usage) / B
var.freq.df <- data.frame(
  variable = var.names,
  frequency = var.freq
)
var.freq.df <- var.freq.df[order(-var.freq.df$frequency),]


# top split variable 
top.var.freq <- table(top.vars.used) / B 
top.var.df <- data.frame(
  variable = names(top.var.freq),
  frequency = as.numeric(top.var.freq)
)
top.var.freq <- top.var.df[order(-top.var.df$frequency),]


# analyze number of variables used 
n.vars.summary <- data.frame(
  mean = mean(n.vars.used),
  median = median(n.vars.used),
  min = min(n.vars.used),
  max = max(n.vars.used), 
  sd = sd(n.vars.used)
)

#-----------------------------------------------------------
cat("\n--- Bootstrap Tree Structure Analysis ---\n")
cat("\nFrequency of Variables Used in Trees:\n")
print(var.freq.df)

cat("\nFrequency of Top Split Variables:\n")
print(top.var.freq)

cat("\nSummary of Number of Variables Used:\n")
print(n.vars.summary)



#-----------------------------------------------------------
# I) Continue with prediction interval calculations as before
#-----------------------------------------------------------
# Calculate bagged predictions from all bootstrap samples
bagged.preds <- rowMeans(Yhat.all)

# Calculate OOB predictions, ignoring NAs
oob.preds <- rep(NA, n.rows)

for(i in 1:n.rows) {
  
  # Get valid predictions for this row
  valid_preds <- Yhat.oob[i, !is.na(Yhat.oob[i, ])]
  
  if(length(valid_preds) > 0) {
    
    # If we have valid predictions, use their mean
    oob.preds[i] <- mean(valid_preds)
  } else {
    
    # If all predictions are NA, use the bagged prediction
    oob.preds[i] <- bagged.preds[i]
    cat("Node", i, "had no OOB samples, using bagged prediction instead\n")
  }
}

# Calculate confidence intervals
lower.ci <- apply(Yhat.all, 1, function(x) quantile(x, alpha.sig / 2))
upper.ci <- apply(Yhat.all, 1, function(x) quantile(x, 1 - alpha.sig / 2))
interval.width <- upper.ci - lower.ci

# Calculate prediction standard errors
pred.se <- apply(Yhat.all, 1, function(x) sd(x))

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

bootstrap.tree.results <- list(
  pred.results = pred.results,
  # tree.structure = tree.structure,
  var.freq.df = var.freq.df,
  top.var.df = top.var.df,
  n.vars.summary = n.vars.summary
)

library(dplyr)
library(ggplot2)
library(scales)
library(gridExtra)
library(viridis)

# Original check.agreement function (unchanged)
check.agreement <- function(data, col1, col2, decimal=4) {

      # round both cols to specific decimal
    round.col.1 <- round(data[[col1]], decimal)
    round.col.2 <- round(data[[col2]], decimal)
    
    # agree? 
    agreement <- round.col.1 == round.col.2
    
    # calc % agreement
    percent.agreement <- mean(agreement, na.rm = TRUE) * 100
    
    # data.frame of disagreements 
    disagreements <- data[!agreement, c("node", col1, col2)] %>%
      mutate(
        diff = abs(!!sym(col1) - !!sym(col2)),
        diff.rounded = abs(round.col.1[!agreement] - round.col.2[!agreement])
      )
    
    return(list(
      agreement.percent = percent.agreement,
      num.agreements = sum(agreement, na.rm = TRUE),
      num.disagreements = sum(!agreement, na.rm = TRUE),
      total.nodes = length(agreement), 
      disagreements = disagreements,
      mean.diff = mean(abs(data[[col1]] - data[[col2]]), na.rm = TRUE),
      max.diff = max(abs(data[[col1]] - data[[col2]]), na.rm = TRUE)
    ))
}






preds.sorted <- pred.results[order(pred.results$bagged.pred), ]
lowest.5.pred <- head(preds.sorted, 5)
highest.5.pred <- tail(preds.sorted, 5)
highest.5.pred <- highest.5.pred[order(-highest.5.pred$bagged.pred), ]

find.common.preds <- function(data) {
  
  preds <- c("sex", "dmar", "mrace15", "mager", "meduc", "precare5", "cig_0")
  common <- c()
  
  for (pred in preds) {
      values <- unique(data[[pred]])
      if (length(values) == 1) {
        common <- c(common, pred)
      }
    }
    
  return(common)
}

lowest.5.common <- find.common.preds(lowest.5.pred)
highest.5.common <- find.common.preds(highest.5.pred)

lowest.vals <- c(); highest.vals <- c()
for (pred in lowest.5.common) {
   lowest.vals <- c(lowest.vals, paste0(pred, " = ", lowest.5.pred[1, pred]))
}
for (pred in highest.5.common) {
   highest.vals <- c(highest.vals, paste0(pred, " = ", highest.5.pred[1, pred]))
}

combined.data <- rbind(
  cbind(Group = "Lowest 5", lowest.5.pred[, c("node", "bagged.pred", "sex", "dmar", "mrace15", "mager", "meduc", "precare5", "cig_0")]),
  cbind(Group = "Highest 5", highest.5.pred[, c("node", "bagged.pred", "sex", "dmar", "mrace15", "mager", "meduc", "precare5", "cig_0")])
)

pred.desc <- function(data){
    apply(data[, c("sex", "dmar", "mrace15", "mager", "meduc", "precare5", "cig_0")], 1, function(row) {
      sex_desc <- ifelse(row["sex"] == 1, "Male", "Female")
      dmar_desc <- ifelse(row["dmar"] == 1, "Married", "Not Married")
      mrace15_desc <- ifelse(row["mrace15"] == 1, "Black", "Not Black")
      mager_desc <- ifelse(row["mager"] == 1, "Age > 33", "Age ≤ 33")
      meduc_desc <- ifelse(row["meduc"] == 1, "High School", "Not High School")
      precare5_desc <- ifelse(row["precare5"] == 1, "Full Prenatal", "Not Full Prenatal")
      cig_0_desc <- ifelse(row["cig_0"] == 1, "Smoker", "Non-smoker")
      
      paste(sex_desc, dmar_desc, mrace15_desc, mager_desc, meduc_desc, precare5_desc, cig_0_desc, sep = ", ")
    })
}

combined.data$Description <- pred.desc(data=combined.data)


create_latex_table <- function(data, common_low, common_high) {

  require(xtable)
  require(knitr)
  
  # Create detailed table
  detailed_table <- data[, c("Group", "node", "bagged.pred", "Description")]
  colnames(detailed_table) <- c("Group", "Node", "Bagged Probability", "Predictor Combination")
  
  # Create summary table of common predictors
  common_table <- data.frame(
    Group = c("Lowest 5", "Highest 5"),
    Common_Predictors = c(paste(common_low, collapse = ", "), paste(common_high, collapse = ", ")),
    Common_Values = c(paste(lowest.vals, collapse = ", "), paste(highest.vals, collapse = ", ")),
    Mean_Bagged_Prob = c(mean(lowest.5.pred$bagged.pred), mean(highest.5.pred$bagged.pred))
  )
  colnames(common_table) <- c("Group", "Common Predictors", "Common Values", "Mean Bagged Probability")
  
  # Convert to LaTeX format
  detailed_xtable <- xtable(detailed_table, 
                            caption = "Detailed Comparison of Lowest and Highest Bagged Probability Classes",
                            label = "tab:detailed_comparison")
  
  common_xtable <- xtable(common_table,
                          caption = "Common Predictors Among Lowest and Highest Bagged Probability Classes",
                          label = "tab:common_predictors")
  
  # Output LaTeX tables
  detailed_latex <- print(detailed_xtable, 
                          include.rownames = FALSE,
                          floating = TRUE,
                          table.placement = "htbp",
                          sanitize.text.function = function(x) {x},
                          caption.placement = "top",
                          print.results = FALSE)
  
  common_latex <- print(common_xtable,
                        include.rownames = FALSE,
                        floating = TRUE, 
                        table.placement = "htbp",
                        sanitize.text.function = function(x) {x},
                        caption.placement = "top",
                        print.results = FALSE)
  
  # Create full LaTeX document
  latex_doc <- c(
    "\\documentclass{article}",
    "\\usepackage{booktabs}",
    "\\usepackage{siunitx}",
    "\\usepackage{caption}",
    "\\usepackage{threeparttable}",
    "\\usepackage{array}",
    "\\newcolumntype{L}{>{\\raggedright\\arraybackslash}X}",
    "\\begin{document}",
    "",
    "% Common predictors table",
    common_latex,
    "",
    "% Detailed comparison table",
    detailed_latex,
    "",
    "\\end{document}"
  )
  
  # Save LaTeX document
  writeLines(latex_doc, "bagged_probability_tables.tex")
  
  # Return both tables
  return(list(detailed_table = detailed_table, common_table = common_table))
}


tables <- create_latex_table(combined.data, lowest.vals, highest.vals)



ggplot(plot.data, aes(x = reorder(as.factor(node), bagged.pred), y = bagged.pred, fill = group)) +
  geom_bar(stat = "identity") +
  labs(title = "Lowest and Highest Bagged Probabilities by Node",
       x = "Node",
       y = "Bagged Probability") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

sum.table <- data.frame(
  Group = c("Lowest 5", "Highest 5"),
  Common_Predictors = c(paste(lowest.vals, collapse = ", "), 
                        paste(highest.vals, collapse = ", ")),
  Mean_Bagged_Prob = c(mean(lowest.5.pred$bagged.pred), mean(highest.5.pred$bagged.pred))
)




#-------------------------------------------------------------
# J. Visualize Bootstrap Results
#-------------------------------------------------------------
library(ggplot2)
library(tidyr)

create_publication_visuals <- function(data, agreement_results, decimal_places) {
  
    # Create a descriptive title based on decimal places
    precision_description <- switch(as.character(decimal_places),
                                    "1" = "one decimal place",
                                    "2" = "two decimal places",
                                    "3" = "three decimal places",
                                    "4" = "four decimal places",
                                    "5" = "five decimal places",
                                    paste(decimal_places, "decimal places"))
    
    # Calculate absolute differences for all points
    comparison_data <- data %>%
      mutate(
        diff = abs(bagged.pred - oob.pred),
        agrees = round(bagged.pred, decimal_places) == round(oob.pred, decimal_places),
        # Create a categorical variable for magnitude of difference
        diff_category = case_when(
          diff < 10^-(decimal_places+1) ~ "Virtually identical",
          diff < 10^-decimal_places ~ "Minor difference",
          diff < 10^-(decimal_places-1) ~ "Notable difference",
          TRUE ~ "Substantial difference"
        ),
        # Convert to factor with ordered levels
        diff_category = factor(diff_category, 
                               levels = c("Virtually identical", "Minor difference", 
                                          "Notable difference", "Substantial difference"))
      )
    
    # 1. Create a scatterplot with enhanced visual elements
    scatter_plot <- ggplot(comparison_data, aes(x = bagged.pred, y = oob.pred)) +
      # Add reference line for perfect agreement
      geom_abline(intercept = 0, slope = 1, color = "darkgray", linetype = "dashed", size = 0.8) +
      
      # Add a subtle confidence band (optional)
      geom_smooth(method = "lm", color = "#3366CC", fill = "#3366CC", alpha = 0.1, se = TRUE) +
      
      # Add points with color based on agreement category
      geom_point(aes(color = diff_category), size = 3, alpha = 0.7) +
      
      # Use a colorblind-friendly palette
      scale_color_viridis_d(option = "plasma", end = 0.9, direction = -1, name = "Difference Category") +
      
      # Format axes as percentages for better readability
      scale_x_continuous(labels = percent_format(accuracy = 0.1)) +
      scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
      
      # Add informative labels
      labs(
        title = "Calibration Plot: Bagged vs. Out-of-Bag Predictions",
        subtitle = paste0(round(agreement_results$agreement.percent, 1), 
                          "% agreement to ", precision_description, 
                          " (", agreement_results$num.agreements, 
                          " of ", agreement_results$total.nodes, " nodes)"),
        x = "Bagged Predictions (Probability)",
        y = "Out-of-Bag Predictions (Probability)",
        caption = paste("Mean absolute difference:", 
                        format(agreement_results$mean.diff, scientific = FALSE, digits = 5),
                        "| Maximum difference:", 
                        format(agreement_results$max.diff, scientific = FALSE, digits = 5))
      ) +
      
      # Use publication-quality theme
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9),
        axis.title = element_text(size = 11, face = "bold"),
        panel.grid.minor = element_line(color = "gray95"),
        panel.grid.major = element_line(color = "gray90"),
        panel.border = element_rect(fill = NA, color = "gray80", size = 0.5)
      ) +
      
      # Force 1:1 aspect ratio
      coord_equal(xlim = c(min(comparison_data$bagged.pred) * 0.95, 
                           max(comparison_data$bagged.pred) * 1.05),
                  ylim = c(min(comparison_data$oob.pred) * 0.95, 
                           max(comparison_data$oob.pred) * 1.05))
    
    # 2. Create a Bland-Altman plot (agreement plot) for detailed analysis
    agreement_plot <- ggplot(comparison_data, 
                             aes(x = (bagged.pred + oob.pred)/2, y = bagged.pred - oob.pred)) +
      # Add reference line at zero
      geom_hline(yintercept = 0, color = "darkgray", linetype = "dashed", size = 0.8) +
      
      # Add points with color based on agreement
      geom_point(aes(color = diff_category), size = 3, alpha = 0.7) +
      
      # Use the same colorblind-friendly palette
      scale_color_viridis_d(option = "plasma", end = 0.9, direction = -1, name = "Difference Category") +
      
      # Format x-axis as percentages
      scale_x_continuous(labels = percent_format(accuracy = 0.1)) +
      
      # Add mean and +/- 1.96 SD reference lines
      geom_hline(yintercept = mean(comparison_data$bagged.pred - comparison_data$oob.pred), 
                 color = "#3366CC", linetype = "solid", size = 1) +
      geom_hline(yintercept = mean(comparison_data$bagged.pred - comparison_data$oob.pred) + 
                   1.96 * sd(comparison_data$bagged.pred - comparison_data$oob.pred), 
                 color = "#3366CC", linetype = "dotted", size = 0.8) +
      geom_hline(yintercept = mean(comparison_data$bagged.pred - comparison_data$oob.pred) - 
                   1.96 * sd(comparison_data$bagged.pred - comparison_data$oob.pred), 
                 color = "#3366CC", linetype = "dotted", size = 0.8) +
      
      # Add informative labels
      labs(
        title = "Bland-Altman Plot: Agreement Between Prediction Methods",
        subtitle = "Differences between bagged and out-of-bag predictions",
        x = "Mean of Bagged and OOB Predictions",
        y = "Difference (Bagged - OOB)",
        caption = "Solid blue line: mean difference | Dotted lines: 95% limits of agreement"
      ) +
      
      # Use publication-quality theme
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9),
        axis.title = element_text(size = 11, face = "bold"),
        panel.grid.minor = element_line(color = "gray95"),
        panel.grid.major = element_line(color = "gray90"),
        panel.border = element_rect(fill = NA, color = "gray80", size = 0.5)
      )
    
    # Return both plots
    return(list(
      calibration_plot = scatter_plot,
      agreement_plot = agreement_plot
    ))
}

# Run the original agreement check
agree.results <- check.agreement(data=bootstrap.tree.results$pred.results,
                                 col1="bagged.pred",
                                 col2="oob.pred",
                                 decimal=2)

# Create publication-quality visualizations
publication_plots <- create_publication_visuals(
  bootstrap.tree.results$pred.results, 
  agree.results, 
  decimal_places = 1
)

# Save the plots in high quality
ggsave("bagged_oob_calibration_plot.pdf", publication_plots$calibration_plot, 
       width = 8, height = 7, device = cairo_pdf)
ggsave("bagged_oob_agreement_plot.pdf", publication_plots$agreement_plot, 
       width = 8, height = 7, device = cairo_pdf)

# Also save PNG versions for quick viewing
ggsave("bagged_oob_calibration_plot.png", publication_plots$calibration_plot, 
       width = 8, height = 7, dpi = 300)
ggsave("bagged_oob_agreement_plot.png", publication_plots$agreement_plot, 
       width = 8, height = 7, dpi = 300)

# Print the agreement statistics
cat("Agreement to", agree.results$decimal, "decimal places:", 
    round(agree.results$agreement.percent, 2), "%\n")
cat("Number of agreements:", agree.results$num.agreements, 
    "out of", agree.results$total.nodes, "nodes\n")
cat("Mean absolute difference:", agree.results$mean.diff, "\n")
cat("Maximum absolute difference:", agree.results$max.diff, "\n")

# Display the plots
grid.arrange(publication_plots$calibration_plot, publication_plots$agreement_plot, ncol = 1)
ggsave("combined_publication_plots.png", 
       path = sprintf("%s/combined_publication_plots_%d.png", boot.pwd, year),
       plot = grid.arrange(publication_plots$calibration_plot, publication_plots$agreement_plot, ncol = 1), 
       width = 8, height = 14, dpi = 300)



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
  geom_line(color = "red") +
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
save(bootstrap.tree.results,file = sprintf("%s/bootstrap_tree_results_%d.RData", results.pwd, year))

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