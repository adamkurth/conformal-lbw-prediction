rm(list=ls())
library(rpart)
library(ggplot2)
library(xtable)

#-----------------------------------------------------------
# A. Load Data & Construct Matrices
#-----------------------------------------------------------
year <- 2021
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

print(alphavec)
length(alphavec)
Y.df <- counts.df; response.cols <- colnames(Y.df)

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
  return(list(goodness = pmax(goodness, 0), direction = direction))

}
 
#----------------- 4) pred() ------------------#
mypred <- function(fit, newdata = NULL, type = c("vector", "prob", "class", "matrix")) {
  type <- match.arg(type)
  
  # Create a modified copy of the tree to avoid changing the original
  fit_copy <- fit
  
  # Check and repair the var field for the root node if needed
  root_row <- which(rownames(fit_copy$frame) == "1")
  if (length(root_row) > 0) {
    if (is.na(fit_copy$frame$var[root_row]) || fit_copy$frame$var[root_row] == "<leaf>") {
      cat("WARNING: Root node var is NA or leaf - attempting to fix\n")
      
      # Try to infer the root variable from splits
      split_vars <- rownames(fit_copy$splits)
      if (length(split_vars) > 0) {
        # Find the first variable used in splits
        root_var <- split_vars[1]
        fit_copy$frame$var[root_row] <- root_var
        cat("Setting root node var to:", root_var, "\n")
      } else {
        cat("ERROR: Could not determine root variable from splits\n")
      }
    }
  }
  
  
  # Determine terminal nodes for each observation
  if (is.null(newdata)) {
    where <- fit_copy$where
  } else {
    newdata <- as.data.frame(newdata)
    where <- integer(nrow(newdata))
    
    for (i in 1:nrow(newdata)) {
      # Start at root
      node <- 1
      
      while (TRUE) {
        node_str <- as.character(node)
        if (!node_str %in% rownames(fit_copy$frame)) break
        
        node_row <- which(rownames(fit_copy$frame) == node_str)
        varname <- fit_copy$frame$var[node_row]
        
        if (is.na(varname) || varname == "<leaf>") break
        
        # Find split information for this variable
        split_rows <- which(rownames(fit_copy$splits) == varname)
        if (length(split_rows) == 0) break
        
        # Get the split value
        split_val <- fit_copy$splits[split_rows[1], "index"]
        
        # Get current observation's value for the split variable
        if (!varname %in% colnames(newdata)) {
          cat("WARNING: Variable", varname, "not found in newdata for observation", i, "\n")
          break
        }
        
        cur_val <- newdata[i, varname]
        
        # Determine direction based on split
        go_left <- !is.na(cur_val) && as.numeric(cur_val) <= split_val
        
        # Move to child node or stop if not available
        child_node <- if (go_left) 2*node else 2*node + 1
        child_str <- as.character(child_node)
        
        if (child_str %in% rownames(fit_copy$frame)) {
          node <- child_node
        } else {
          break
        }
      }
      
      where[i] <- node
    }
  }
  
  # Create result matrix
  result_matrix <- matrix(0, nrow=length(where), ncol=ncol(Y.df))
  colnames(result_matrix) <- colnames(Y.df)
  
  # Extract counts for each observation
  unique_nodes <- unique(where)
  # cat("Number of unique terminal nodes:", length(unique_nodes), "\n")
  
  for (i in 1:length(where)) {
    node_str <- as.character(where[i])
    row_idx <- which(rownames(fit_copy$frame) == node_str)
    
    if (length(row_idx) > 0) {
      if (!is.null(fit_copy$frame$yval2)) {
        # Try to get counts from yval2
        if (is.matrix(fit_copy$frame$yval2) && ncol(fit_copy$frame$yval2) == ncol(Y.df)) {
          result_matrix[i,] <- fit_copy$frame$yval2[row_idx, ]
        }
      } else if ("label" %in% names(fit_copy$frame)) {
        # Try to get label from frame
        node_label <- fit_copy$frame$label[row_idx]
        if (is.list(node_label)) node_label <- unlist(node_label)
        if (length(node_label) == ncol(Y.df)) {
          result_matrix[i, ] <- node_label
        }
      }
    }
  }
  
  # Apply Dirichlet smoothing and return appropriate format
  if (type == "matrix") {
    return(result_matrix)
  }
  
  probs <- t(apply(result_matrix, 1, function(counts) {
    (counts + alphavec) / sum(counts + alphavec)
  }))
  colnames(probs) <- colnames(Y.df)
  
  unique_probs <- nrow(unique(round(probs, 8)))
  # cat("Number of unique probability patterns:", unique_probs, "\n")
  
  if (type == "vector" || type == "prob") {
    return(probs)
  } else if (type == "class") {
    return(apply(probs, 1, which.max))
  }
}

#-----------------------------------------------
# original
#-----------------------------------------------
#----------------- 4) pred() ------------------#
# mypred <- function(fit, newdata = NULL, type = c("vector", "prob", "class", "matrix")) {
#   # 1. argument setup
#   type <- match.arg(type) # ensure type is valid
#   
#   # 2. determine terminal nodes i.e. "where" for each observation 
#   if (is.null(newdata)){
#     # no new data provided, use terminal node assignment from fitting
#     where <- fit$where
#   } else {
#     newdata <- as.data.frame(newdata) 
#     where <- integer(nrow(newdata)) # initialize node assignments
#     
#     for (i in 1:nrow(newdata)){
#       # start at root
#       node <- 1
#       
#       while (TRUE){ # loop until leaf
#         if (fit$frame$var[node] == "<leaf>") break # reached leaf
#         
#         split.var <- as.character(fit$frame$var[node]) # split on this variable
#         split.rows <- which(rownames(fit$splits) == split.var) # find split info
#         
#         if(length(split.rows) == 0) break # no split info: treat as leaf
#         
#         cur.val <- newdata[i, split.var] # get current value 
#         split.row <- split.rows[1] # first split info
#         split.val <- fit$splits[split.row, "index"] # get split value
#         
#         
#         # decide direction based on split
#         go.left <- !is.na(cur.val) && as.numeric(cur.val) <= split.val
#         
#         # traverse to child node of stop if not available 
#         if (go.left){
#           left.child <- 2*node
#           if (left.child %in% as.numeric(rownames(fit$frame))){
#             node <- left.child 
#           } else break # left child doesn't exist
#         } else {
#           right.child <- 2*node + 1
#           if (right.child %in% as.numeric(rownames(fit$frame))){
#             node <- right.child
#           } else break # right child doesn't exist
#         }
#       }
#       
#       where[i] <- node # assign final node
#     }
#   }
#   
#   # 3. construct result matrix w/ label counts from terminal nodes
#   frame <- fit$frame
#   result.matrix <- matrix(0, nrow=length(where), ncol=ncol(Y.df))
#   colnames(result.matrix) <- colnames(Y.df)
#   
#   for (i in 1:length(where)) { # for each observation 
#     node.idx <- where[i]
#     
#     if ("label" %in% names(frame)){ # try to get label from frame
#       
#       node.label <- frame$label[node.idx]
#       
#       if (is.list(node.label)) node.label <- unlist(node.label) # unlist
#       
#       if (length(node.label) == ncol(Y.df)) {
#         result.matrix[i, ] <- node.label # assign label counts
#       }
#       
#     } else if (!is.null(frame$yval2) && ncol(frame$yval2) == ncol(Y.df)) { # try to get yval2 for counts
#       # if no label, try yval2
#       result.matrix[i, ] <- frame$yval2[node.idx, ]
#     }
#   }
#   
#   # 4. convert counts to probabilities or predicted class
#   if (type == "matrix") {
#     return(result.matrix) # return raw counts
#   }
#   
#   # use dirichlet smoothing to convert counts to probabilities
#   probs <- t(apply(result.matrix, 1, function(counts){
#     (counts + alphavec) / sum(counts + alphavec)
#   }))
#   colnames(probs) <- colnames(Y.df)
#   
#   # 5. return based on type 
#   if (type == "vector" || type == "prob") {
#     return(probs) # return smoothed probabilities
#   } else if (type == "class") {
#     result <- apply(probs, 1, which.max) # return index of highest prob birth weight category
#     return(result)
#   }
#   
#   
# }

dm.method <- list(init=myinit, eval=myeval, split=mysplit, pred=mypred, method="dm")
# pred <- mypred(dm.tree, newdata = data.frame(X.matrix), type = "prob")


#-----------------------------------------------------------
# D. Fit rpart Tree Using Our Custom DM Method
#-----------------------------------------------------------
cat("\n--- Fitting the DM-based rpart tree ---\n")

dm.control <- rpart.control(minsplit=2, cp=0, maxdepth=8, xval=0, usesurrogate = 0)
dm.tree <- rpart(
  formula = Y.df~.,
  data = data.frame(X.matrix),
  method = dm.method,
  control = dm.control
)


# preds <- mypred(dm.tree, newdata = data.frame(X.matrix), type = "prob")
# head(preds, 10)

# write.csv(preds, file = sprintf("/Users/adamkurth/Downloads/preds_%d.csv", year), row.names = FALSE)
# dim(mypreds) # 128 x 11
#   - hits "where <- fit$where" and outputs pred. prob. vector for 128 rows
#   - retrieves terminal node indices for each row
#   - "result.matrix" builds one row per terminal node (i.e. 128 total)
#   - "frame$label[node.idx]" contains label counts at leaf
#   - each i-th row in result.matrix represents raw counts of label occurrences at term. node that i-th obs falls in.
#    -- i.e. node 12 has 3 examples w/ labels [1,0,0,...][1,1,0,...][0,1,0,...] then counts might look like [2,2,0,...]
#   - "probs" smooths the counts into prob. using alphavec, sum to 1
#   - sum of each row = 1, each column = category


#-----------------------------------------------------------
# E. Save and Inspect Results
#-----------------------------------------------------------
save(dm.tree, file=sprintf("%s/dm_tree_rebin_%d.RData", data.pwd, year))
print(dm.tree)
printcp(dm.tree)

plot(dm.tree, main="DM Tree",uniform=TRUE)
text(dm.tree,use.n=TRUE,cex=1,font=3)

path.rpart(dm.tree, node = 2)     # left child  => mrace15 = 0
path.rpart(dm.tree, node = 3)     # right child => mrace15 = 1

#-----------------------------------------------------------
# F. Multinomial Parametric Bootstrap Sampling
#-----------------------------------------------------------

B <- 10000  # Number of bootstrap samples
# B <- 100
n.rows <- nrow(counts.df)

# category definitions
get.category.cols <- function(type = "lbw"){

    if (type == "lbw") {
      return(1:10)  # Indices of LBW categories 
    } else if (type == "nbw") {
      return(11:ncol(Y.df))  # Adjust based on your data structure
    } else if (type == "all") {
      return(1:ncol(Y.df))  # All categories
    } else {
      stop("Invalid category type. Use 'lbw', 'nbw', or 'all'")
    }

}

# define category cols 
lbw.cols <- get.category.cols(type = "lbw")
nbw.cols <- get.category.cols(type = "nbw")

# create one structure to store all probability vectors 
# each element in list is bootstrap sample, of matrix n.rows x ncol(Y.df)
Yhat <- vector("list", length = B)


# Dirichlet-smoothed cell probabilities
multinomial.probs <- matrix(0, nrow = n.rows, ncol = ncol(Y.df))

for(i in 1:n.rows) {
  
  # Get counts for this cell/row
  cell.counts <- as.numeric(Y.df[i,])
  
  # Calculate probability with Dirichlet smoothing
  cell.probs <- (cell.counts + alphavec) / (sum(cell.counts) + sum(alphavec))
  
  multinomial.probs[i, ] <- cell.probs
}


# row-probability vector for total-count bootstrap
row.probs <- rowMeans(counts.df) # positive weights
row.probs <- row.probs / sum(row.probs)  # Normalize to sum to 1

barplot(row.probs, main="Row Probabilities for Bootstrap Sampling", 
        xlab="Row Index", ylab="Probability")

# initialize list for tree structure information
tree.structure <- list()
top.vars.used <- character(B)
n.vars.used <- integer(B)

cat("Multinomial parametric bootstrap     (inter-/intra- resampling)\n")
pb <- txtProgressBar(0, B, style = 3)


#-----------------------------------------------------------
# G. Main Loop over B Bootstrap Samples
#-----------------------------------------------------------
for (b in seq_len(B)) {
    # Note: Y.df == counts.df
    #---------------------------------------
    # (A) 2-stage parametric bootstrap of the 128×11 count table
    #---------------------------------------
    # 1st stage: generate new counts from multinom

    n.star <- as.vector( rmultinom(1, size = sum(Y.df[]), prob = row.probs) )
    
    # 2nd stage: generate counts for each row using multinom
    boot.counts <- t(
      sapply(1:n.rows, function(j) rmultinom(1, size = n.star[j], prob = multinomial.probs[j, ]))
    )
    
    colnames(boot.counts) <- colnames(Y.df)
    boot.df <- as.data.frame(boot.counts)

    
    #---------------------------------------------------------
    # (B) fit the full DM tree on the bootstrap counts
    #---------------------------------------------------------
    model.data <- cbind(boot.df, X.matrix)
    r.formula <- as.formula( paste0("cbind(", paste(colnames(boot.df), collapse = ","), ") ~ .") )
    
    # fit the tree on the bootstrap sample
    dm.tree.b <- rpart(
      formula = r.formula,
      data = model.data,
      method  = dm.method,
      control = dm.control
    )
    
    #---------------------------------------------------------
    # (C) record tree statistics
    #--------------------------------------------------------
    vars.used <- unique(dm.tree.b$frame$var[dm.tree.b$frame$var != "<leaf>"])
    n.vars.used[b] <- length(vars.used)
    top.vars.used[b] <- if (length(vars.used)) vars.used[1] else "none" # store first split variable
    
    tree.structure[[b]] <-list(
      variables = vars.used,
      frame = dm.tree.b$frame,
      splits = dm.tree.b$splits
    )


    #---------------------------------------------------------
    # (D) Get predicted probability vectors
    #--------------------------------------------------------
    # custom predict function to get probability vectors
    all.preds <- mypred(dm.tree.b, newdata = data.frame(X.matrix), type = "prob")
    
    # Store the complete probability matrix for this bootstrap sample
    Yhat[[b]] <- all.preds
  
    
    setTxtProgressBar(pb, b)
}
close(pb)
cat("\n--- Bootstrap sampling completed ---\n")


#-----------------------------------------------------------
# Visualize the bootstrap results
#-----------------------------------------------------------
is.high.risk.indices <- which(
  X.matrix$mrace15 == 1 &
    X.matrix$dmar == 0 & 
    X.matrix$cig_0 == 1 &
    X.matrix$sex == 0 &  # female 
    X.matrix$mager == 0 &  # young
    X.matrix$precare5 == 0 &  # not adequate prenatal
    X.matrix$meduc == 0  # not high school grad
)

is.low.risk.indices <- which(
  X.matrix$mrace15 == 0 &
    X.matrix$dmar == 1 & 
    X.matrix$cig_0 == 0 &
    X.matrix$sex == 1 &  # male 
    X.matrix$mager == 1 &  # older
    X.matrix$precare5 == 0 &  # adequate prenatal
    X.matrix$meduc == 1  # high school grad
)


direct.bootstrap.calculation <- function(indicies, yhat.list = Yhat) {
  B <- length(yhat.list)
  K <- ncol(yhat.list[[1]])
  
  # for each b and category k, calculate the mean probability across B
  prob.matrix <- matrix(0, nrow = B, ncol = K)
  
  for (b in 1:B) {
    if (length(indicies) > 0) {
      # Get probabilities for the selected indices in this bootstrap sample
      prob.matrix[b, ] <- colMeans(yhat.list[[b]][indicies, , drop = FALSE])
    }
  }
  
  # mean probabilities by category
  mean.probs <- colMeans(prob.matrix)
  
  # 95% percentile bootstrap conf int
  ci.lwr <- apply(prob.matrix, 2, function(x) quantile(x, 0.025))
  ci.upr <- apply(prob.matrix, 2, function(x) quantile(x, 0.975))

  return(list(
    means = mean.probs,
    lwr = ci.lwr,
    upr = ci.upr
  ))
}


high.risk.result <- direct.bootstrap.calculation(indicies = is.high.risk.indices, yhat.list = Yhat)
high.risk.probs <- high.risk.result$means
high.risk.lwr <- high.risk.result$lwr
high.risk.upr <- high.risk.result$upr

low.risk.result <- direct.bootstrap.calculation(indicies = is.low.risk.indices, yhat.list = Yhat)
low.risk.probs <- low.risk.result$means
low.risk.lwr <- low.risk.result$lwr
low.risk.upr <- low.risk.result$upr

K <- length(high.risk.probs)



# 
if (type == 1) {
  result.tab <- read.csv(sprintf("%s/result_table_1.csv", results.pwd))
  result.tab[-1]
} else {
  result.tab <- read.csv(sprintf("%s/result_table_2.csv", results.pwd))
  result.tab[-1]
}

# 
result.tab.1 <- read.csv("~/Downloads/result_table_1.csv")[-1]
result.tab.2 <- read.csv("~/Downloads/result_table_2.csv")[-1]
diff.pi.hat <- result.tab.1$high.risk.prob[1:10] - result.tab.1$low.risk.prob[1:10]



cat.labels <- paste0("C", 1:K)
pdf(file = sprintf("%s/high_low_risk_pred_%d.pdf", boot.pwd, year), width = 8.5, height = 8)
par(mfrow = c(2, 1), 
    mar = c(3, 4.5, 3, 2),    # More space around individual plots
    oma = c(2, 0, 4, 0))      # Larger outer margin at top for title and subtitle

# top
high.bars <- barplot(high.risk.probs, 
                     main = "",
                     xlab = "", 
                     ylab = "Probability",
                     ylim = c(0, max(high.risk.upr) * 1.1),
                     names.arg = cat.labels,  
                     col = "red",
                     cex.names = 1.0,
                     cex.axis = 1.0,
                     border = "darkred")  # Add border for better definition
title("High Risk Subgroup", line = 1, cex.main = 1.2)

arrows(high.bars, high.risk.lwr, 
       high.bars, high.risk.upr, 
       angle = 90, code = 3, length = 0.05, 
       col = "black", lwd = 1.5)


# bottom
low.bars <- barplot(low.risk.probs,
                    main = "",  # Remove individual title
                    xlab = "", 
                    ylab = "Probability",
                    ylim = c(0, max(high.risk.upr) * 1.1),
                    names.arg = cat.labels,
                    col = "blue",
                    cex.names = 1.0,
                    cex.axis = 1.0,
                    border = "darkblue")  # Add border

title("Low Risk Subgroup", line = 1, cex.main = 1.2)

arrows(low.bars, low.risk.lwr, 
       low.bars, low.risk.upr, 
       angle = 90, code = 3, length = 0.05, 
       col = "black", lwd = 1.5)

mtext(sprintf("Birth Weight Category Probabilities (%d)", year), 
      side = 3, line = 2, outer = TRUE, cex = 1.4, font = 2)

# Add subtitle
mtext("Predicted probabilities with 95% confidence intervals", 
      side = 3, line = 0.5, outer = TRUE, cex = 1.1)

# Add note about categories if needed
mtext("Categories represent different birth weight ranges", 
      side = 1, line = 0, outer = TRUE, cex = 0.9, col = "darkgray")

dev.off()
par(mfrow = c(1, 1))




#-----------------------------------------------------------
# Quick check for pred. prob. vectors 
#-----------------------------------------------------------
analyze.bootstrap.samples <- function(Y.hat) {
  B <- length(Y.hat)
  
  lbw.cols <- get.category.cols(type = "lbw")
  nbw.cols <- get.category.cols(type = "nbw")
  
  # --- Category-level variation ---
  category.stats <- lapply(1:ncol(Y.hat[[1]]), function(i) {
    probs <- sapply(Y.hat, function(b) b[, i])
    list(
      mean = rowMeans(probs),
      sd = apply(probs, 1, sd),
      cv = apply(probs, 1, sd) / rowMeans(probs)
    )
  })
  
  # --- LBW probability analysis ---
  lbw.probs <- sapply(Y.hat, function(b) rowSums(b[, lbw.cols]))
  lbw.means <- rowMeans(lbw.probs)
  lbw.sds <- apply(lbw.probs, 1, sd)
  lbw.cvs <- lbw.sds / lbw.means
  
  # # --- NBW probability analysis ---
  # nbw.probs <- rows)
  # nbw.means <- rowMeans(nbw.probs)
  # nbw.sds <- apply(nbw.probs, 1, sd)
  # nbw.cvs <- nbw.sds / nbw.means
  
  cat("LBW mean prob:", mean(lbw.means),
      "| SD:", sd(lbw.means),
      "| Mean CV:", mean(lbw.cvs), "\n")
  
  # cat("NBW mean prob:", mean(nbw.means),
  #     "| SD:", sd(nbw.means),
  #     "| Mean CV:", mean(nbw.cvs), "\n")
  
  # --- Category proportion consistency ---
  category.props <- sapply(Y.hat, function(b) colMeans(b))
  category.means <- rowMeans(category.props)
  category.sds <- apply(category.props, 1, sd)
  category.cvs <- category.sds / category.means
  
  proportion.summary <- data.frame(
    # category = colnames(Y.hat[[1]]),
    mean.prob = category.means,
    sd = category.sds,
    cv = category.cvs
  )
  print(proportion.summary)
  
  # --- Terminal node pattern consistency ---
  pattern.counts <- sapply(Y.hat, function(b) {
    length(unique(apply(round(b, 6), 1, paste, collapse = ",")))
  })
  cat("Unique terminal node patterns per sample:\n")
  print(summary(pattern.counts))
  
  list(
    category.stats = category.stats,
    lbw.stats = list(means = lbw.means, sds = lbw.sds, cvs = lbw.cvs),
    category.proportions = proportion.summary,
    pattern.counts = pattern.counts
  )
}

compare.category.predictions <- function(Y.hat, category.idx) {
  B <- length(Y.hat)
  category.preds <- sapply(Y.hat, function(b) b[, category.idx])
  
  means <- rowMeans(category.preds)
  sds <- apply(category.preds, 1, sd)
  cvs <- sds / means
  
  par(mfrow = c(3, 1))
  
  hist(means,
       main = paste("Distribution of Mean Probabilities -", colnames(Y.hat[[1]])[category.idx]),
       xlab = "Mean Probability", col = "lightblue", border = "white")
  
  hist(sds,
       main = "Distribution of Standard Deviations",
       xlab = "Standard Deviation", col = "lightblue", border = "white")
  
  hist(cvs,
       main = "Distribution of Coefficients of Variation",
       xlab = "CV", col = "lightblue", border = "white")
  # 
  # plot(means, sds,
  #      main = "Standard Deviation vs. Mean",
  #      xlab = "Mean Probability", ylab = "Standard Deviation",
  #      pch = 19, col = "blue")
  # 
  par(mfrow = c(1, 1))
  
  list(
    means = means,
    sds = sds,
    cvs = cvs,
    overall.mean = mean(means),
    overall.sd = sd(means),
    mean.cv = mean(cvs)
  )
}

validate.full.probability.vectors <- function(Yhat) {
  validations <- lapply(seq_along(Yhat), function(i) {
    probs <- Yhat[[i]]
    
    row.sums <- rowSums(probs)
    all.rows.sum.to.1 <- all(abs(row.sums - 1) < .Machine$double.eps^0.5)
    all.values.in.01 <- all(probs >= 0 & probs <= 1)
    
    if (!all.rows.sum.to.1) {
      cat(sprintf("Warning (bootstrap %d): Not all rows sum to 1. Min: %.5f | Max: %.5f\n",
                  i, min(row.sums), max(row.sums)))
    }
    if (!all.values.in.01) {
      cat(sprintf("Warning (bootstrap %d): Some values outside [0,1] range.\n", i))
    }
    
    faulty.rows <- which(abs(row.sums - 1) > .Machine$double.eps^0.5 |
                           probs < 0 | probs > 1, arr.ind = TRUE)
    
    list(
      bootstrap.index = i,
      all.rows.sum.to.1 = all.rows.sum.to.1,
      all.values.in.01 = all.values.in.01,
      faulty.rows = faulty.rows
    )
  })
  
  return(validations)
}

bootstrap.analysis <- analyze.bootstrap.samples(Yhat)
print(bootstrap.analysis$category.proportions)

cat.1.analysis <- compare.category.predictions(Yhat, 1)
cat.6.analysis <- compare.category.predictions(Yhat, 6)
cat.10.analysis <- compare.category.predictions(Yhat, 10)
# cat.11.analysis <- compare.category.predictions(Yhat, 11)

# category.comparison <- data.frame(
#   category = c(colnames(Yhat[[1]])[1], colnames(Yhat[[1]])[6], colnames(Yhat[[1]])[11]),
#   mean.prob = c(cat.1.analysis$overall.mean, cat.6.analysis$overall.mean, cat.11.analysis$overall.mean),
#   sd.prob = c(cat.1.analysis$overall.sd, cat.6.analysis$overall.sd, cat.11.analysis$overall.sd),
#   mean.cv = c(cat.1.analysis$mean.cv, cat.6.analysis$mean.cv, cat.11.analysis$mean.cv)
# )
category.comparison <- data.frame(
  category = c(colnames(Yhat[[100]])[1], colnames(Yhat[[100]])[6], colnames(Yhat[[100]])[10]),
  mean.prob = c(cat.1.analysis$overall.mean, cat.6.analysis$overall.mean, cat.10.analysis$overall.mean),
  sd.prob = c(cat.1.analysis$overall.sd, cat.6.analysis$overall.sd, cat.10.analysis$overall.sd),
  mean.cv = c(cat.1.analysis$mean.cv, cat.6.analysis$mean.cv, cat.10.analysis$mean.cv)
)

print(category.comparison)

hist(bootstrap.analysis$category.proportions$mean.prob,
     main = "Distribution of Mean Predicted Probabilities",
     xlab = "Mean Probability", col = "lightblue")

vector.validations <- validate.full.probability.vectors(Yhat)
print(vector.validations[1:5])





compare.probability.matrices <- function(Yhat, alphavec, return.full.diff = FALSE) {
  comparisons <- lapply(seq_along(Yhat), function(i) {
    probs <- Yhat[[i]]
    
    if (!all(dim(probs) == dim(alphavec))) {
      stop(sprintf("Dimension mismatch in bootstrap %d: Yhat has dim %s but alphavec has dim %s",
                   i,
                   paste(dim(probs), collapse = "x"),
                   paste(dim(alphavec), collapse = "x")))
    }
    
    diff.matrix <- probs - alphavec
    max.abs.diff <- max(abs(diff.matrix))
    mean.abs.diff <- mean(abs(diff.matrix))
    
    list(
      bootstrap.index = i,
      max.abs.diff = max.abs.diff,
      mean.abs.diff = mean.abs.diff,
      diff.matrix = if (return.full.diff) diff.matrix else NULL
    )
  })
  
  return(comparisons)
}


comparison.results <- compare.probability.matrices(Yhat, alphavec)
str(comparison.results)




#-----------------------------------------------------------
# Depth distributions visualization
#-----------------------------------------------------------
library(reshape2)
library(dplyr)
library(ggplot2)

# get all possible predictor variables
all.vars <- colnames(X.matrix)
n.vars <- length(all.vars)
co.occur.mat <- matrix(0, nrow = n.vars, ncol = n.vars, dimnames = list(all.vars, all.vars))
depth.sums <- setNames(numeric(n.vars), all.vars)
depth.counts <- setNames(numeric(n.vars), all.vars)
depth.values <- vector("list", length = n.vars)
names(depth.values) <- all.vars

for (b in seq_len(B)){ 
  tree <- tree.structure[[b]]
  frame <- tree$frame
  vars.used <- tree$variables
  
  # update co-coccurence matrix
  if (length(vars.used) >= 1) {
      idx <- which(all.vars %in% vars.used)
      co.occur.mat[idx, idx] <- co.occur.mat[idx, idx] + 1
  }
  
  # compute depth per node
  node.ids <- as.numeric(rownames(frame))
  node.depths <- floor(log2(node.ids))
  vars <- frame$var
  internal.nodes <- which(vars != "<leaf>")
  
    for (i in internal.nodes) {
      var <- vars[i]
      depth <- node.depths[i]
        if (var %in% names(depth.sums)) {
            depth.sums[var] <- depth.sums[var] + depth
            depth.counts[var] <- depth.counts[var] + 1
            depth.values[[var]] <- c(depth.values[[var]], depth)
        }
    }
}    

co.occur.prop <- co.occur.mat / B
avg.depths <- depth.sums / depth.counts
depth.df <- data.frame(
  Variable = names(avg.depths),
  AvgDepth = avg.depths
)

depth.long <- bind_rows(
  lapply(names(depth.values), function(var) {
    depths <- depth.values[[var]]
    if (length(depths) > 0) {
      data.frame(Variable = var, Depth = depths)
    }
  })
)
depth.long$Variable <- factor(depth.long$Variable, levels = depth.df$Variable[order(depth.df$AvgDepth)])

means <- depth.long %>%
  group_by(Variable) %>%
  summarise(mean_depth = mean(Depth))

ggplot(depth.long, aes(x = Depth, fill = Variable)) +
  geom_histogram(binwidth = 1, alpha = 0.6) +
  geom_vline(data = means, aes(xintercept = mean_depth),
             color = "black", linetype = "dashed", size = 0.5) +
  facet_wrap(~ Variable, scales = "free_y") +
  labs(title = "Distribution of Depths for Each Predictor Variable",
       x = "Depth",
       y = "Number of Splits") +
  theme_minimal()

ggsave(sprintf("%s/depth_distributions_%d.png", boot.pwd, type), width = 12, height = 8, dpi = 400)

#-----------------------------------------------------------
# H. Decision Tree Depth Analysis for cig_0 Predictor
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
# J) Comparative Analysis of Trees
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

#-----------------------------------------------------------
# K) Analyze Bootstrap Tree Structure Results
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


bootstrap.tree.results <- list(
  var.usage = var.usage,
  var.freq.df = var.freq.df,
  top.var.df = top.var.df,
  n.vars.summary = n.vars.summary
)

#-----------------------------------------------------------
cat("\n--- Bootstrap Tree Structure Analysis ---\n")
cat("\nFrequency of Variables Used in Trees:\n")
print(var.freq.df)

cat("\nFrequency of Top Split Variables:\n")
print(top.var.freq)

cat("\nSummary of Number of Variables Used:\n")
print(n.vars.summary)

save(bootstrap.tree.results,file = sprintf("%s/bootstrap_tree_results_%d.RData", results.pwd, year))