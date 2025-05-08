rm(list=ls())
library(rpart)
library(ggplot2)
library(xtable)

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
  
  # If no new data, use data from the fit
  if (is.null(newdata)) {
    where <- fit$where
  } else {
    # Ensure newdata has the right format
    newdata <- as.data.frame(newdata)
    
    # Custom tree traversal to find terminal nodes
    where <- integer(nrow(newdata))
    
    # For each observation in newdata
    for (i in 1:nrow(newdata)) {
      # Start at root node
      node <- 1
      
      # Traverse until reaching a leaf or no valid split found
      while (TRUE) {
        # Check if current node is a leaf
        if (fit$frame$var[node] == "<leaf>") {
          break
        }
        
        # Get split variable name
        split_var <- as.character(fit$frame$var[node])
        
        # Find split information in splits matrix
        split_rows <- which(rownames(fit$splits) == split_var)
        
        if (length(split_rows) == 0) {
          # No split information available, treat as leaf
          break
        }
        
        # Get current value for this variable
        cur_val <- newdata[i, split_var]
        
        # Find which row in splits to use
        # For binary variables, use the first split info
        split_row <- split_rows[1]
        split_val <- fit$splits[split_row, "index"]
        
        # Go left if value <= threshold
        go_left <- !is.na(cur_val) && as.numeric(cur_val) <= split_val
        
        # Move to next node
        if (go_left) {
          # Go left
          left_child <- 2 * node
          if (left_child %in% as.numeric(rownames(fit$frame))) {
            node <- left_child
          } else {
            # Left child doesn't exist, stay at current node
            break
          }
        } else {
          # Go right
          right_child <- 2 * node + 1
          if (right_child %in% as.numeric(rownames(fit$frame))) {
            node <- right_child
          } else {
            # Right child doesn't exist, stay at current node
            break
          }
        }
      }
      
      where[i] <- node
    }
  }
  
  # Frame data
  frame <- fit$frame
  
  # Initialize result matrix
  result_matrix <- matrix(0, nrow = length(where), ncol = ncol(Y.df))
  colnames(result_matrix) <- colnames(Y.df)
  
  # Get the label values (counts) for each node
  for (i in 1:length(where)) {
    node_idx <- where[i]
    
    # Try to get label from the frame
    if ("label" %in% names(frame)) {
      node_label <- frame$label[node_idx]
      
      # Handle case where label is a list
      if (is.list(node_label)) {
        node_label <- unlist(node_label)
      }
      
      if (length(node_label) == ncol(Y.df)) {
        result_matrix[i, ] <- node_label
      }
    } else if (!is.null(frame$yval2) && ncol(frame$yval2) == ncol(Y.df)) {
      # If no label, try yval2
      result_matrix[i, ] <- frame$yval2[node_idx, ]
    }
  }
  
  # Return appropriate output based on type
  if (type == "matrix") {
    # Return raw counts
    return(result_matrix)
  } else if (type == "vector" || type == "prob") {
    # Convert counts to probabilities using Dirichlet smoothing
    probs <- t(apply(result_matrix, 1, function(counts) {
      (counts + alphavec) / sum(counts + alphavec)
    }))
    colnames(probs) <- colnames(Y.df)
    return(probs)
  } else if (type == "class") {
    # Return index of highest probability category
    probs <- t(apply(result_matrix, 1, function(counts) {
      (counts + alphavec) / sum(counts + alphavec)
    }))
    result <- apply(probs, 1, which.max)
    return(result)
  }
}




dm.method <- list(init=myinit, eval=myeval, split=mysplit, pred=mypred, method="dm")
# dm.method <- list(init=myinit, eval=myeval, split=mysplit, method="dm")

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


preds <- mypred(dm.tree, newdata = data.frame(X.matrix), type = "prob")
# head(preds)  # Check the first few predictions

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

# B <- 10000  # Number of bootstrap samples
B <- 100
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
    
    # alternative implementation using vapply
    # boot.counts <- t( vapply(seq_len(n.rows), function(j)
    #   rmultinom(1, n.star[j], multinomial.probs[j, ]),
    #   integer(ncol(counts.df)))
    # ) 
    
    colnames(boot.counts) <- colnames(counts.df)
    boot.df <- as.data.frame(boot.counts)
  
    #---------------------------------------------------------
    # (B) fit the full DM tree on the bootstrap counts
    #---------------------------------------------------------
    model.data <- cbind(boot.df, X.matrix)
    r.formula <- as.formula( paste0("cbind(", paste(colnames(boot.df), collapse = ","), ") ~ .") )
    
    
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










analyze_bootstrap_samples <- function(Yhat, lbw.cols) {
  B <- length(Yhat)  # Number of bootstrap samples
  
  # 1. Analyze probability variation across bootstrap samples
  cat("\n===== Probability Variation Across Bootstrap Samples =====\n")
  
  # Function to analyze a specific category across all bootstrap samples
  analyze_category <- function(cat_idx) {
    # Extract probabilities for this category from all bootstrap samples
    cat_probs <- sapply(1:B, function(b) {
      # Extract column for this category from each bootstrap sample
      Yhat[[b]][, cat_idx]
    })
    
    # Calculate statistics
    mean_probs <- rowMeans(cat_probs)
    sd_probs <- apply(cat_probs, 1, sd)
    cv_probs <- sd_probs / mean_probs
    
    return(list(
      mean = mean_probs,
      sd = sd_probs,
      cv = cv_probs,
      min = apply(cat_probs, 1, min),
      max = apply(cat_probs, 1, max)
    ))
  }
  
  # Analyze all categories
  category_stats <- lapply(1:ncol(Yhat[[1]]), analyze_category)
  
  # Summarize results
  cat("Category-level statistics across bootstrap samples:\n")
  for (i in 1:length(category_stats)) {
    cat("Category", i, "- Mean CV:", mean(category_stats[[i]]$cv), 
        ", Max CV:", max(category_stats[[i]]$cv), "\n")
  }
  
  # 2. Calculate LBW probability for each bootstrap sample
  cat("\n===== LBW Probability Analysis =====\n")
  
  # Calculate LBW probability for each observation across bootstrap samples
  lbw_probs <- matrix(0, nrow = nrow(Yhat[[1]]), ncol = B)
  for (b in 1:B) {
    lbw_probs[, b] <- rowSums(Yhat[[b]][, lbw.cols, drop = FALSE])
  }
  
  # Calculate statistics
  lbw_means <- rowMeans(lbw_probs)
  lbw_sds <- apply(lbw_probs, 1, sd)
  lbw_cvs <- lbw_sds / lbw_means
  
  # Summarize results
  cat("LBW probability statistics:\n")
  cat("Mean LBW probability:", mean(lbw_means), "\n")
  cat("SD of mean LBW probabilities:", sd(lbw_means), "\n")
  cat("Mean CV of LBW probabilities:", mean(lbw_cvs), "\n")
  cat("Number of unique mean LBW probabilities:", length(unique(round(lbw_means, 4))), "\n")
  
  # 3. Check consistency of category proportions
  cat("\n===== Category Proportion Consistency =====\n")
  
  # Calculate category proportions for each bootstrap sample
  cat_props <- matrix(0, nrow = B, ncol = ncol(Yhat[[1]]))
  colnames(cat_props) <- colnames(Yhat[[1]])
  
  for (b in 1:B) {
    # Calculate mean probability for each category
    cat_props[b, ] <- colMeans(Yhat[[b]])
  }
  
  # Calculate statistics
  cat_prop_means <- colMeans(cat_props)
  cat_prop_sds <- apply(cat_props, 2, sd)
  cat_prop_cvs <- cat_prop_sds / cat_prop_means
  
  # Create summary data frame
  prop_summary <- data.frame(
    category = colnames(cat_props),
    mean_prob = cat_prop_means,
    sd = cat_prop_sds,
    cv = cat_prop_cvs
  )
  
  # Print summary
  print(prop_summary)
  
  # 4. Analyze terminal node assignment consistency
  cat("\n===== Terminal Node Assignment Consistency =====\n")
  
  # Check if we can identify terminal nodes from probability patterns
  # Create a function to identify unique probability vectors
  identify_patterns <- function(b) {
    # Round probabilities to reduce numerical noise
    rounded_probs <- round(Yhat[[b]], 6)
    
    # Identify unique patterns
    unique_patterns <- unique(apply(rounded_probs, 1, paste, collapse = ","))
    
    return(list(
      num_patterns = length(unique_patterns),
      patterns = unique_patterns
    ))
  }
  
  # Analyze patterns across bootstrap samples
  pattern_counts <- sapply(1:B, function(b) identify_patterns(b)$num_patterns)
  
  # Summarize results
  cat("Number of unique probability vectors per bootstrap sample:\n")
  print(summary(pattern_counts))
  
  # Return the comprehensive analysis results
  return(list(
    category_stats = category_stats,
    lbw_stats = list(means = lbw_means, sds = lbw_sds, cvs = lbw_cvs),
    category_proportions = prop_summary,
    pattern_counts = pattern_counts
  ))
}
compare_category_predictions <- function(Yhat, category_idx) {
  # Get total number of bootstrap samples
  B <- length(Yhat)
  
  # Extract predictions for this category across all bootstrap samples
  category_preds <- matrix(0, nrow = nrow(Yhat[[1]]), ncol = B)
  for (b in 1:B) {
    category_preds[, b] <- Yhat[[b]][, category_idx]
  }
  
  # Calculate statistics
  means <- rowMeans(category_preds)
  sds <- apply(category_preds, 1, sd)
  cvs <- sds / means
  
  # Plot results
  par(mfrow = c(2, 2))
  
  # Plot 1: Distribution of mean probabilities
  hist(means, 
       main = paste("Distribution of Mean Probabilities -", colnames(Yhat[[1]])[category_idx]),
       xlab = "Mean Probability", 
       col = "lightblue", border = "white")
  
  # Plot 2: Distribution of standard deviations
  hist(sds, 
       main = "Distribution of Standard Deviations",
       xlab = "Standard Deviation", 
       col = "lightblue", border = "white")
  
  # Plot 3: Distribution of CVs
  hist(cvs, 
       main = "Distribution of Coefficients of Variation",
       xlab = "CV", 
       col = "lightblue", border = "white")
  
  # Plot 4: Scatter plot of SD vs. Mean
  plot(means, sds, 
       main = "Standard Deviation vs. Mean",
       xlab = "Mean Probability", ylab = "Standard Deviation",
       pch = 19, col = "blue")
  
  # Reset plot parameters
  par(mfrow = c(1, 1))
  
  # Return statistics
  return(list(
    means = means,
    sds = sds,
    cvs = cvs,
    overall_mean = mean(means),
    overall_sd = sd(means),
    mean_cv = mean(cvs)
  ))
}

# Run for each category (example for categories 1, 6, and 11)
cat1_analysis <- compare_category_predictions(Yhat, 1) # First category
cat6_analysis <- compare_category_predictions(Yhat, 6) # Middle category
cat11_analysis <- compare_category_predictions(Yhat, 11) # Last category (above_2.5kg)

# Compare statistics across categories
category_comparison <- data.frame(
  category = c(colnames(Yhat[[1]])[1], colnames(Yhat[[1]])[6], colnames(Yhat[[1]])[11]),
  mean_prob = c(cat1_analysis$overall_mean, cat6_analysis$overall_mean, cat11_analysis$overall_mean),
  sd_prob = c(cat1_analysis$overall_sd, cat6_analysis$overall_sd, cat11_analysis$overall_sd),
  mean_cv = c(cat1_analysis$mean_cv, cat6_analysis$mean_cv, cat11_analysis$mean_cv)
)

print(category_comparison)
# Run the comprehensive analysis
bootstrap_analysis <- analyze_bootstrap_samples(Yhat, lbw.cols)

# Create visualizations to show results
par(mfrow = c(2, 1))

# all category mean predicted probabilities
hist(bootstrap_analysis$category_stats[[1]],
     main = "Distribution of Mean Predicted Probabilities",
     xlab = "Mean Probability", ylab = "Frequency",
     col = "lightblue", border = "white")

# 2. Plot CV of category proportions
barplot(bootstrap_analysis$category_proportions$cv, 
        main = "Coefficient of Variation by Category",
        xlab = "Category", ylab = "CV",
        col = "lightblue",
        names.arg = 1:nrow(bootstrap_analysis$category_proportions), 
        las = 2, cex.names = 0.7)

# Reset plotting parameters
par(mfrow = c(1, 1))


# 1. Validate full probability vectors
cat("\n===== Full Probability Vector Validation =====\n")
vector_validation <- validate_full_probability_vectors(dm.tree.b, data.frame(X.matrix))
print(vector_validation)

# 4. Compare with original data
cat("\n===== Comparison with Original Data =====\n")
data_comparison <- compare_with_original_data(dm.tree.b, data.frame(X.matrix), Y.df)
print(data_comparison)






















#-----------------------------------------------------------
# Depth distributions visualization

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

#-----------------------------------------------------------
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

### USE table-boot-freq.R to generate tables from this data





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