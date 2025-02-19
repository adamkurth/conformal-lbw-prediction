# plot rowwise count data
year <- 2021
data.rebin.pwd <- "/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/birthweight_data/rebin"
counts.df <- read.csv(sprintf("%s/counts_df_rebin_%d.csv", data.rebin.pwd, year))

plot.rowwise.intervals <- function(counts.df, save.path = "row_wise_plots_rebin.pdf") {
  # 1. Extract relevant columns and prepare data
  counts.sub <- counts.df[, 1:(ncol(counts.df)-1)]
  
  # 2. Clean up column names (remove prefixes, suffixes, etc.)
  #    and parse the numeric boundaries.
  intervals <- colnames(counts.sub)
  
  # A helper function to strip out "count_btw_" / "count_above_" and "kg"
  clean_interval_name <- function(x) {
    x <- gsub("^count_btw_", "", x)
    x <- gsub("^count_above_", "", x)
    x <- gsub("kg$", "", x)
    x
  }
  
  cleaned_intervals <- sapply(intervals, clean_interval_name)
  # cleaned_intervals should now look like "0_0.1", "0.1_0.2", ..., "2.5"
  
  # A helper function to split the string on "_" and convert to numeric.
  # If there's only one numeric value, we interpret it as "above X kg."
  parse_interval <- function(x) {
    parts <- strsplit(x, "_")[[1]]
    if (length(parts) == 1) {
      # e.g., "2.5" means above 2.5
      left <- as.numeric(parts[1])
      right <- Inf
    } else {
      # e.g., c("0", "0.1")
      left <- as.numeric(parts[1])
      right <- as.numeric(parts[2])
    }
    c(left, right)
  }
  
  # Build a numeric matrix of left/right boundaries
  interval.bounds <- t(sapply(cleaned_intervals, parse_interval))
  colnames(interval.bounds) <- c("left", "right")
  
  # 3. Create a more readable label for each bar
  #    If the right boundary is Inf, label as "2.5+"
  custom.abbrev <- ifelse(
    is.infinite(interval.bounds[, "right"]),
    paste0(interval.bounds[, "left"], "+"), 
    paste0(interval.bounds[, "left"], "-", interval.bounds[, "right"])
  )
  
  # 4. Create output directory if needed
  if (!dir.exists(dirname(save.path))) {
    dir.create(dirname(save.path), recursive = TRUE)
  }
  
  # 5. Configure PDF output
  pdf(save.path, width = 16, height = 12)
  par(mfrow = c(3, 2), mar = c(8, 5, 3, 1), oma = c(2, 2, 2, 1))
  
  # 6. Define colors for intervals
  interval.colors <- terrain.colors(length(intervals))
  
  # 7. Plot each row
  for (i in 1:nrow(counts.sub)) {
    # Prepare data
    row.counts <- unlist(counts.sub[i, ])
    y.max <- max(row.counts) * 1.5
    
    # Create main plot
    bp <- barplot(
      row.counts,
      main       = paste("Row/Class", i),
      ylab       = "Counts",
      col        = interval.colors,
      border     = NA,
      las        = 2,
      names.arg  = custom.abbrev,
      cex.names  = 0.9,
      ylim       = c(0, y.max)
    )
    
    # Add features
    grid(NA, NULL, col = "gray90", lty = 3)
    box()
    
    # Add value labels for top 3 categories in that row
    sorted.counts <- sort(row.counts, decreasing = TRUE)
    top.counts <- head(sorted.counts, n = 3)
    
    for (j in seq_along(row.counts)) {
      if (row.counts[j] %in% top.counts) {
        text(bp[j], row.counts[j],
             labels = format(row.counts[j], big.mark = ","),
             pos = 3, cex = 0.7, col = "darkred")
      }
    }
    
    # Add legend on first plot of each page
    if (i %% 6 == 1) {
      legend("topright",
             legend = custom.abbrev,
             fill   = interval.colors,
             title  = "Weight Intervals (kg)",
             cex    = 0.7,
             ncol   = 2,
             bg     = "white")
    }
  }
  
  dev.off()
  cat("Successfully generated plots at:", save.path, "\n")
}

plot.rowwise.intervals(counts.df, save.path = sprintf("%s/row_wise_plots_rebin_%d.pdf", data.rebin.pwd, year))
