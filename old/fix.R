rm(list=ls())
library(rpart)
library(rpart.plot)

# Load your data and tree
year <- 2021
data.pwd <- sprintf("/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/birthweight_data/rebin/")
load(sprintf("%s/dm_tree_%d.RData", data.pwd, year))
load(sprintf("%s/informed_prior_%d.RData", data.pwd, 2020))
plots.pwd <- sprintf("/Users/adamkurth/Documents/RStudio/conformal-lbw-prediction/plots")


visualize_dm_tree <- function(tree, file_path=NULL) {
  if(!is.null(file_path)) {
    pdf(file_path, width = 12, height = 8)
    on.exit(dev.off())
  }
  
  par(xpd=TRUE)  # Allow plotting outside the figure region
  plot(tree, 
       uniform=TRUE, 
       # branch=0,
       # compress=FALSE, 
       # margin=0.15, 
       main="Dirichlet-Multinomial Decision Tree"
 )
  
  # Add node text with better formatting
  text(tree, use.n=FALSE, pretty=TRUE, all=FALSE, cex=1, xpd=TRUE, font=3)

  # Return success message
  if(!is.null(file_path)) {
    return(paste("Tree visualization saved to", file_path))
  } else {
    return("Tree visualization complete")
  }
}

# Test the visualization function in console
message(visualize_dm_tree(dm.tree))

# Save to PDF file
pdf_path <- file.path(plots.pwd, paste0("dm_tree", year, ".pdf"))
message(visualize_dm_tree(dm.tree, pdf_path))