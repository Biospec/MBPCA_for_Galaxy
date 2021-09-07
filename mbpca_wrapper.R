#! Rscript

##########################################
# Multiblock PCA wrapper script for Galaxy
##########################################
#Startup log
sink("startup_log.txt")

pkgs=c("RSpectra", "batch", "MASS")

for(pkg in pkgs) {
  suppressPackageStartupMessages( stopifnot( library(pkg, quietly=TRUE, logical.return=TRUE, character.only=TRUE)))
  cat(pkg,"\t",as.character(packageVersion(pkg)),"\n",sep="")
}

listArguments <- parseCommandArgs(evaluate = FALSE) #interpretation of arguments given in command line as an R list of objects
sink()

#Redirect all stdout to the log file
sink(listArguments$information)



# ----- PACKAGE -----
cat("\tPACKAGE INFO\n")
sessionInfo()

source_local <- function(fname) {
  argv <- commandArgs(trailingOnly = FALSE)
  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  source(paste(base_dir, fname, sep="/"))
}

# Load functions
source_local("multiblock.R")
print("Initial loading successful")



cat('\n\nRunning multiblock PCA\n');
options(warn=-1);
#remove rgl warning
options(rgl.useNULL = TRUE);
model_type <- "multiblock PCA" ## module name

cat("\nStart of the '", model_type, "' Galaxy module call: ",
    format(Sys.time(), "%a %d %b %Y %X"), "\n", sep="")
if (listArguments[["instrument"]] == "no" | listArguments[["instMetadata_in"]] == "NA") instrument = FALSE else  instrument = TRUE
# If it is data source blocking, obtain data source info from a different meta file
if (instrument){
  meta <- list()
  meta[[1]] <- read.csv(listArguments[["instMetadata_in"]], header = TRUE)
  tmp <- read.csv(listArguments[["sampleMetadata_in"]])
  if (dim(tmp)[2] > 1) {
    cat("\n Meta data contains more than two columns, only first column used for blocking. \n")
    meta[[2]] <- tmp[, 1]
  } else{
    meta[[2]] <- tmp
  }
  data_mat <- read.csv(listArguments[["dataMatrix_in"]], 
                       header = TRUE,
                       row.names = NULL)
  xaxis <- colnames(data_mat)
  data_mat <- as.matrix(data_mat)} else {
  meta <- read.csv(listArguments[["sampleMetadata_in"]], header = TRUE)
  data_mat <- read.csv(listArguments[["dataMatrix_in"]], 
                       header = TRUE,
                       row.names = NULL)
  col_id <- colnames(data_mat)
  if (col_id[1] == "xaxis"){
    xaxis <- data_mat[, 1]
    data_mat <- t(as.matrix(data_mat[, -1]))
    col_id <- col_id[-1]
  } else {
    xaxis <- NULL
    data_mat <- t(as.matrix(data_mat))
  }
}






if (!is.null(xaxis)) {
  var_id <- xaxis
} else {
  var_id <- 1:dim(data_mat)[2]
}

no_vars <- length(var_id)
no_pcs <- as.numeric(listArguments[["no_pcs"]])

model_collection <- mbpca_host(data_mat, meta, instrument, listArguments[["model_name"]], no_pcs)


##saving
filename_superscores <- listArguments[["file_super"]]
filename_blockscores <- listArguments[["file_block"]]
filename_blockloadings <- listArguments[["file_loadings"]]
filename_figures <- listArguments[["file_figures"]]

if (exists("model_collection")) {
  ## writing output files
  if (instrument){
    cat("\n\nWriting output data files\n\n")
    resultNow <- model_collection
    superscores <- data.frame(cbind(as.character(resultNow$block_label), resultNow$models$super_scores))
    colnames(superscores) <- c("Sample names", "PC 1 super scores", "PC 2 super scores")
    write.table(superscores,
                file = filename_superscores,
                quote = FALSE,
                row.names = FALSE,
                sep = ",")
    blockloadings <- data.frame(cbind(as.character(var_id), resultNow$models$block_loadings))
    colnames(blockloadings) <- c("Variable names", "PC 1 loadings", "PC 2 loadings")
    write.table(blockloadings,
                file = filename_blockloadings,
                quote = FALSE,
                row.names = FALSE,
                sep = ",")
    blockscores <- data.frame(cbind(as.character(resultNow$block_label), resultNow$models$block_scores))
    blockscores_colnames <- list()
    blockscores_colnames[1] <- "Sample names"
    no_blks <- length(resultNow$block_id)
    for (i in 1:no_pcs){
      for (ii in 1:no_blks)
        blockscores_colnames <- c(blockscores_colnames, 
                                  paste(resultNow$block_id[[ii]], ": PC", i, " block scores", sep = ""))
    }
    colnames(blockscores) <- blockscores_colnames
    write.table(blockscores,
                file = filename_blockscores,
                quote = FALSE,
                row.names = FALSE,
                sep = ",")
    
    # Graphical display for each significant parameter
    cat("\n\nWriting figures\n\n")
    pdf(filename_figures, onefile = TRUE)
    plot.mbpca(resultNow$models, resultNow$model_name, resultNow$block_label, resultNow$block_id)
  } else{
    no_models <- length(model_collection)
    for (i in 1:no_models){
      cat("\n\nWriting output data files\n\n");
      resultNow <- model_collection[[i]]
      superscores <- data.frame(cbind(as.character(resultNow$block_label), resultNow$models$super_scores))
      colnames(superscores) <- c("Sample names", "PC 1 super scores", "PC 2 super scores")
      if (i == 1){
        write.table(superscores,
                    file = filename_superscores,
                    quote = FALSE,
                    row.names = FALSE,
                    sep = ",")
        no_blks <- length(resultNow$block_id)
        blockloadings <- data.frame(as.character(rep(var_id, no_blks)), resultNow$models$block_loadings)
        colnames(blockloadings) <- c("Variable names", "PC 1 loadings", "PC 2 loadings")
        write.table(blockloadings,
                    file = filename_blockloadings,
                    quote = FALSE,
                    row.names = FALSE,
                    sep = ",")
        blockscores <- data.frame(cbind(as.character(resultNow$block_label), resultNow$models$block_scores))
        blockscores_colnames <- list()
        blockscores_colnames[1] <- "Sample names"
        
        for (ii in 1:no_pcs){
          for (iii in 1:no_blks)
            blockscores_colnames <- c(blockscores_colnames, 
                                      paste(resultNow$block_id[[iii]], ": PC", ii, " block scores", sep = ""))
        }
        colnames(blockscores) <- blockscores_colnames
        write.table(blockscores,
                    file = filename_blockscores,
                    quote = FALSE,
                    row.names = FALSE,
                    sep = ",")
    
        cat("\n\nWriting figures\n\n");
        pdf(filename_figures, onefile = TRUE)
        plot.mbpca(resultNow$models, resultNow$model_name, resultNow$block_label, resultNow$block_id)
      } else {
        write.table(superscores,
                    file = filename_superscores,
                    append = TRUE,
                    quote = FALSE,
                    row.names = FALSE,
                    sep = ",")
        no_blks <- length(resultNow$block_id)
        blockloadings <- data.frame(cbind(as.character(rep(var_id, no_blks)), resultNow$models$block_loadings))
        colnames(blockloadings) <- c("Variable names", "PC 1 loadings", "PC 2 loadings")
        write.table(blockloadings,
                    file = filename_blockloadings,
                    append = TRUE,
                    quote = FALSE,
                    row.names = FALSE,
                    sep = ",")
        blockscores <- data.frame(cbind(as.character(resultNow$block_label), resultNow$models$block_scores))
        blockscores_colnames <- list()
        blockscores_colnames[1] <- "Sample names"
        
        for (ii in 1:no_pcs){
          for (iii in 1:no_blks)
            blockscores_colnames <- c(blockscores_colnames, 
                                      paste(resultNow$block_id[[iii]], ": PC", ii, " block scores", sep = ""))
        }
        colnames(blockscores) <- blockscores_colnames
        write.table(blockscores,
                    file = filename_blockscores,
                    append = TRUE,
                    quote = FALSE,
                    row.names = FALSE,
                    sep = ",")
        cat("\n\nWriting figures\n\n")
        plot.mbpca(resultNow$models, resultNow$model_name, resultNow$block_label, resultNow$block_id)
      }
       
     } 
  }
  dev.off()
}

tryCatch({
  save.image(file="mbpca.RData");
}, warning = function(w) {
  print(paste("Warning: ", w));
}, error = function(err) {
  stop(paste("ERROR saving result RData object:", err));
});

## ending
##-------

cat("\nEnd of the '", model_type, "' Galaxy module call: ",
    format(Sys.time(), "%a %d %b %Y %X"), "\n", sep = "")

sink()

rm(list = ls())
