#!#!/usr/bin/env Rscript

plot_coverage<-function(input_folder, nchr){
if (!dir.exists(input_folder)){
  stop("Error: The input folder does not exist. Please check the path.")
}
fileslist<-list.files(input_folder, pattern="\\.coverage$", full.names=TRUE)
if (length(fileslist)==0){
  stop("Error: No coverage files found in the specified input folder.Please ensure the folder contains the expected files with extention <.coverage> .")
}

pdf(file.path(input_folder,"ALLcov.pdf"))
par(mfrow = c(2, 2))
for (eachfile in fileslist){
  SAMP_COV<-read.table(eachfile)
  chromosomes <- paste("CHR", 1:nchr, sep="")
  chromosomes <- c(chromosomes, "MEAN")
  barplot(c(SAMP_COV[c(1:nchr),7],mean(SAMP_COV[c(1:nchr),7])), 
          names.arg=chromosomes, las=2, main=basename(eachfile))
  axis(4)
  }
dev.off()

pdf(file.path(input_folder, "ALLcovP.pdf"))
par(mfrow = c(2, 2))
for (eachfile in fileslist){
  SAMP_COV<-read.table(eachfile)
  chromosomes <- paste("CHR", 1:nchr, sep="")
  chromosomes <- c(chromosomes, "MEAN")
  barplot(c(SAMP_COV[c(1:nchr),6],mean(SAMP_COV[c(1:nchr),6])), 
          names.arg=chromosomes, las=2, main=basename(eachfile))
  axis(4)
}
dev.off()
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Please provide at least an input folder as an argument, and optionally, the number of chromosomes.\nUsage: plot_coverage.R <input folder> <number of chromosomes>\nExample: plot_coverage.R /Yuma/coverage/ 24\n\n")
}

Rinput_folder <- args[1]
n_chr <- if (length(args) > 1) args[2] else 24
plot_coverage(Rinput_folder, n_chr)
