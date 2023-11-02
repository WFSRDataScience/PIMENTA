#!/usr/bin/env Rscript

##Script for making taxonomic summaries across DNA barcodes
##M.Staats, 4 February 2020

#optional: clear workspace
rm(list = ls())

#check agruments
args <- commandArgs(trailingOnly=TRUE)

if(length(args)<1)
{
	stop("Please specify input directory: <PATH to InFile>")
} else if (length(args)==1) {
	# default output file
	args[4] = "out.txt"
}

InDir <- args[1]
SampleID <- basename(InDir)
print(SampleID)

dir.create(path = Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)

list.of.packages=c("gplots", "dplyr", "data.table", "stringr", "purrr", "readr", "Hmisc", "gdata")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, lib = Sys.getenv("R_LIBS_USER"), repos='http://cran.us.r-project.org', INSTALL_opts = c('--no-lock'))


library(gplots,quietly = T)
#library(plyr)
library(dplyr)
library(data.table, quietly = T)
library(stringr, quietly = T)
library(readr, quietly = T)
library(Hmisc, quietly = T)
library(gdata, quietly = T)
library("purrr", quietly = T)
###Find taxonomy files
file.names <- dir(InDir, pattern = "*.tsv", full.names = TRUE, ignore.case = TRUE, recursive = TRUE)
print(file.names)
print("length")
for(i in 1:length(file.names)){
  file <- read.table(file.names[i], header=FALSE, sep="\t", stringsAsFactors=FALSE, fill=TRUE)
#old: #qseqid TAXID   Taxonomic Level Taxonomic Name  Number of Top Hits      Lineage Reads   Length  Marker
#new: #DNA_barcode	Taxonomiclevel    Taxonomy        Lineage Readcount
##Extract Barcode Name
  NAME <- basename(file.names[i])
  DirTable <- dirname(file.names[i]) 
  BARCODE <- gsub(pattern = "\\.tsv$", "", NAME)
  BARCODE2 <- gsub(pattern = paste(SampleID, ".", sep=""), "", BARCODE)
#  BARCODE2 <- paste0("BARCODE_", BARCODE1)
  print(BARCODE2)  
  ##Extract read number from OTU-name
  out <- file$qseqid
  out2 <- cbind(file[c(2:6)])
  head(out2)  
  #Rename column names
  colnames(out2) <- c("Taxid","Taxonomic.Level", "Taxonomic.Name" ,"Lineage","Total.Reads")
  out2$Lineage <- gsub(",",":",as.character(out2$Lineage))
  print(head(out2))
  #rename family_low_hits to family
  out2$Taxonomic.Level <- gsub('_low_hits', '', out2$Taxonomic.Level)
  ##Calculate Total Reads
  ReadsSum <- colSums(out2[ , 5, drop = FALSE])
  ##Remove rows with unidentified taxa <NA>
  out3 <- out2[complete.cases(out2), ]
  ##Calcule Total Reads without unidentifed taxa
  ReadsSum2 <- colSums(out3[ , 5, drop = FALSE])
  ##Add row for taxa without taxonomic ID
  UnidNo <- ReadsSum - ReadsSum2
  out3 <- rbind(out3, c("Unidentified", "Unidentified", "Unidentified", "Unidentified", UnidNo))
  ##Group by Taxonomic.Name to remove redundant entries, while adding to total sum
  out3$Total.Reads <- as.numeric(as.character(out3$Total.Reads))
  out4 <- aggregate(out3$Total.Reads, by=list(Taxonomic.Name=out3$Taxonomic.Name, Taxonomic.Level=out3$Taxonomic.Level, Lineage=out3$Lineage), FUN=sum)
  ##Set barcode name to Total.Reads column
  colnames(out4)[4] <- BARCODE2
  write.csv(out4, file=paste0(DirTable, "/", BARCODE, ".Taxonomic.summary.csv"),row.names=FALSE, quote=FALSE)
}

##Create summary across Barcodes
data_files <- dir(InDir, pattern = "*Taxonomic.summary.csv", full.names = TRUE, ignore.case = TRUE, recursive = TRUE)
print(data_files)
list_of_data_sets <- map(data_files, read.csv)
#list_of_data_sets <- read_list(data_files, read_csv)
multi_join <- function(list_of_loaded_data, join_func, ...){
    output <- Reduce(function(x, y) {join_func(x, y, ...)}, list_of_loaded_data)
    return(output)
}
merged_data <- multi_join(list_of_data_sets, full_join)
##sort Table based on Taxonomic.Level
merged_data$Taxonomic.Level <- as.factor(merged_data$Taxonomic.Level)
merged_data$Taxonomic.Level <- reorder(merged_data$Taxonomic.Level, new.order=c("species", "genus", "family", "order", "class", "phylum", "kingdom", "superkingdom", "Unidentified"))
merged_data2 <- merged_data[order(merged_data$Taxonomic.Level),]
write.table(merged_data2, file=paste0(InDir, "/", SampleID, ".Taxonomic.summary.all_barcodes.tsv"),row.names=FALSE, quote=FALSE, sep="\t")
