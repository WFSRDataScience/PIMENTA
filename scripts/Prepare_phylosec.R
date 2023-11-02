#---
#title: "Prepare output of Pimenta to analyse with Phyloseq"
#author: "Margot Maathuis"
#Edited by: "Valerie van der Vorst"
#date: '2023-10'
#---

dir.create(path = Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)

list.of.packages=c("stringr", "dplyr", "data.table", "tidyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, lib = Sys.getenv("R_LIBS_USER"), repos='http://cran.us.r-project.org', INSTALL_opts = c('--no-lock'))


##Clear workspace
rm(list = ls())

#Check arguments
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("Please specifiy workdir clustercontent_file marker_ps_fasta_file marker", call.=FALSE)
}

#Inputs
workdir=args[1]
cc_file=args[2]
ps_fasta=args[3]
marker=args[4]

# Load the required libraries
library(data.table)
library(stringr)
library(dplyr)
library(tidyr)

setwd(workdir)
##### Make otu_table #####

# Read the readcountfile using fread() from data.table       FILE CAN BE FOUND IN GENERAL OUTPUT FOLDER OutDir/runname
clustcont <- fread(cc_file)
clustcont <- as.data.frame(clustcont)
#sum(clustcont$Size)                                         #gives numbers of sequences within clusters, but does not tell how many reads are in a consensus
clustcont <- subset(clustcont, DNA_barcode == marker)      #only needed if more primers present
data <- clustcont[,c("ClusterID", "Size", "Clusters_from_samples")]

# Function to split and create new rows
split_and_expand <- function(ClusterID, Clusters_from_samples) 
  {
  split_strings <- unlist(strsplit(Clusters_from_samples, ".combined.fasta:"))
  data.frame(
    ClusterID = rep(ClusterID, length(split_strings)),
    String = split_strings)
  }

# Apply the function to each row and bind the results & CHECK IF IT WORKED!
exp_data <- data %>%
  rowwise() %>%
  do(split_and_expand(.$ClusterID, .$Clusters_from_samples))
table(exp_data$ClusterID)                    #Check by comparing counts to data$Size

# READCOUNT. Extract the number after 'size=' 
Readcount <- str_extract(exp_data$String, "(?<=size=)[0-9]+")

# BARCODE. Extract the strings after the hyphen and starting with 'barcode'
Barcode <- sub(".*-(barcode\\d+).*", "\\1", exp_data$String) 

# SAMPLENAME. Change this for the different runs. Want to check the names? See file: Lab_DNA_Barcodes_SWIMSP_6runs.xlxs
Sample <- sub(".*([r]6\\w+).*", "\\1", exp_data$String)
unique(Sample)        # Check if names of controls were captured correctly

# Make df
exp_data <- cbind(exp_data, Readcount, Barcode, Sample) 
str(exp_data)
exp_data$Readcount <- as.numeric(exp_data$Readcount)
#table(exp_data$Readcount)    # >5, as I specified I want at least 5 reads in a cluster (first clustering step)

# Make OTU table --> Also possible to use Barcodes instead of sample names (just replace Sample for Barcode)
otustart <- aggregate(Readcount ~ ClusterID + Sample , exp_data, FUN = sum)
otucount <- aggregate(Readcount ~ ClusterID + Barcode , exp_data, FUN = length)       #QUESTION! why is it possible to have multiple rows for one barcode + clusterID combination
otu_table <- spread(otustart, Sample, Readcount)
otu_table[is.na(otu_table)] <- 0

otu_table

#PHYLOSEC requirement: table needs sample names and OTU names --> otumatrix
otumat <- otu_table
rownames(otumat) <- paste0("OTU", otu_table$ClusterID)
colnames(otumat)
new_colnames <- gsub("^r", "Sampler", colnames(otumat))
colnames(otumat) <- new_colnames
otumat <- otumat[,-1]       #remove clusterID column after checking if names are rownames are correct
class(otumat)               #should be matrix
otumat <- as.matrix(otumat)
class(otumat)               #Okay? matrix array

otumat

##### Make table with sequence info ##### Not directly needed for phylosec, but want to check if ClusterIDs in PS.fasta are similar to clustercontent (assumption)

# Read the sequence info file using fread() from data.table       MAKE SURE RIGHT PS.fasta FILE IS LOADED! From Barcodes folder
clust <- fread(ps_fasta, header = F)
clust <- as.data.frame(clust)

# CLUSTERID. Extract the number after '>' and remove other empty lines
ID <- str_extract(clust$V1,  "(?<=\\>)\\d+")
keep <- seq(1, nrow(clust), by = 2)
ID <- ID[keep]
unique(ID)         #check: this still works fine

# SEQUENCE. Only keep the sequence information (every 2nd line)
keep <- seq(2, nrow(clust), by = 2)
seqs <- clust[keep,]

# Check length of sequences
hist(nchar(seqs), 300)
table(nchar(seqs))
length(which(nchar(seqs) < 300))
length(which(nchar(seqs) > 650))

# Make df and check
toblast <- data.frame(ClusterID = ID, Sequence = seqs)
length(unique(toblast$ClusterID))
length(unique(toblast$Sequence))   
notunique <- table(toblast$Sequence)
names(notunique[notunique > 1])                              

# Combine into one df IF CLUSTERIDS ARE EQUAL (important assumption, otherwise data will be swapped)
str(toblast$ClusterID)           #character --> Change
str(otu_table$ClusterID)         #integer
toblast$ClusterID <- as.integer(toblast$ClusterID)
length(unique(toblast$ClusterID))      #check if length is equal to x obs (number of rows)
length(unique(otu_table$ClusterID))    #check if length is equal to x obs (number of rows)
toblast$ClusterID
otu_table$ClusterID
identical(sort(toblast$ClusterID), otu_table$ClusterID)     #ONLY CONTINUE IF TRUE!

otu_summary <- merge(otu_table, toblast, by = "ClusterID") 
unique(otu_summary$ClusterID)



#export
readcountcsv_file=paste0(workdir, "/", marker,"_OTU_readcounts.csv", sep="")
readcount_summary_file=paste0(workdir, "/", marker,"_OTU_readcounts_summary.csv", sep="")
write.csv(otu_table, readcountcsv_file, row.names = FALSE)    #export
write.csv(otu_summary, readcount_summary_file, row.names = FALSE) 




