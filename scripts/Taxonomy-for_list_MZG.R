#! /usr/bin/env Rscript

## Valerie van der Vorst 22 october 2022
# outputs list with taxid per sequence ID in stdout, alternative for get_taxids_mzg.py
# requires taxonomizr sql database, can be installed with this script.
# input is a txt file with GenBank ids

# How to run:
# Rscript Taxonomy-for_list_MZG.R accessionTaxa.sql listGenbank 


dir.create(path = Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)

list.of.packages=c("stringr", "dplyr", "readr", "data.table", "taxonomizr", "stringi")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
#old_packages <- list.of.packages[(list.of.packages %in% old.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, lib = Sys.getenv("R_LIBS_USER"), repos='http://cran.us.r-project.org', INSTALL_opts = c('--no-lock'))
#if(length(old_packages)) install.packages(old_packages, lib = Sys.getenv("R_LIBS_USER"), repos='http://cran.us.r-project.org', INSTALL_opts = c('--no-lock'))


##Clear workspace
rm(list = ls())


#Setting up database for taxonomizr
#prepareDatabase('accessionTaxa.sql')
#loading taxdmp
#taxaNodes<-read.nodes('nodes.dmp')
#taxaNames<-read.names('names.dmp')


#Requirements:
library('dplyr')
library('data.table')
library('readr')
library("stringi")
library("taxonomizr")
library("rlist")
args = commandArgs(trailingOnly=TRUE)

#check if taxonomizr database exists
accessionTaxa=args[1]
accessionList=args[2]
if (!file.exists(accessionTaxa)){
print("preparing taxonomizr database")
accessionTaxa="accessionTaxa.sql"
prepareDatabase(accessionTaxa)
}

 
sallac=read.table(accessionList, header=FALSE)
head(sallac)
List = list()
for(y in sallac$V1)
	{
	        TAX <- accessionToTaxa(y,accessionTaxa, version='base')
        	List[[length(List)+1]] = TAX
			#print(TAX)
		cat(y,"\t",TAX,"\n")
	}
t <- unlist(List)

print(length(t))
