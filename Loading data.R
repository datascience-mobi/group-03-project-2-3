#LOADING DATA

#Download necessary package
library(rstudioapi)

#Set wd to the working directory where this Rscript is located
wd = dirname(rstudioapi::getSourceEditorContext()$path)

#Load data
basalexpression <- readRDS(paste0(wd,"/CCLE_basalexpression.rds"))
copynumber <- readRDS(paste0(wd,"/CCLE_copynumber.rds"))
mutations <- readRDS(paste0(wd,"/CCLE_mutations.rds"))
treated <- readRDS(paste0(wd,"/NCI_TPW_gep_treated.rds"))
untreated <- readRDS(paste0(wd,"/data/NCI_TPW_gep_untreated.rds"))
drugsensitivity <- readRDS(paste0(wd,"/data/NegLogGI50.rds"))

cannotations <- read.table(paste0(wd,"/cellline_annotation.tsv"), header = TRUE, sep ="\t", stringsAsFactors = TRUE)
dannotations <- read.table(paste0(wd,"/drug_annotations.tsv"), header = TRUE, sep ="\t", stringsAsFactors = TRUE)
metadata <- read.table(paste0(wd,"/data/NCI_TPW_metadata.tsv"), header = TRUE, sep ="\t", stringsAsFactors = TRUE)
