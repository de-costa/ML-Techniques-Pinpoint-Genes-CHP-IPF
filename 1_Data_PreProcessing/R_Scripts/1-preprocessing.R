# i want to remove the low expressed genes first. But we cant remove a arbitrary value from subjects.
# if we remove we have to remove entire rows.
# so that we can take total gene count in each row. 
# then we can plot them in a histogram or any better way to see outliers
# so we move those outliers

# But there is one problem. classes are not balanced. low expression of a gene may be due that imbalance.
# what if chp is the only thing that expressed the gene. then we will get a low expression too.


# Bioconductor provide software tool for bioinformatics analysis. It suggest that for they have a package called
# GEOquerry that can be used to work with GSE data. 
# for that we need to install 'BiocManager' and install the 'GEOquery' using their package manager.

# this {install.packages("BiocManager")} is anoying. so wrap it out . it only install when it is not already installed

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("GEOquery", "DESeq2", "limma"),
                     update = FALSE)

install.packages(c("dplyr", "ggplot2", "pheatmap"),
                 dependencies = TRUE)

library(GEOquery)
library(DESeq2)
library(limma)
library(dplyr)
library(ggplot2)
library(pheatmap)
#sp: I added the a directory path for downloading the data {getGEO()}. otherwise. it is downloaded to temp derectory. i have to download it
#sp: every time i re-open r
#sp: reference for getGEO : https://www.rdocumentation.org/packages/GEOquery/versions/2.38.4/topics/getGEO




# 
# wanna see the what data is present in a sample. transpose it so that easier to read
t(metadata[1,])

# find how many samples are there in each category
print(table(metadata$`diagnosis:ch1`))


all(colnames(g_data) %in% meta_clean$sample_id)
# Create clean metadata table with relevant columns
meta_clean <- data.frame(
  sample_id = metadata$title,
  diagnosis = metadata$`diagnosis:ch1`,
  batch     = metadata$`batch:ch1`,
  sex       = metadata$`Sex:ch1`,
  age       = metadata$`age:ch1`,
  row.names = metadata$title
)

head(meta_clean)

# table(meta_clean$diagnosis)

table(meta_clean$batch) # according to the data, if we use the raw data. we will face batch effect. Better to use preprocessed count

###@@@@@@@@issue@@@@@@@###  count kiyanne mokadda

cat("All samples found in metadata:", all(colnames(counts) %in% rownames(meta_clean)), "\n") 

# Reorder metadata to match count matrix column order
meta_clean <- meta_clean[colnames(counts), ]

# Verify perfect alignment
cat("Perfect alignment:",
    all(colnames(counts) == rownames(meta_clean)), "\n")

##metadata is loaded correctly and check whether it is matching with g_data
##----------------------------------------------------->


###############################################################################################################
# Start the process of genotype data and metadata
###############################################################################################################

# lets import our dataset
g_data <- read.csv("../../DATASET/GSE150910_gene-level_count_file.csv", header = TRUE, sep = ",", row.names = 1) 
# we have a comma seperrated file

head(g_data)
# columns = subjects
cat("  Genes:   ", nrow(g_data), "\n") # rows = genes (18838)
cat("  Samples: ", ncol(g_data), "\n") # column = samples (288)

# chp = 82
# ipf = 103
# control = 103 

# lets see null values
anyNA(g_data)
sum(is.na(g_data))
# colSums(is.na(g_data)) # check for each columns


gse_150910 <- getGEO("GSE150910",GSEMatrix  = TRUE,AnnotGPL = FALSE, destdir="../../DATASET")
# Extract metadata
metadata <- pData(gse_150910[[1]])
# to see how many data we have in header.
dim(metadata)
cat("\nDiagnosis groups:\n")
print(table(metadata$`diagnosis:ch1`))


###############################################################################################################
# testing the data of the header and counts are matches
###############################################################################################################
meta_clean <- data.frame(
  sample_id = metadata$title,
  diagnosis = metadata$`diagnosis:ch1`,
  batch     = metadata$`batch:ch1`,
  sex       = metadata$`Sex:ch1`,
  age       = metadata$`age:ch1`,
  row.names = metadata$title
)

print(head(meta_clean))

print(table(meta_clean$diagnosis))

meta_clean <- meta_clean[colnames(g_data), ]

#check the alignment is correct in both datasets
all(colnames(counts) == rownames(meta_clean))

###############################################################################################################
# PLOT noises
###############################################################################################################

# i want to add the every count of a gene of all samples. That means i want to get the summation of values along rows.
# then from all rows we can plot histogram
# to see the out liers.

# i am using 'apply()' function for this and default 'plot()'
# refere the links to understand
# https://www.geeksforgeeks.org/r-language/apply-lapply-sapply-and-tapply-in-r/
# https://youtu.be/_V8eKsto3Ug 44:44

gene_counts = apply(g_data, MARGIN= 1 , FUN=sum)
print(gene_counts)

# this is not good way to watch them. it is messy

type(gene_counts) # int - means this is a pairlist of integers and names. we can seperate them and make a datafreame

# refere 
# for integer separation : https://www.geeksforgeeks.org/r-language/as-numeric-function-in-r/
# for name seperation : https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/names

# lets make a dataframe with coulumns 'gene_name' and 'count'
gene_counts_df = data.frame(
  gene_name = names(gene_counts),
  count = as.numeric(gene_counts)
)

print(gene_counts_df) # better


###############################################################################################################
# LOW EXPRESSION DATA
###############################################################################################################
