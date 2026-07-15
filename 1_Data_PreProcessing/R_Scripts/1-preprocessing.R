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
anyNA(g_data) #no null values found
sum(is.na(g_data))
# colSums(is.na(g_data)) # check for each columns

#load the ExpressionSet
#sp: I added the a directory path for downloading the data {getGEO()}. otherwise. it is downloaded to temp derectory. i have to download it
#sp: every time i re-open r
#sp: reference for getGEO : https://www.rdocumentation.org/packages/GEOquery/versions/2.38.4/topics/getGEO

gse_150910 <- getGEO("GSE150910",GSEMatrix  = TRUE,AnnotGPL = FALSE, destdir="../../DATASET")
# Extract metadata
metadata <- pData(gse_150910[[1]])
# to see how many data we have in header.
dim(metadata)

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

# meta_clean <- meta_clean[colnames(g_data), ]

# aligning
meta_clean <- meta_clean[colnames(g_data), ]
#check the alignment is correct in both datasets
all(colnames(g_data) == rownames(meta_clean))

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


# befor plotting the data. we need to know the minimum and maximum values of the data
print(min(gene_counts_df$count))
print(max(gene_counts_df$count))

# lets try to plot using Base R
# i want to divide the range to 100 bins.
# that means bin width is  = {max-min}/100
# it is arround = 2775142.76 
# so a bar represent a huge range
# and in this dataset most of the data is at the begining
# so it is hard to visulize the whole data using this.
# but lets try

# here i want to see the pobabilty distribution. not the frequency.
# that why probability = TRUE,
# breaks = 100, for dividing the range
hist(gene_counts_df$count, 
     probability = TRUE, 
     breaks = 100, 
     col = "lightblue", 
     main = "Histogram with Normal Curve", 
     xlab = "Gene Counts")

# now i want to plot a curve to see how it would be if this is a normal distribution
# we can do that using 'cruve(dnorm())' 
# we need to enter the mean and standard deviation of our data as parameters
mean_val <- mean(gene_counts_df$count)
sd_val <- sd(gene_counts_df$count)

curve(dnorm(x, mean = mean_val, sd = sd_val), 
      add = TRUE, 
      col = "red", 
      lwd = 2)


# now we see that this is not a good representation.
# lets rescale the x axis to logorithms

# to do that we can use 'ggplot' 's scale_x_log10 attribute
# refere this to understand how it works : https://www.r-bloggers.com/2021/08/beginning-a-ggplot2-series-logarithmize-your-scales/
# please watch the video and understand the easthetic function : https://www.youtube.com/watch?v=CLUnPk_1oaE&t=149s
# in aesthetic we have add 1 for every gene counts. we call it 'Pseudo-count' 
# we are doing it because there is no definition for log0. but we have 0 values

install.packages("ggplot2")
library(ggplot2)

ggplot(gene_counts_df, aes(x = count + 1)) +
  geom_histogram(bins = 100, fill = "steelblue", color = "black") +
  scale_x_log10() +
  labs(title = "Log10 Transformed Gene Counts",
       x = "Gene Counts (Log10 scale, count + 1)",
       y = "Frequency (Number of Genes)") +
  theme_minimal() # to have a classic backgroud for the graph

# now we can see thar in some specific ranges in left tail gene counts are very low. there are many genes near zero expression counts


# we can cleary see that using a 'x intercept'
# lets make a density curve for better smoothness
# we cut the genes that has expression for all samples under 10. for 288 samples. those are noise

ggplot(gene_counts_df, aes(x = log10(count + 1))) +
  geom_density(fill = "purple", alpha = 0.5, color = "black") +
  geom_vline(xintercept = log10(80 + 1), color = "red", linetype = "dashed", size = 1) +
  labs(title = "Density Plot of Gene Counts",
       x = "Log10(Count + 1)",
       y = "Density") +
  theme_minimal()


###############################################################################################################
# LOW EXPRESSION DATA
###############################################################################################################

# we want to remove the genes that are expression level under the 80 , from our study.
# why 80 ???
# *** samaples are chp=82, ipf=103, control=103.
# so that ,to show a class realted feature, the total expression count for a perticular gene for all subjects should be aleast 82
# so that we cutoff the noise off from 80

# we take the count of 'gene_count*' dataframe's 'count' under 80
# here we get and idea about how much we can remove
low_expr_count <- sum(gene_counts_df$count < 80)
print(paste("Low expression count : ", low_expr_count))

# then we take those gene_names as a list so that we can remove them from our 'gene_count_df'
low_expr_gene_names <- gene_counts_df$gene_name[gene_counts_df$count < 80]

# just see what we are removing.
head(low_expr_gene_names)

# we can make the filetered dataset by taking the all the rows that has not 'gene_name' in 'low_expr_gene_names', with all the columns

g_data_rmvd_low <- g_data[!(rownames(g_data) %in% low_expr_gene_names), ]

# lets check the dimensions of the data so that we can verify whether the cutting correct or wrong
dim(g_data_rmvd_low) # it has 17931 entries we cutted 907 entrie . so total is 18838. data frame is correct

# we can use write.csv to save the dataframe as a csv file. if we want enter the data row names , we need to specify 'row.names' attribute as true
write.csv(g_data_rmvd_low, file = "../../DATASET/g_data_rmvd_low.csv", row.names = TRUE)


#** LETS PLOT THE DENSITY CURVE NOW **#

# as the way did before.
# take the count across rows
# get a list of summatinons
# make a dataframe only contain gene_name and summations
# plot using 'ggplot'

gene_counts_rmvd_low = apply(g_data_rmvd_low, MARGIN= 1 , FUN=sum)
print(gene_counts_rmvd_low)

gene_counts_rmvd_low_df = data.frame(
  gene_name = names(gene_counts_rmvd_low),
  count = as.numeric(gene_counts_rmvd_low)
)

ggplot(gene_counts_rmvd_low_df, aes(x = log10(count + 1))) +
  geom_density(fill = "yellow", alpha = 0.5, color = "black") +
  geom_vline(xintercept = log10(80 + 1), color = "red", linetype = "dashed", size = 1) +
  labs(title = "Density Plot of Gene Counts",
       x = "Log10(Count + 1)",
       y = "Density") +
  theme_minimal()

# now we can see that the low noises has been removed


###############################################################################################################
# DESeq2 Median of Ratios
###############################################################################################################

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

library(DESeq2)

# Before we want to continue with DESeq2 we need to make a DESeqDataSet object.
# To make that object we need to enter the count data and coldata. for that we use 'g_data_rmvd_low' and 'meta_clean'
# If we are using them to combined object, those columns of both dataframes should match each other.
# so that we need to do some checkings

# check the sample names are available in both dataframes
print(rownames(g_data_rmvd_low))
print(rownames(meta_clean))

length(colnames(g_data_rmvd_low))
length(rownames(meta_clean))
# we can see that all the names are available. 
                    
# Now that we need to align the both dataframes by their names.
# We can reassign the head_clean in same order with the g_data
meta_clean <- meta_clean[colnames(g_data), ]

# lets check is it aligned
head(colnames(g_data_rmvd_low), 10)
head(rownames(meta_clean), 10)
# ok it is aligned

# lets make the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(
  countData = g_data_rmvd_low,   # our filtered raw count matrix (genes x samples)
  colData   = meta_clean,        # our sample metadata (must match column order of countData)
  design    = ~ diagnosis        # we will eventually compare CHP vs IPF vs Control
)
dds
#this gives a warning because of the variables in the design formula represent experimental groups. DESeq2 performs differential expression analysis by comparing these groups, so they must be categorical variables (factors) rather than plain text (characters).
class(dds)

meta_clean$diagnosis <- factor(meta_clean$diagnosis)#converts the diagnosis column in meta_clean from a character variable to a factor (categorical variable)
class(meta_clean$diagnosis) #checks the data type (class) of the diagnosis column in the meta_clean
#recreate the object
dds <- DESeqDataSetFromMatrix(
  countData = g_data_rmvd_low,
  colData   = meta_clean,
  design    = ~ diagnosis
)

#now also give a warning but it can be ignore becuase it says the chp and control is not in correct order.
dds #now we can obtain correct data 

levels(colData(dds)$diagnosis)


# 2026-07-15 start in here

###############################################################################################################
# NORMALIZATION & VST TRANSFORMATION
###############################################################################################################

#calculates a normalization factor for each sample
dds <- estimateSizeFactors(dds)
sizeFactors(dds)

#Apply Variance Stabilizing Transformation (VST)
#because of the low expression and high expression noises we do that. this is necessary for the PCA and heatmap 
vsd <- vst(dds, blind = TRUE) # blind =true means this doing with not considering the sample variation.
vsd
vst_matrix <- assay(vsd) #Only the transformed expression values(rpw= genes, columns=samples)
cat("VST matrix dimensions:", nrow(vst_matrix), "genes x", ncol(vst_matrix), "samples\n")


###############################################################################################################
# VISUALIZATION - PCA PLOT
###############################################################################################################


pca_result <- prcomp(t(vst_matrix), scale. = FALSE)
pca_var <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)

pca_df <- data.frame(
  PC1       = pca_result$x[, 1],
  PC2       = pca_result$x[, 2],
  diagnosis = meta_clean$diagnosis,
  batch     = meta_clean$batch
)

pca_plot <- ggplot(pca_df, aes(PC1, PC2, color = diagnosis)) +
  geom_point(size = 2.5, alpha = 0.8) +
  labs(
    title = "PCA - GSE150910 (Gene-Level Counts)",
    x     = paste0("PC1 (", pca_var[1], "% variance)"),
    y     = paste0("PC2 (", pca_var[2], "% variance)")
  ) +
  scale_color_manual(values = c("chp" = "#E74C3C",
                                "control" = "#2ECC71",
                                "ipf" = "#3498DB")) +
  theme_bw() +
  theme(legend.title = element_text(face = "bold"))
print(pca_plot)

ggsave("../../Plots/PCA_plot.png",
       plot = pca_plot, width = 8, height = 6, dpi = 300)
