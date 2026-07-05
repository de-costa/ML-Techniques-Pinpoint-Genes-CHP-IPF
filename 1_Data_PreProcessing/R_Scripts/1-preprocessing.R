# lets import our dataset
g_data <- read.csv("../../DATASET/GSE150910_gene-level_count_file.csv", header = TRUE, sep = ",") # we have a comma seperrated file

head(g_data)
# columns = subjects
  # chp = 82
  # ipf = 103
  # control = 103 
# rows = genes (18838)
# total_data = (18838)*(82+103+103) = 5425344

# lets see null values
anyNA(g_data)
sum(is.na(g_data))
# colSums(is.na(g_data)) # check for each columns

# i want to remove the low expressed genes first. But we cant remove a arbitrary value from subjects.
# if we remove we have to remove entire rows.
# so that we can take total gene count in each row. 
# then we can plot them in a histogram or any better way to see outliers
# so we move those outliers

# But there is one problem. classes are not balanced. low expression of a gene may be due that imbalance.
# what if chp is the only thing that expressed the gene. then we will get a low expression too.
sapply(g_data, sum, 1)
