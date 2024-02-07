# Load necessary libraries
library(DESeq2)

# Get input and output file paths from command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]
#input_file <- "/home/alireza/02-03/Bio/project/BIO_Project/tasks/the-songbird-study/dataset/normalized_dataset.csv"
#output_file <- "/home/alireza/02-03/Bio/project/BIO_Project/tasks/the-songbird-study/result_2.csv"


# Read the input CSV file
data <- read.csv(input_file)

# Find columns with NaN values
nan_columns <- colnames(data)[colSums(is.na(data)) > 0]
# Remove columns with NaN values
data <- data[, !(colnames(data) %in% nan_columns)]

# Convert 'social_setting' and 'study_group' and 'tissue_id' to factors
data$social_settting <- factor(data$social_settting)
data$study_group <- factor(data$study_group)

# Extract metadata columns as a dataframe
metadata_df <- data.frame(study_group = data$study_group, 
                       social_settting = data$social_settting)

a <- data$sample_id
rownames(metadata_df) <- a 

# Remove metadata columns from data
data <- data[, !(names(data) %in% c("sample_id", "study_group", "social_settting", "tissue_id"))]

# Transpose data matrix
data <- as.matrix(data)
data <- t(data)
data <- matrix(as.integer(data), nrow = nrow(data), ncol = ncol(data))
#colnames(data) <- a
epsilon <- matrix(1, nrow = nrow(data), ncol = ncol(data))
data <- data + epsilon
# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = metadata_df,
                              design = ~ study_group + social_settting)

# Run DESeq2 analysis
dds <- DESeq(dds, fitType = 'mean')

# Get normalized counts

# Write results to CSV
dds <- results(dds)
write.csv(dds, file = output_file, row.names = TRUE)