#Title: "SCA1 Excercise vs Control Mice"
#Authors: Sharan Srinivasan, Mariana Sierra
#date: "05/20/2024"

# script to visualize gene expression data
setwd("/Users/marianasierra/Desktop/srinivasan_bulk_rna_seq")

# load libraries
library(tidyverse)
library(ggplot2)

# data Upload

# Read in the file. Don't set row names yet
# Note if using R < 4.0.0, set stringsAsFactors = FALSE in read.delim
data_1 <- read.delim("/Users/marianasierra/Desktop/srinivasan_bulk_rna_seq/data/gene_expected_count.annot.txt", row.names = NULL)
# Deal with genes that don't have annotated gene symbols (external_gene_name)
# Use ENSEMBL ID if gene symbol not available
data_1$external_gene_name <- ifelse(
  data_1$external_gene_name == ".",
  data_1$gene_id,
  data_1$external_gene_name
)
# Deal with duplicated gene symbols
# Combine gene symbol with ENSEMBL ID if non-unique
data_1$external_gene_name <- ifelse(
  duplicated(data_1$external_gene_name),
  paste(data_1$external_gene_name, data_1$gene_id, sep="_"),
  data_1$external_gene_name
)
# Then we can use the gene symbol column as the row names,
# Set gene names as a separate column, not as row names
# Set gene names as a separate column, not as row names
#gene_names <- data$external_gene_name

# Subset the count data for further analysis
rownames(data_1) <- data_1$external_gene_name
counts.data <- data_1[,5:ncol(data_1)] # All columns after 4 are count data

#I had this code here to change the column names so they were easier to interpret but I am just leaving it off for now 
#columns = c("1","2","3","4","5","6","7","9","10","11","12","13","14","15","16","17","19")

#This would be to change the column name using the vector we made above
#colnames(counts.data) = columns

# Add a column named "gene" to count.data containing gene names
#count.data <- cbind(gene = gene_names, count.data)

# Now count.data has "gene" as the first column containing gene names
summary(counts.data)
###0.Plotting Data ----

#basic format for ggplot
# ggplot(data, aes(x = variable, y = variable1)) +
#   geom_col()

### 1. Barplot ----
par(mar=c(7,4,4,1)+.01) #to fix margins 
barplot(colSums(counts.data)/1e6, las=3) #las is to get the labels to fit


### 2. Histogram ----
#og2() histograms are a valuable visualization tool in RNA-seq data analysis for exploring the distributional properties of gene expression levels
hist(counts.data$X10653.SS.1, br=100) #br is so number of bins

#taking the histogram of our log data instead since it is a much better way to visualize the date - we will have some reads that have high counts and thats why that first histogram is so ugly since we have a lot without counts that are mudying the data so to speak 
logCountData = log2(1+counts.data)
hist(logCountData$X10653.SS.1, br=100)

##this code is to put all the histograms against each other to see how they fair - maybe will not be super useful##
## still editing 
# Check the column names of counts.data
#colnames(counts.data)

# Define the sample numbers
#sample_numbers <- c("X10653.SS.1", "X10653.SS.10", "X10653.SS.11", "X10653.SS.12", "X10653.SS.13", "X10653.SS.14", "X10653.SS.15", "X10653.SS.16", "X10653.SS.17", "X10653.SS.19", "X10653.SS.2", "X10653.SS.3", "X10653.SS.4", "X10653.SS.5", "X10653.SS.6", "X10653.SS.7", "X10653.SS.9")

# Define the number of samples
#num_samples <- length(sample_numbers)

# Set up the layout for the plots
#num_rows <- 4  # Adjust as needed
#num_cols <- 5  # Adjust as needed
#layout(matrix(1:(num_rows * num_cols), nrow = num_rows, byrow = TRUE))

# Define colors for each sample
#sample_colors <- rainbow(num_samples)

#par(mar = c(2, 2, 1, 1))

# Iterate through all combinations of samples
#for (i in 1:num_samples) {
  #for (j in 1:num_samples) {
    # Create a new plot with different colors
    #hist(counts.data[, sample_numbers[i]], main = paste("Sample", sample_numbers[i], "vs Sample", sample_numbers[j]),
         #xlab = "", ylab = "", col = sample_colors[j], xlim = c(0, max(counts.data)))
  #}
#}

# Print column names of counts.data
colnames(counts.data)

# Print sample_numbers
sample_numbers

# Check the structure of counts.data
str(counts.data)

#plotting replicates against each other
#perfect replicate: meaning replicate against itself as a check for what is considered "ideal"
#plot(logCountData[,1], logCountData[,1])

plot(logCountData[,1], logCountData[,3])

# Define the sample pairs
sample_pairs <- c("X10653.SS.1", "X10653.SS.10", "X10653.SS.11", "X10653.SS.12", "X10653.SS.13", "X10653.SS.14", "X10653.SS.15", "X10653.SS.16", "X10653.SS.17", "X10653.SS.19", "X10653.SS.2", "X10653.SS.3", "X10653.SS.4", "X10653.SS.5", "X10653.SS.6", "X10653.SS.7", "X10653.SS.9")

focus_sample <- "X10653.SS.19"

# Create scatter plots for the focus sample against all others
par(mfrow = c(4, 4))  # Adjust the layout as needed
for (sample2 in sample_pairs) {
  if (sample2 != focus_sample) {
    # Print selected columns
    print(paste("Selected columns:", focus_sample, sample2))
    
    # Plot scatter plot for the focus sample against the current sample
    plot(logCountData[, focus_sample], logCountData[, sample2], 
         main = paste("Sample", focus_sample, "vs Sample", sample2),
         xlab = paste("Sample", focus_sample), ylab = paste("Sample", sample2),
         col = "blue", pch = 16)
  }
}

### 3. Box Plot -----
#boxplots are a useful tool for exploring the distributional properties of gene expression data and identifying potential outliers or patterns of interest.

par(mar=c(8,4,4,1)+.01) #to fix margins 

boxplot(logCountData, las=3)

### 4. Density PLots ----
#density plots offer a more nuanced understanding of the distributional properties of gene expression data, making them a valuable complement to other visualization techniques like boxplots. They are especially useful for exploring the shape, smoothness, and comparative aspects of expression level distributions in RNA-seq datasets.

col1 =rgb(0,0,1,0.5)
col2 = rgb(1,0,0,0.5)

d1 = density(logCountData[,1])
d2 = density(logCountData[,2])
plot(d1, col='blue')
polygon(d1, col=col1, border = 'blue')
lines(d2, col='red')
polygon(d2, col=col2, border = 'red')

##loop to make the density plots for all the permutations of one sample i.e. sample 1 to 2, sample 1 to 3, sample 1 to 4, ect.
focus_sample <- "X10653.SS.19"

# Create density plots for the focus sample against all others
par(mfrow = c(4, 4))  # Adjust the layout as needed
for (sample2 in sample_pairs) {
  if (sample2 != focus_sample) {
    # Generate density plots for the focus sample and the current sample
    d_focus <- density(logCountData[, focus_sample])
    d_other <- density(logCountData[, sample2])
    
    # Determine plot range
    plot_range <- range(d_focus$x, d_other$x)
    
    # Plot density plots for the focus sample and the current sample
    plot(d_focus$x, d_focus$y, type = 'l', col = 'blue', main = paste("Density plot for", focus_sample, "vs", sample2), xlim = plot_range)
    lines(d_other$x, d_other$y, col = 'red')
    # Add polygons under the density curves
    polygon(d_focus, col = col1, border = NA)
    polygon(d_other, col = col2, border = NA)
  }
}

