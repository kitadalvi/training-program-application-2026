# ---------------------------------------------------------

# Melbourne Bioinformatics Training Program

# This exercise to assess your familiarity with R and git. Please follow
# the instructions on the README page and link to your repo in your application.
# If you do not link to your repo, your application will be automatically denied.

# Leave all code you used in this R script with comments as appropriate.
# Let us know if you have any questions!


# You can use the resources available on our training website for help:
# Intro to R: https://mbite.org/intro-to-r
# Version Control with Git: https://mbite.org/intro-to-git/

# ----------------------------------------------------------

# Load libraries -------------------
# You may use base R or tidyverse for this exercise

# ex. library(tidyverse)

library(dplyr)
library(tidyr)
library(ggplot2)
library(edgeR)

# Load data here ----------------------
# Load each file with a meaningful variable name.
metadata <- read.csv('data/GSE60450_filtered_metadata.csv', header=T, row.names = 1)
expr_matrix <- read.csv('data/GSE60450_GeneLevel_Normalized(CPM.and.TMM)_data.csv', header = T, row.names = 1)

# Inspect the data -------------------------

# What are the dimensions of each data set? (How many rows/columns in each?)
# Keep the code here for each file.

## Expression data
dim(expr_matrix)
cat(paste("The dimensions of the expression matrix are:","\n","Rows:", dim(expr_matrix)[1],"\n","Columns:", dim(expr_matrix)[2]))

## Metadata
dim(metadata)
cat(paste("The dimensions of the metadata are:","\n","Rows:", dim(metadata)[1],"\n","Columns:", dim(metadata)[2]))

# Prepare/combine the data for plotting ------------------------
# How can you combine this data into one data.frame?

### Creating a DGEList object to store count and sample information
y <- DGEList(counts = expr_matrix[,2:13], genes = expr_matrix[,1])

### Log transforming CPM values
y$counts <- cpm(y$counts, log = T)

### adding metadata to sample information by matching rownames
y$samples <- y$samples %>%
  merge(metadata, by=0)
  
# Plot the data --------------------------
## Plot the expression by cell type
## Can use boxplot() or geom_boxplot() in ggplot2

### Transpose data
expr_matrix_t <- t(data.frame(y))
colnames(expr_matrix_t) <- expr_matrix_t[1,]
expr_matrix_t <- expr_matrix_t[-1,]

### select cell types 
cell_types <- c(y$samples$immunophenotype)
names(cell_types) <- rownames(expr_matrix_t)

### bind cell types to expr matrix dataframe 
df <- cbind(expr_matrix_t,cell_types)

### change to long format for plotting
df_long <- pivot_longer(data.frame(df),
                                 cols= -cell_types,
                                 names_to = "gene",
                                 values_to = "expression")

### checking variable types 
df_long$cell_types <- as.factor(df_long$cell_types)
df_long$expression <- as.numeric(df_long$expression)
df_long$expression <- log2(df_long$expression)

### plot expression by cell type
plot <- ggplot(df_long, aes(x = cell_types, y = expression, fill = cell_types)) +
  stat_boxplot(geom ='errorbar') + 
  geom_boxplot()+
  labs(title = "Total Gene Expression Distribution per Cell Type",
       x = "Cell Type",
       y = "log2(CPM) expression") +
  theme_classic()

plot

## Save the plot
### Show code for saving the plot with ggsave() or a similar function

ggsave(plot = plot, 
       path = 'results/',
       filename = 'boxplot.pdf')
