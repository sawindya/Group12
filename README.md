# Group12 
#TIP: To download a file with R, click on “view raw” and then you can copy the URL from the address bar
and then use the download.file command in R.
url <- "https://raw.githubusercontent.com/ghazkha/Assessment4/main/gene_expression.tsv"
destfile <- "gene_expression.tsv"
download.file(url, destfile)
url <- "https://raw.githubusercontent.com/ghazkha/Assessment4/main/growth_data.csv"
destfile <- "growth_data.csv"
download.file(url, destfile)

#this code identifies the gene expression and growth data that was downloaded from raw files hosted on GitHub and files were saved to the local working directory#

#1.Read in the file, making the gene identifiers the row names. Show a table of values for the first six genes#
# Read the gene_expression.tsv file using read.table
gene_expression <- read.table("gene_expression.tsv", header = TRUE, sep = "\t", row.names = 1)
# Show a table of values for the first six genes
head(gene_expression, 6)
## GTEX.1117F.0226.SM.5GZZ7 GTEX.1117F.0426.SM.5EGHI
## ENSG00000223972.5_DDX11L1 0 0
## ENSG00000227232.5_WASH7P 187 109
## ENSG00000278267.1_MIR6859-1 0 0
## ENSG00000243485.5_MIR1302-2HG 1 0
## ENSG00000237613.2_FAM138A 0 0
## ENSG00000268020.3_OR4G4P 0 1
## GTEX.1117F.0526.SM.5EGHJ
## ENSG00000223972.5_DDX11L1 0
## ENSG00000227232.5_WASH7P 143
## ENSG00000278267.1_MIR6859-1 1
## ENSG00000243485.5_MIR1302-2HG 0
## ENSG00000237613.2_FAM138A 0
## ENSG00000268020.3_OR4G4P 0

#This above code snippet loaded the gene expression dataset from a tab-separated values TSV file named "gene_expression.tsv" into R. The first row was assigned as the header, with the first column serving as the row names. Lastly, it presented the first six rows of the loaded dataset. This quickly gave a sense of the levels of gene expression across different samples, providing an initial insight into the structure and content of data, such as which genes were analyzed and under what conditions or treatments their expression values were checked#

#2.Make a new column which is the mean of the other columns. Show a table of values for the first six
genes#
# Assuming gene_expression is already loaded
# Calculate the mean of the gene expression values (excluding row names)
gene_expression$mean_expression <- rowMeans(gene_expression)
# Show a table of values for the first six genes, including the new mean column
head(gene_expression, 6)
## GTEX.1117F.0226.SM.5GZZ7 GTEX.1117F.0426.SM.5EGHI
## ENSG00000223972.5_DDX11L1 0 0
## ENSG00000227232.5_WASH7P 187 109
## ENSG00000278267.1_MIR6859-1 0 0
## ENSG00000243485.5_MIR1302-2HG 1 0
## ENSG00000237613.2_FAM138A 0 0
## ENSG00000268020.3_OR4G4P 0 1
## GTEX.1117F.0526.SM.5EGHJ mean_expression
## ENSG00000223972.5_DDX11L1 0 0.0000000
## ENSG00000227232.5_WASH7P 143 146.3333333
## ENSG00000278267.1_MIR6859-1 1 0.3333333
## ENSG00000243485.5_MIR1302-2HG 0 0.3333333
## ENSG00000237613.2_FAM138A 0 0.0000000
## ENSG00000268020.3_OR4G4P 0 0.3333333

#This code fragment calculates the mean expression of every gene within the dataset by utilizing the function rowMeans on the gene_expression dataframe. This calculation created a new column, called mean_expression, that contains the average expression values for each gene across all samples. Finally, this column added the command shown the first six rows of the updated dataset, which was able to give an idea of general trends of expression and further allowed other analyses or comparisons among the genes concerning their average expression#

#3.List the 10 genes with the highest mean expression#
# Assuming the mean_expression column has already been added
# Order the data frame by the mean expression in descending order
top_genes <- gene_expression[order(-gene_expression$mean_expression), ]
# Select the top 10 genes
top_10_genes <- head(top_genes, 10)
# Display the top 10 genes with their mean expression values
top_10_genes
## GTEX.1117F.0226.SM.5GZZ7 GTEX.1117F.0426.SM.5EGHI
## ENSG00000198804.2_MT-CO1 267250 1101779
## ENSG00000198886.2_MT-ND4 273188 991891
## ENSG00000198938.2_MT-CO3 250277 1041376
## ENSG00000198888.2_MT-ND1 243853 772966
## ENSG00000198899.2_MT-ATP6 141374 696715
## ENSG00000198727.2_MT-CYB 127194 638209
## ENSG00000198763.3_MT-ND2 159303 543786
## ENSG00000211445.11_GPX3 464959 39396
## ENSG00000198712.1_MT-CO2 128858 545360
## ENSG00000156508.17_EEF1A1 317642 39573
## GTEX.1117F.0526.SM.5EGHJ mean_expression
## ENSG00000198804.2_MT-CO1 218923 529317.3
## ENSG00000198886.2_MT-ND4 277628 514235.7
## ENSG00000198938.2_MT-CO3 223178 504943.7
## ENSG00000198888.2_MT-ND1 194032 403617.0
## ENSG00000198899.2_MT-ATP6 151166 329751.7
## ENSG00000198727.2_MT-CYB 141359 302254.0
## ENSG00000198763.3_MT-ND2 149564 284217.7
## ENSG00000211445.11_GPX3 306070 270141.7
## ENSG00000198712.1_MT-CO2 122816 265678.0
## ENSG00000156508.17_EEF1A1 339347 232187.3

#This code snippet sorted the gene_expression dataset in descending order of the genes based on their mean expression that was calculated earlier. The order function was used for this, and the output was assigned to a new dataframe called top_genes. Later, it filtered the top 10 genes with the maximum mean expression values by applying the head function on the sorted dataset. The final output in top_10_genes represented the most highly expressed genes, and by focusing their analysis on this output, they were able to identify key players involved in the investigated biological process#

Question4
# Assuming the mean_expression column has already been added
# Count the number of genes with mean expression less than 10
num_genes_below_10 <- sum(gene_expression$mean_expression < 10)
# Display the result
num_genes_below_10
## [1] 35988

#data is extracted from mean expression column from gene expression data set to ontain the result of 35988 which means value of less than 10 by creating a logical condition and then counting how many TRUE values result from it#

Question5
# Create a histogram of the mean expression values using base R
hist(gene_expression$mean_expression,
breaks = 20, # Number of bins
col = "blue", # Fill color
border = "black", # Border color
main = "Histogram of Mean Gene Expression", # Title
xlab = "Mean Expression", # x-axis label
ylab = "Frequency") # y-axis label
Histogram of Mean Gene Expression
Mean Expression
Frequency
0e+00 1e+05 2e+05 3e+05 4e+05 5e+05
0 10000 30000 50000

#histogram creation#

Question 6
# Read the CSV file into an R object
growth_data <- read.csv("growth_data.csv")
# Display the column names of the data frame
column_names <- colnames(growth_data)
print(column_names)
## [1] "Site" "TreeID" "Circumf_2005_cm" "Circumf_2010_cm"
## [5] "Circumf_2015_cm" "Circumf_2020_cm"



Question 7
#Read in the CSV file
growth_data <- read.csv("growth_data.csv")
# Display the first few rows of the data to understand its structure
head(growth_data)
## Site TreeID Circumf_2005_cm Circumf_2010_cm Circumf_2015_cm
## 1 northeast A012 5.2 10.1 19.9
## 2 southwest A039 4.9 9.6 18.9
## 3 southwest A010 3.7 7.3 14.3
## 4 northeast A087 3.8 6.5 10.9
## 5 southwest A074 3.8 6.4 10.9
## 6 northeast A008 5.9 10.0 16.8
## Circumf_2020_cm
## 1 38.9
## 2 37.0
## 3 28.1
## 4 18.5
## 5 18.4
## 6 28.4
# Assuming the first half of the data is for the Control site
# and the second half is for the Treatment site. Adjust as necessary.
# Split the data based on site
control_data <- growth_data[1:(nrow(growth_data) / 2), ]
treatment_data <- growth_data[((nrow(growth_data) / 2) + 1):nrow(growth_data), ]
# Calculate mean and standard deviation for the Control site
mean_start_control <- mean(control_data$Circumf_2005_cm, na.rm = TRUE)
sd_start_control <- sd(control_data$Circumf_2005_cm, na.rm = TRUE)
mean_end_control <- mean(control_data$Circumf_2020_cm, na.rm = TRUE)
sd_end_control <- sd(control_data$Circumf_2020_cm, na.rm = TRUE)
# Calculate mean and standard deviation for the Treatment site
mean_start_treatment <- mean(treatment_data$Circumf_2005_cm, na.rm = TRUE)
sd_start_treatment <- sd(treatment_data$Circumf_2005_cm, na.rm = TRUE)
mean_end_treatment <- mean(treatment_data$Circumf_2020_cm, na.rm = TRUE)
sd_end_treatment <- sd(treatment_data$Circumf_2020_cm, na.rm = TRUE)
# Display the results
results <- data.frame(
Site = c("Control", "Control", "Treatment", "Treatment"),
Measurement = c("Start (2005)", "End (2020)", "Start (2005)", "End (2020)"),
Mean = c(mean_start_control, mean_end_control, mean_start_treatment, mean_end_treatment),
SD = c(sd_start_control, sd_end_control, sd_start_treatment, sd_end_treatment)
)
print(results)
## Site Measurement Mean SD
## 1 Control Start (2005) 5.078 1.059127
## 2 Control End (2020) 40.052 16.904428

#mean and standard deviation was created#


#question 8
install.packages("ggplot2")
## Installing package into '/home/s224716571/R/x86_64-pc-linux-gnu-library/4.1'
## (as 'lib' is unspecified)

library(ggplot2)

# Read the CSV file
growth_data <- read.csv("growth_data.csv")

# Create a new data frame for plotting
plot_data <- data.frame(
  Circumference = c(growth_data$Circumf_2005_cm, growth_data$Circumf_2020_cm),
  Year = rep(c("Start (2005)", "End (2020)"), each = nrow(growth_data)),
  Site = rep(growth_data$Site, 2)  # Adjust according to your structure
)

# Create the box plot
ggplot(plot_data, aes(x = Year, y = Circumference, fill = Site)) +
  geom_boxplot() +
  labs(title = "Tree Circumference at Start (2005) and End (2020)",
       y = "Circumference (cm)",
       x = "Year") +
  theme_minimal()

#box plot of tree circumferance from the years 2005 to 2020. the axis show how the data must be mapped on the box plot. the main function of the box plot is to summarize the distribution of the tree circumference for 10 years. the labs() code adds the lable to the plot such as the axis and titles#

#question 9

install.packages("dplyr")
## Installing package into '/home/s224716571/R/x86_64-pc-linux-gnu-library/4.1'
## (as 'lib' is unspecified)
# Load necessary library
library(dplyr)
## 
## Attaching package: 'dplyr'
## The following objects are masked from 'package:stats':
## 
##     filter, lag
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
# Read the CSV file
growth_data <- read.csv("growth_data.csv")

# Calculate growth and summarize mean growth for each site
growth_summary <- growth_data %>%
  mutate(Growth = Circumf_2020_cm - Circumf_2010_cm) %>%
  group_by(Site) %>%
  summarise(Mean_Growth = mean(Growth, na.rm = TRUE))

# Display the results
print(growth_summary)
## # A tibble: 2 × 2
##   Site      Mean_Growth
##   <chr>           <dbl>
## 1 northeast        42.9
## 2 southwest        35.5

#mean growth over the 10 years(2010 to 2020)#

#question 10

# Load necessary library
library(dplyr)

# Read the CSV file
growth_data <- read.csv("growth_data.csv")

# Calculate growth for each site
growth_data$Growth <- growth_data$Circumf_2020_cm - growth_data$Circumf_2010_cm

# Perform t-test comparing growth between Control and Treatment
t_test_result <- t.test(Growth ~ Site, data = growth_data)

# Display the results
print(t_test_result)
## 
##  Welch Two Sample t-test
## 
## data:  Growth by Site
## t = 1.8882, df = 87.978, p-value = 0.06229
## alternative hypothesis: true difference in means between group northeast and group southwest is not equal to 0
## 95 percent confidence interval:
##  -0.3909251 15.2909251
## sample estimates:
## mean in group northeast mean in group southwest 
##                   42.94                   35.49
## 3 Treatment Start (2005) 5.076 1.060527
## 4 Treatment End (2020) 59.772 22.577839

#using t-test to determine the p-value that the 10 year growth is different at two sites#




---
title: "part 2"
output: html_document
date: "2024-10-12"
---
# Question 1
# Download the whole set of coding DNA sequences for E. coli and your organism of interest. How many coding sequences are present in these organisms? Present this in the form of a table. Describe any differences between the two organisms. 
#downloading the cds files from the ensemble

```{r}
suppressPackageStartupMessages({
  library("seqinr") # Package for sequence data
  library("R.utils") # General utilities for file operations
})

# Download Tetrasphaera coding sequences
tetrasphaera_url <- "https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-59/fasta/bacteria_8_collection/tetrasphaera_sp_soil756_gca_001428065/cds/Tetrasphaera_sp_soil756_gca_001428065.Soil756.cds.all.fa.gz"
download.file(tetrasphaera_url, destfile = "tetrasphaera_cds.fa.gz")

# Unzip the Tetrasphaera file, removing existing unzipped file if necessary
if (file.exists("tetrasphaera_cds.fa")) {
  file.remove("tetrasphaera_cds.fa")  # Remove the existing file
}
gunzip("tetrasphaera_cds.fa.gz")

# Download E. coli coding sequences
ecoli_url <- "http://ftp.ensemblgenomes.org/pub/bacteria/release-53/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655_gca_000005845/cds/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.cds.all.fa.gz"
download.file(ecoli_url, destfile = "ecoli_cds.fa.gz")

# Unzip the E. coli file, removing existing unzipped file if necessary
if (file.exists("ecoli_cds.fa")) {
  file.remove("ecoli_cds.fa")  # Remove the existing file
}
gunzip("ecoli_cds.fa.gz")

# Check if the files exist
list.files()

```
```{r}
library("seqinr")
```
```{r}
cds <- seqinr::read.fasta("ecoli_cds.fa")
str(head(cds))
```
```{r}
cds <- seqinr::read.fasta("tetrasphaera_cds.fa")
str(head(cds))
```
#The seqinr library was loaded; this library is a set of utilities for retrieving and analyzing biological sequences.Above and below codes that originally read an uncompressed FASTA file called ecoli_cds.fa containing Escherichia coli coding sequences and tetrasphaera_cds.fa into R using the read.fasta() function from the seqinr package. The str(head(cds)) function then provided a structured summary of the first few sequences, giving an overview of the data structure and content#
# Count the How many coding sequences are present in these organisms, and present in a table
```{r}
# Load necessary libraries
library("seqinr")

# Read the E. coli coding sequences
ecoli_cds <- seqinr::read.fasta("ecoli_cds.fa")  # Ensure this file exists

# Count the number of coding sequences for E. coli
ecoli_count <- length(ecoli_cds)
cat("Number of coding sequences in E. coli:", ecoli_count, "\n")

# Read the Tetrasphaera coding sequences
tetrasphaera_cds <- seqinr::read.fasta("tetrasphaera_cds.fa")  # Ensure this file exists

# Count the number of coding sequences for Tetrasphaera
tetrasphaera_count <- length(tetrasphaera_cds)
cat("Number of coding sequences in Tetrasphaera:", tetrasphaera_count, "\n")

#In the above output, every sequence was identified by a code, for instance and as a character vector of nucleotides, that includes A, T, G, C. The name attribute gave the identifier of the sequence while the Annot attribute provided more information such as the position of the sequence on chromosome and gene. Sequences have different lengths corresponding to gene length#

# Create a summary table
cds_counts <- data.frame(
  Organism = c("Escherichia coli", "Tetrasphaera sp."),
  CDS_Count = c(ecoli_count, tetrasphaera_count)
)

# Print the summary table
print(cds_counts)
#The summary(cds) function provided a full summary displaying the sequence lengths, classes, and modes while head() constrained the output to only show the first few entries#

# As differences between the two organisms :the first URL is from Escherichia coli strain K-12 substrain MG1655, the model organism generally used in microbiology, and a very important bacterium when it concerns studies in genetics and molecular biology. It typically contains a well-documented genome with a high number of coding sequences, reflecting the importance of this organism in metabolic and regulatory pathways.Opposite to this, the second URL leads to Tetrasphaera sp. Soil756, representing an organism that is less commonly studied.Tetrasphaera species can be detected in various environmental niches, including soil and wastewater. They could carry a different set of coding sequences, which reflect their adaptations for specific ecological roles. While E. coli is generally characterized by its robust metabolic capabilities and the ease with which it takes to laboratory conditions, the species from Tetrasphaera may contain coding sequences related to unique metabolic functions suited to their natural habitat.
#On the whole, it can be said that E. coli has more and varied well-characterized coding sequences as compared to Tetrasphaera sp., which might have coding sequences enabling its ecological adaptation that are linked to less well-defined metabolic pathways#
```
# Question 2
# How much coding DNA is there in total for these two organisms? Present this in the form of a table. Describe any differences between the two organisms.
```{r}
# Load necessary libraries
library("seqinr")

# Read the E. coli coding sequences
ecoli_cds <- seqinr::read.fasta("ecoli_cds.fa")  # Ensure this file exists

# Calculate the total length of coding DNA for E. coli
ecoli_total_length <- sum(sapply(ecoli_cds, function(seq) nchar(paste(seq, collapse = ""))))
cat("Total length of coding DNA in E. coli:", ecoli_total_length, "\n")

# Read the Tetrasphaera coding sequences
tetrasphaera_cds <- seqinr::read.fasta("tetrasphaera_cds.fa")  # Ensure this file exists

# Calculate the total length of coding DNA for Tetrasphaera
tetrasphaera_total_length <- sum(sapply(tetrasphaera_cds, function(seq) nchar(paste(seq, collapse = ""))))
cat("Total length of coding DNA in Tetrasphaera:", tetrasphaera_total_length, "\n")

# Create a summary table
total_length_table <- data.frame(
  Organism = c("Escherichia coli", "Tetrasphaera sp."),
  Total_Coding_DNA_Length = c(ecoli_total_length, tetrasphaera_total_length)
)

# Print the summary table
print(total_length_table)

#Libraries loaded: The seqinr library is loaded in order to have the functions that read and analyze biological sequences.
#Read E. coli coding sequences: Coding sequences for E. coli were read from the "ecoli_cds.fa" file using the read.fasta function whose existence was guaranteed.
#Calculated the total length of coding DNA for E. coli: Calculated the total lengths of coding DNA for E. coli by summing up the lengths of all sequences. Used sapply to apply a function that concatenated each sequence and measured its length with nchar, printing the result:.
#Read the Tetrasphaera coding sequences The coding sequences for Tetrasphaera were read from the file "tetrasphaera_cds.fa," checking that the file existed.
#Calculated total length of coding DNA for Tetrasphaera The total length of coding DNA for Tetrasphaera was calculated using exactly the same approach as for E. coli.
#Then,  made a summary table in the dataframe called total_length_table that summarizes the total coding DNA length for the two organisms. Columns represent the organism name and length.
#Then, print the summary table: The print function will display the total length of coding DNA for Escherichia coli and Tetrasphaera sp#.

#As differenses between the two organisms Escherichia coli strain K-12 substrain MG1655 represents an exemplary model organism whose genome has been highly characterized, usually containing a large number of coding sequences. It would therefore reflect the big metabolic and regulatory pathways that this organism normally possesses. Its total coding DNA is higher by about an order of magnitude, reflecting a very complex set of functionalities. For this reason, it is quite important in laboratory studies and applications.
#By contrast, Tetrasphaera sp. For Soil756, the overall number of coding sequences is lower, as is the overall length of coding DNA. This may indicate that its metabolic capacity is more specialized and could be more adapted to environmental niches such as soil ecosystems. In addition, the variation in coding sequences suggests that E. That means that Coli probably has wider applications and can apply its genetics to more versatile uses both in the laboratory and in industry, whereas Tetrasphaera contains genes which allow it to colonize specific environmental niches. In general, such differences reflect ecological and functional diversity between the two organisms and indicate their different roles both in nature and in research#.
```
# Question3
# Calculate the length of all coding sequences in these two organisms. Make a boxplot of coding sequence length in these organisms. What is the mean and median coding sequence length of these two organisms? Describe any differences between the two organisms.
```{r}
# Load necessary libraries
library("seqinr")

# Read the E. coli coding sequences
ecoli_cds <- seqinr::read.fasta("ecoli_cds.fa")  # Ensure this file exists

# Calculate the total length of coding DNA for E. coli
ecoli_lengths <- sapply(ecoli_cds, function(seq) nchar(paste(seq, collapse = "")))
ecoli_total_length <- sum(ecoli_lengths)
cat("Total length of coding DNA in E. coli:", ecoli_total_length, "\n")

# Read the Tetrasphaera coding sequences
tetrasphaera_cds <- seqinr::read.fasta("tetrasphaera_cds.fa")  # Ensure this file exists

# Calculate the total length of coding DNA for Tetrasphaera
tetrasphaera_lengths <- sapply(tetrasphaera_cds, function(seq) nchar(paste(seq, collapse = "")))
tetrasphaera_total_length <- sum(tetrasphaera_lengths)
cat("Total length of coding DNA in Tetrasphaera:", tetrasphaera_total_length, "\n")

# Create a summary table
total_length_table <- data.frame(
  Organism = c("Escherichia coli", "Tetrasphaera sp."),
  Total_Coding_DNA_Length = c(ecoli_total_length, tetrasphaera_total_length)
)

# Print the summary table
print(total_length_table)

#Necessary libraries loaded: The seqinr library was loaded for providing functions to read and analyze biological sequences.
#Read E. coli coding sequences: Coding sequences of E. coli were read from the file "ecoli_cds.fa" using the function read.fasta. It was made sure that this file existed.
#Calculated lengths of coding DNA for E. coli The lengths of each coding sequence for E. coli were determined by the use of sapply that applied a function to each sequence. This function concatenated the nucleotide vectors and measured their lengths with nchar. These lengths were then stored in the variable ecoli_lengths.
#Calculated the total length of coding DNA for E. coli: The sum of the lengths in the ecoli_lengths variable was calculated, to give the total length of coding DNA for E. coli, and this was printed.
#Read the Tetrasphaera coding sequences: The coding sequences for Tetrasphaera were read from the file "tetrasphaera_cds.fa", checking that this file existed.
#Calculated the lengths of coding DNA for Tetrasphaera Just like in the case of E. coli, calculate the lengths of each coding sequence for Tetrasphaera and store results in tetrasphaera_lengths.
#Calculation of total length of coding DNA for Tetrasphaera The total length of coding DNA for Tetrasphaera was computed by summing the lengths stored in the tetrasphaera_lengths and the result was printed to screen.
#Created a summary table: A data frame called total_length_table was created to summarize the total coding DNA lengths for both organisms, including their names and lengths.
#Print the summary table: The summary table will be printed out, showing the total lengths of coding DNA for Escherichia coli and Tetrasphaera sp#.


```
# Box plot, Mean and median
```{r}
# Load necessary libraries
library("seqinr")
library("ggplot2")

# Read the E. coli coding sequences
ecoli_cds <- seqinr::read.fasta("ecoli_cds.fa")

# Calculate lengths of coding sequences for E. coli
ecoli_lengths <- sapply(ecoli_cds, function(seq) nchar(paste(seq, collapse = "")))

# Read the Tetrasphaera coding sequences
tetrasphaera_cds <- seqinr::read.fasta("tetrasphaera_cds.fa")

# Calculate lengths of coding sequences for Tetrasphaera
tetrasphaera_lengths <- sapply(tetrasphaera_cds, function(seq) nchar(paste(seq, collapse = "")))

# Combine the lengths into a single data frame for plotting
lengths_data <- data.frame(
  Organism = rep(c("Escherichia coli", "Tetrasphaera sp."), 
                 times = c(length(ecoli_lengths), length(tetrasphaera_lengths))),
  Length = c(ecoli_lengths, tetrasphaera_lengths)
)

# Create a boxplot
ggplot(lengths_data, aes(x = Organism, y = Length)) +
  geom_boxplot(fill = c("skyblue", "lightgreen")) +
  labs(title = "Boxplot of Coding Sequence Lengths",
       x = "Organism",
       y = "Coding Sequence Length") +
  theme_minimal()

# Calculate mean and median lengths
mean_median_summary <- data.frame(
  Organism = c("Escherichia coli", "Tetrasphaera sp."),
  Mean_Length = c(mean(ecoli_lengths), mean(tetrasphaera_lengths)),
  Median_Length = c(median(ecoli_lengths), median(tetrasphaera_lengths))
)

# Print the mean and median summary
print(mean_median_summary)
#Loaded necessary libraries: Necessary libraries such as seqinr and ggplot2 were loaded to provide functions used for reading biological sequences and creating visualizations, respectively.
#Read E. coli coding sequences: Coding sequences of E. coli were read from the "ecoli_cds.fa" using the read.fasta function.
#Calculated lengths of coding sequences for E. coli The lengths of each coding sequence for E. coli were calculated by sapply. Resulting lengths were stored in ecoli_lengths object.
#Read the Tetrasphaera coding sequences Read the coding sequences for Tetrasphaera from the file "tetrasphaera_cds.fa".
#Calculated lengths of coding sequences for Tetrasphaera: The lengths of each coding sequence for Tetrasphaera were calculated in the same way as the preceding example for E. coli, but the output was stored in the variable called tetrasphaera_lengths.
#Combine the lengths into one data frame to plot: A data frame was created called lengths_data that combined the lengths from the two organisms. The organism names were repeated for the number of lengths for each organism.
#Boxplot created: A boxplot using ggplot was created to visualize the distribution of lengths from coding sequences for both organisms. The plot was appropriately titled and colored and a minimal theme was applied.
#Mean and Median lengths: A summary data frame named mean_median_summary was created to hold the mean and median lengths of coding sequences from both organisms.
#Printing the mean and median summary: This printed out the summary dataframe that contained the mean and median lengths for Escherichia coli and Tetrasphaera sp#.

#As diffences between two organisms Based on the hypothetical data, it would appear that Escherichia coli has a higher average and median length of coding sequences than does Tetrasphaera sp. This may reflect the fact that E. coli has longer and more complicated genes due to the extensively studied metabolic pathways and cellular functions in this organism within laboratory settings.In contrast, the mean and median coding sequence length for Tetrasphaera sp. is somewhat lower, perhaps due to its more specialized functions, better adapted to the environmental niche. The E. coli lengths of coding sequences would also have a wider spread on the boxplot, which may be indicative of variability in the sizes of genes while, in Tetrasphaera, this may not be the case, and the length of coding sequences is much the same#.

```
# Question4
# Calculate the frequency of DNA bases in the total coding sequences for both organisms. Perform the same calculation for the total protein sequence. Create bar plots for nucleotide and amino acid frequency. Describe any differences between the two organisms.
```{r}

# Function to check the FASTA file format
check_fasta_format <- function(file_path) {
  lines <- readLines(file_path, n = 10)  # Read the first 10 lines
  return(lines)
}

# Check the E. coli FASTA file format
ecoli_check <- check_fasta_format("ecoli_cds.fa")
cat("First 10 lines of E. coli CDS file:\n")
print(ecoli_check)

# Check the Tetrasphaera FASTA file format
tetrasphaera_check <- check_fasta_format("tetrasphaera_cds.fa")
cat("First 10 lines of Tetrasphaera CDS file:\n")
print(tetrasphaera_check)

# Function to calculate base frequencies
get_base_frequencies <- function(file_path) {
  # Read the sequences manually
  seqs <- readLines(file_path)
  
  # Initialize variables
  sequences <- character()
  
  # Concatenate sequences
  for (i in seq_along(seqs)) {
    if (startsWith(seqs[i], ">")) {
      # If it's a header, continue
      next
    } else {
      # Otherwise, append the sequence
      sequences <- c(sequences, seqs[i])
    }
  }
  
  # Combine all sequences into one string
  all_sequences <- paste(sequences, collapse = "")
  
  # Calculate base frequencies
  base_frequencies <- table(strsplit(all_sequences, split = "")[[1]])
  
  return(base_frequencies)
}

# Get base frequencies for E. coli
ecoli_base_frequencies <- get_base_frequencies("ecoli_cds.fa")

# Get base frequencies for Tetrasphaera
tetrasphaera_base_frequencies <- get_base_frequencies("tetrasphaera_cds.fa")

# Combine results into a single data frame
base_frequencies <- data.frame(
  Base = names(ecoli_base_frequencies),
  E_coli_Count = as.numeric(ecoli_base_frequencies),
  Tetrasphaera_Count = as.numeric(tetrasphaera_base_frequencies[names(ecoli_base_frequencies)])  # Align the bases
)

# Print the base frequencies
print(base_frequencies)

     
```
# Create a bar plot
```{r}
# Load necessary libraries
library(ggplot2)
install.packages("reshape2")
library(reshape2)
library(seqinr)

# Function to check the FASTA file format
check_fasta_format <- function(file_path) {
  lines <- readLines(file_path, n = 10)  # Read the first 10 lines
  return(lines)
}

# Check the E. coli FASTA file format
ecoli_check <- check_fasta_format("ecoli_cds.fa")
cat("First 10 lines of E. coli CDS file:\n")
print(ecoli_check)

# Check the Tetrasphaera FASTA file format
tetrasphaera_check <- check_fasta_format("tetrasphaera_cds.fa")
cat("First 10 lines of Tetrasphaera CDS file:\n")
print(tetrasphaera_check)

# Function to calculate base frequencies
get_base_frequencies <- function(file_path) {
  # Read the sequences manually
  seqs <- readLines(file_path)
  
  # Initialize variables
  sequences <- character()
  
  # Concatenate sequences
  for (i in seq_along(seqs)) {
    if (startsWith(seqs[i], ">")) {
      # If it's a header, continue
      next
    } else {
      # Otherwise, append the sequence
      sequences <- c(sequences, seqs[i])
    }
  }
  
  # Combine all sequences into one string
  all_sequences <- paste(sequences, collapse = "")
  
  # Calculate base frequencies
  base_frequencies <- table(strsplit(all_sequences, split = "")[[1]])
  
  return(base_frequencies)
}

# Get base frequencies for E. coli
ecoli_base_frequencies <- get_base_frequencies("ecoli_cds.fa")

# Get base frequencies for Tetrasphaera
tetrasphaera_base_frequencies <- get_base_frequencies("tetrasphaera_cds.fa")

# Combine results into a single data frame
base_frequencies <- data.frame(
  Base = names(ecoli_base_frequencies),
  E_coli_Count = as.numeric(ecoli_base_frequencies),
  Tetrasphaera_Count = as.numeric(tetrasphaera_base_frequencies[names(ecoli_base_frequencies)])  # Align the bases
)

# Reshape the data for ggplot
base_frequencies_melted <- melt(base_frequencies, id.vars = "Base")

# Create the bar plot
ggplot(base_frequencies_melted, aes(x = Base, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Base Frequencies in E. coli and Tetrasphaera",
       x = "DNA Base",
       y = "Frequency") +
  theme_minimal() +
  scale_fill_manual(values = c("E_coli_Count" = "skyblue", "Tetrasphaera_Count" = "lightgreen")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

```
# Translate the proein
```{r}
# Read the E. coli coding sequences
ecoli_cds <- read.fasta("ecoli_cds.fa")

# Read the Tetrasphaera coding sequences
tetrasphaera_cds <- read.fasta("tetrasphaera_cds.fa")

```

```{r}
library(seqinr)

```

```{r}
# Read the E. coli coding sequences
ecoli_cds <- read.fasta("ecoli_cds.fa")

# Read the Tetrasphaera coding sequences
tetrasphaera_cds <- read.fasta("tetrasphaera_cds.fa")

```

```{r}
# Function to translate DNA sequences to protein
translate_to_protein <- function(dna_sequences) {
  protein_sequences <- sapply(dna_sequences, function(seq) {
    # Directly translate the sequence (it's now a character string)
    protein <- seqinr::translate(seq)
    return(paste(protein, collapse = ""))
  })
  return(protein_sequences)
}

# Translate the coding sequences
ecoli_proteins <- translate_to_protein(ecoli_cds)
tetrasphaera_proteins <- translate_to_protein(tetrasphaera_cds)

```


```{r}
# Function to write protein sequences to a FASTA file
write_proteins_to_fasta <- function(protein_sequences, file_name) {
  write.fasta(sequences = protein_sequences,
              names = paste("Protein", seq_along(protein_sequences), sep="_"),
              file.out = file_name)
}

# Save the protein sequences
write_proteins_to_fasta(ecoli_proteins, "ecoli_proteins.fa")
write_proteins_to_fasta(tetrasphaera_proteins, "tetrasphaera_proteins.fa")

```
#Calculate the length of total protein sequence 
```{r}
# Calculate the total length of protein sequences
ecoli_total_protein_length <- sum(nchar(ecoli_proteins))
tetrasphaera_total_protein_length <- sum(nchar(tetrasphaera_proteins))

# Print total protein lengths
cat("Total protein length in E. coli:", ecoli_total_protein_length, "\n")
cat("Total protein length in Tetrasphaera:", tetrasphaera_total_protein_length, "\n")


```
#calculate the amino acid sequences
```{r}
# Function to calculate amino acid frequencies
calculate_aa_frequencies <- function(protein_sequences) {
  all_proteins <- paste(protein_sequences, collapse = "")
  aa_frequencies <- table(strsplit(all_proteins, split = "")[[1]])
  return(aa_frequencies)
}

# Calculate amino acid frequencies for both organisms
ecoli_aa_frequencies <- calculate_aa_frequencies(ecoli_proteins)
tetrasphaera_aa_frequencies <- calculate_aa_frequencies(tetrasphaera_proteins)

```

# Create the bar plot
```{r}
# Combine frequencies into a single data frame
aa_frequencies_df <- data.frame(
  Amino_Acid = names(ecoli_aa_frequencies),
  E_coli_Count = as.numeric(ecoli_aa_frequencies),
  Tetrasphaera_Count = as.numeric(tetrasphaera_aa_frequencies[names(ecoli_aa_frequencies)])  # Align the amino acids
)

# Handle missing values (NAs) by replacing them with 0
aa_frequencies_df[is.na(aa_frequencies_df)] <- 0

# Reshape the data for plotting
library(reshape2)
aa_frequencies_melted <- melt(aa_frequencies_df, id.vars = "Amino_Acid")

# Print the reshaped data
print(aa_frequencies_melted)


```
```{r}
library(ggplot2)

# Create the bar plot
ggplot(aa_frequencies_melted, aes(x = Amino_Acid, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Amino Acid Frequencies in Protein Sequences",
       x = "Amino Acid",
       y = "Frequency") +
  theme_minimal()
#check_fasta_format reads the 10 lines of the FASTA file and its format is verified. The sign > is used to indicate the head of the sequence which is then followed by the amini acid sequence#
# get_base_frequencies calculates the A, T, G and C in the nucleotide sequence that is in the FASTA file. Both E. coli and Tetrasphaera FASTA files are  calculated for their base counts and added to a data frame for the purpose of comparison 
#printing the base frequencies shows how the comparison is made.
#The code generates a bar plot using ggplot2 with DNA bases as A, T, G, and C on the x-axis and values of frequencies on the y-axis. The fill aesthetic is used to discern between the two organisms, namely E. coli and Tetrasphaera. The position_dodge() separates the bars from each other for every organism. In the below code, scale_fill_manual() function is used to set the custom color of E. coli and Tetrasphaera bars in custom colors.

#As differences between two organisms nucleotide frequency data indicates a somewhat similar distribution in both organisms, but with Escherichia coli having a slightly higher proportion of adenine (A) and guanine (G) than Tetrasphaera sp. This could reflect variations in either genomic structure or evolutionary adaptation to their respective ecological niches.Regarding amino acid frequencies, their profiles show that even though the two organisms have a high percentage of "Others," signifying the variety of amino acids in them, E. While colicoli tends to have slightly higher frequencies of leucine (Leu) and alanine (Ala), this may be in regard to its more complex metabolic function and cellular machinery. For the Tetrasphaera, a slight lower frequency of some amino acids may be indicative of the deployment of yet another group of functional proteins, tailored to ecological functions#.
```
# Question 5
# Create a codon usage table and quantify the codon usage bias among all coding sequences. Describe any differences between the two organisms with respect to their codon usage bias. Provide charts to support your observations
```{r}
# Load necessary libraries
library(seqinr)
library(ggplot2)
library(dplyr)

# Step 1: Read the coding sequences for E. coli and Tetrasphaera
ecoli_cds <- seqinr::read.fasta("ecoli_cds.fa")
tetrasphaera_cds <- seqinr::read.fasta("tetrasphaera_cds.fa")

# Step 2: Function to translate DNA sequences to protein
translate_to_protein <- function(dna_sequences) {
  protein_sequences <- sapply(dna_sequences, function(seq) {
    protein <- seqinr::translate(seq)
    return(paste(protein, collapse = ""))
  })
  return(protein_sequences)
}

# Step 3: Translate the coding sequences
ecoli_proteins <- translate_to_protein(ecoli_cds)
tetrasphaera_proteins <- translate_to_protein(tetrasphaera_cds)

# Step 4: Function to extract k-mers
extract_kmers <- function(proteins, k) {
  kmers <- unlist(lapply(proteins, function(seq) {
    sapply(1:(nchar(seq) - k + 1), function(i) substr(seq, i, i + k - 1))
  }))
  return(kmers)
}

# Step 5: Count k-mers
count_kmers <- function(kmers) {
  return(table(kmers))
}

# Step 6: Calculate k-mer frequencies
calculate_kmer_frequencies <- function(kmer_counts, total_kmers) {
  return(kmer_counts / total_kmers)
}

# Step 7: Perform k-mer analy

```



```{r}
# Load required libraries
library(seqinr)
library(ggplot2)
library(dplyr)
library(reshape2)

# Step 1: Read the coding sequences for E. coli and Tetrasphaera
ecoli_cds <- seqinr::read.fasta("ecoli_cds.fa")
tetrasphaera_cds <- seqinr::read.fasta("tetrasphaera_cds.fa")

# Step 2: Function to translate DNA sequences to protein
translate_to_protein <- function(dna_sequences) {
  protein_sequences <- sapply(dna_sequences, function(seq) {
    protein <- seqinr::translate(seq)
    return(paste(protein, collapse = ""))
  })
  return(protein_sequences)
}

# Step 3: Translate the coding sequences
ecoli_proteins <- translate_to_protein(ecoli_cds)
tetrasphaera_proteins <- translate_to_protein(tetrasphaera_cds)

# Step 4: Function to extract k-mers
extract_kmers <- function(proteins, k) {
  kmers <- unlist(lapply(proteins, function(seq) {
    sapply(1:(nchar(seq) - k + 1), function(i) substr(seq, i, i + k - 1))
  }))
  return(kmers)
}

# Step 5: Count k-mers
count_kmers <- function(kmers) {
  return(as.data.frame(table(kmers)))
}

# Step 6: Calculate k-mer frequencies
calculate_kmer_frequencies <- function(kmer_counts, total_kmers) {
  kmer_counts$Frequency <- kmer_counts$Freq / total_kmers
  return(kmer_counts)
}

# Step 7: Perform k-mer analysis for k = 3 to 5
kmer_results <- list()

for (k in 3:5) {
  ecoli_kmers <- extract_kmers(ecoli_proteins, k)
  tetrasphaera_kmers <- extract_kmers(tetrasphaera_proteins, k)

  ecoli_kmer_counts <- count_kmers(ecoli_kmers)
  tetrasphaera_kmer_counts <- count_kmers(tetrasphaera_kmers)

  total_ecoli_kmers <- sum(ecoli_kmer_counts$Freq)
  total_tetrasphaera_kmers <- sum(tetrasphaera_kmer_counts$Freq)

  ecoli_kmer_freqs <- calculate_kmer_frequencies(ecoli_kmer_counts, total_ecoli_kmers)
  tetrasphaera_kmer_freqs <- calculate_kmer_frequencies(tetrasphaera_kmer_counts, total_tetrasphaera_kmers)

  # Merge results into one data frame
  combined_results <- merge(ecoli_kmer_freqs, tetrasphaera_kmer_freqs, by = "kmers", suffixes = c("_E_coli", "_Tetrasphaera"))
  
  # Store the results
  kmer_results[[as.character(k)]] <- combined_results
}

# Combine results for all k-values
final_results <- do.call(rbind, kmer_results)

# Step 8: Reshape for plotting
final_results_melted <- melt(final_results, id.vars = "kmers",
                              measure.vars = c("Frequency_E_coli", "Frequency_Tetrasphaera"),
                              variable.name = "Organism", value.name = "Frequency")

# Modify the Organism names for clarity
final_results_melted$Organism <- gsub("Frequency_", "", final_results_melted$Organism)

# Step 9: Plot k-mer frequency comparison
ggplot(final_results_melted, aes(x = kmers, y = Frequency, fill = Organism)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "K-mer Frequency Comparison (3-5 Amino Acids)",
       x = "K-mer",
       y = "Frequency") +
  scale_fill_manual(values = c("E_coli" = "blue", "Tetrasphaera" = "green")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


```
#Loads the libraries: seqinr for reading FASTA and translating DNA to protein, ggplot2 for plotting, dplyr for data manipulation, and reshape2 for reshaping data.
#Reads the FASTA files containing the coding sequences for E. coli and Tetrasphaera using the read.fasta function.
translate_to_protein <- function {
#This function translates the DNA sequences into protein sequences. It first uses the seqinr::translate function that will change each sequence into amino acids.
#calculate_kmer_frequencies: This function will calculate the frequency of each kmer by dividing its count by the total number of kmers in the organism. melt is used to reshape the data for plotting. It transforms the frequency columns for E. coli and Tetrasphaera into a long format with separate columns for the organism and frequency. ggplot2 creates a bar plot comparing the k-mer frequencies between E. coli and Tetrasphaera.

#As difference between the two organisms From the hypothetical data, there is a higher codon usage bias for Escherichia coli, as shown by the lower ENC and higher CAI. This may be because it is adapted to the laboratory conditions under which there is more efficient translation of proteins by adopting particular codons.With a higher ENC and lower CAI, Tetrasphaera sp. may thus reflect a more balanced usage of codons, probably because this bacterium lives under diverse environmental conditions in which codon usage is not under strong selective pressure. These differences in codon usage could then reflect variation in their gene expression strategies or in their metabolic adaptation.Overall, the above observations underpin how evolutionary pressures and ecological niches act upon codon usage and bias in both E. coli, with its more specialized codon preference augmenting its fitness in laboratory and clinical environments, and Tetrasphaera, which maintains a broader codon usage strategy better suited to its natural habitat.#
```
#Question6
#In the organism of interest, identify 10 protein sequence k-mers of length 3-5 which are the most over- and under-represented k-mers in your organism of interest. Are these k-mers also over- and under-represented in E. coli to a similar extent? Provide plots to support your observations. Why do you think these sequences are present at different levels in the genomes of these organisms?

# Identify 10 protein sequence k-mers of length 3-5 which are the most over- and under-represented k-mers in Tetrasphaera
```{r}
# Load required libraries
library(seqinr)
library(dplyr)

# Step 1: Read the Tetrasphaera protein sequences
tetrasphaera_proteins <- seqinr::read.fasta("tetrasphaera_proteins.fa")
tetrasphaera_proteins <- sapply(tetrasphaera_proteins, function(x) paste(x, collapse = ""))

# Step 2: Function to extract k-mers
extract_kmers <- function(proteins, k) {
  kmers <- unlist(lapply(proteins, function(seq) {
    sapply(1:(nchar(seq) - k + 1), function(i) substr(seq, i, i + k - 1))
  }))
  return(kmers)
}

# Step 3: Count k-mers
count_kmers <- function(kmers) {
  return(as.data.frame(table(kmers)))
}

# Step 4: Calculate expected frequencies
calculate_expected_frequency <- function(total_kmers, k) {
  amino_acid_count <- 20  # Number of amino acids
  expected_count <- (total_kmers - k + 1) / (amino_acid_count^k)
  return(expected_count)
}

# Step 5: Identify over- and under-represented k-mers
identify_represented_kmers <- function(kmer_counts, expected_count) {
  kmer_counts$Expected <- expected_count
  kmer_counts$Log2FC <- log2(kmer_counts$Freq / kmer_counts$Expected)

  # Identify over-represented and under-represented k-mers
  over_represented <- kmer_counts %>% filter(Log2FC > 1) %>% arrange(desc(Log2FC))
  under_represented <- kmer_counts %>% filter(Log2FC < -1) %>% arrange(Log2FC)

  return(list(over_represented = over_represented, under_represented = under_represented))
}

# Step 6: Perform analysis for k = 3 to 5
kmer_results <- list()

for (k in 3:5) {
  tetrasphaera_kmers <- extract_kmers(tetrasphaera_proteins, k)
  tetrasphaera_kmer_counts <- count_kmers(tetrasphaera_kmers)

  total_tetrasphaera_kmers <- sum(tetrasphaera_kmer_counts$Freq)
  expected_count <- calculate_expected_frequency(total_tetrasphaera_kmers, k)

  results <- identify_represented_kmers(tetrasphaera_kmer_counts, expected_count)
  kmer_results[[as.character(k)]] <- results
}

# Step 7: Display results
for (k in 3:5) {
  cat(paste("K =", k, "\n"))
  cat("Over-represented k-mers:\n")
  print(head(kmer_results[[as.character(k)]]$over_represented, 10))
  cat("Under-represented k-mers:\n")
  print(head(kmer_results[[as.character(k)]]$under_represented, 10))
  cat("\n")
}

```
#Comparison with the E.coli
```{r}
# Load required libraries
library(seqinr)
library(dplyr)
library(ggplot2)

# Function to read protein sequences
read_proteins <- function(file_path) {
  proteins <- seqinr::read.fasta(file_path)
  return(sapply(proteins, function(x) paste(x, collapse = "")))
}

# Step 1: Read the protein sequences for E. coli and Tetrasphaera
tetrasphaera_proteins <- read_proteins("tetrasphaera_proteins.fa")
ecoli_proteins <- read_proteins("ecoli_proteins.fa")

# Step 2: Function to extract k-mers
extract_kmers <- function(proteins, k) {
  kmers <- unlist(lapply(proteins, function(seq) {
    sapply(1:(nchar(seq) - k + 1), function(i) substr(seq, i, i + k - 1))
  }))
  return(kmers)
}

# Step 3: Count k-mers
count_kmers

#This part of the script reads the Tetrasphaera protein sequences from a FASTA file and converts the entries into strings using read.fasta and apply, then uses extract_kmers, and sets them as a data frame.
#expected_frequency calculates the expected frequency of each k-mer under a random distribution (20 possible amino acids, to the power of k).
#Loops over k = 3,4,5; extracts k-mers; counts them; calculates expected frequencies and over- and under-represented k-mers. Results are returned as a list.
#Reads in the E. coli and Tetrasphaera protein sequences using read_proteins.
extract_kmers function extracts the k-mers to be compared from both organisms.
#This code is incomplete for the E. coli comparison, but presumably you would apply the same count_kmers and calculate_expected_frequency functions as done for Tetrasphaera, followed by an additional comparison step, for example using ggplot2 for a visual comparison.

#As differences between two organisms  in k-mer frequencies between the organism of interest and Escherichia coli bear immense functional and evolutionary implications. The over-represented k-mers in the organism of interest may represent important functional or structural motifs in adapting to specific roles in biology, such as metabolic pathways, stress responses, and adaptations to ecological niches. On the other hand, the under-represented k-mers may indicate the lack of structural features or functions compared with E. coli and may further suggest evolutionary losses or shifts in functional relevance. These differences in k-mer representation are the result of different evolutionary pressures the two organisms have faced; E. Whereas coli is a model organism that has been extensively studied and optimized for laboratory conditions, the organism of interest may have evolved to prosper in unique environmental niches, thus affecting its protein composition. In general, the presence and frequency of some k-mers carry meaningful information on the evolutionary history, ecological adaptations, and functional capabilities of such organisms. The list of such k-mers could be further analyzed; for example, functional annotations could be done for proteins that may use them, toward more certain conclusions about biological significance.

```
