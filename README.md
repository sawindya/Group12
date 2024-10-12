# Group12 
#TIP: To download a file with R, click on “view raw” and then you can copy the URL from the address bar
and then use the download.file command in R.
url <- "https://raw.githubusercontent.com/ghazkha/Assessment4/main/gene_expression.tsv"
destfile <- "gene_expression.tsv"
download.file(url, destfile)
url <- "https://raw.githubusercontent.com/ghazkha/Assessment4/main/growth_data.csv"
destfile <- "growth_data.csv"
download.file(url, destfile)
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
