# Dataset 1
# 1) Use PCA to reduce dimensions.  How many components do you need to keep
# to reproduce the digits reasonably well? what is your final matrix ?
data <- read.csv("/Users/sarahferaidoon/Desktop/train.csv.gz")
data

# Number 2
n <- c(0:9)
# keep track of all digits corresponding to a handwriting from
# 0 to 9 in variable n
for (x in n){
  # loop through each digit
  # assign puts all data we want to grep in a variable
  # paste0 creates variable names dynamically
  # grep takes digit from first column following pattern
  # put entire row matching digit into a variable
  # loop moves onto next digit
  assign(paste0("label_",x),data[grep(x,data[1]),])
}
label_0
label_1
label_2
label_3
label_4
label_5
label_6
label_7
label_8
label_9
data[grep("0",data[,1]),]


length(label_1)
# scale FALSE to take out columns with variance of 0
pca_label_1 <- prcomp(label_1[,c(2:784)],center=TRUE,scale=FALSE)
summary(pca_label_1)
# St dev goes below 0.05 after PC496. It is above 0.05 from PC1-PC496
plot(pca_label_1$x[,1],pca_label_1$x[,2])
biplot(pca_label_1, scale=0)
# using screeplot, find the pc cut off that contains 90% of the variance
# and convert it back to the data to see how the digits look
screeplot(pca_label_1,type="lines")
new_label_1 <- label_1[,c(2:496)]

# get standard deviation
standard_dev <- pca_label_1$sdev
# get variance of pc
var_label_1 <- standard_dev^2
var_label_1[1:10]

# explain proportion of variance
prop_var <- var_label_1/sum(var_label_1)
prop_var[1:35]
# the first 35 pc (principal components) consist of 90% of the variance
sum(prop_var[1:35])
#convert the first 35 pc back into data and see how digit looks
pc.use <- 35
???reconstructed <- pca_label_1$x[,1:pc.use] %*% t(pca_label_1$rotation[,1:pc.use])
reconstructed
image(reconstructed, axes=FALSE, col=grey.colors(256, start=0, end=1))
pca_label_1$x[,1:pc.use]


# Dataset 2
# 1) Build hierarchical trees based on the columns and for
# the rows (exclude rows that are "low" expression)
# read data file as variable
countdata <- read.csv(file="/Users/sarahferaidoon/Desktop/Mnemiopsis_count_data.csv", header=TRUE)
# read data file as dataframe
countdata <- as.data.frame(countdata)
countdata
# use grep function to extract all cases of "aboral"
print(grep("aboral",countdata))
# use grep function to extract all cases of "oral"
print(grep("oral",countdata))
# create dataframe for aboral
dfaboral <- countdata[,c(2,3,4,5)]
# create dataframe for oral
dforal <- countdata[,c(6,7,8,9)]
dfcounts <- countdata[,-1]
rownames(dfcounts) <- countdata[,1]
print(dfcounts)
# Create a dataframe called expgroup, with 1 column labeled condition (assigning healthy column). column names of readcount as expgroup's rownames
# Create a df with 8 values (corresponding to the columns in readcount)
expgroup <- data.frame(condition=1:8)
# use grep info and assign the first 20 values as healthy
expgroup$condition[c(1:4)] <- "Aboral"
# use grep info and assign the last 20 values as CF
expgroup$condition[c(5:8)] <- "Oral"
# use column names of read count as row names of expgroup
rownames(expgroup) <- colnames(countdata[,c(2:9)])
# print the df
print(expgroup)

# Step 3
# Load SeSeq2 library
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
library("DESeq2")
cds <- DESeqDataSetFromMatrix(countData = dfcounts, colData=expgroup, design = ~condition)
cds

#Step 4
# Make sure DESeq2 can correct for library size and dispersion estimates
# Use functions estimateSizeFactors and estimateDispersions
# Plot the dispersion using plotDispEsts
# Analyze graph
# Function for estimating size factors
cds <- estimateSizeFactors(cds)
# Function for estimating dispersion
cds <- estimateDispersions((cds))
# Function for plotting dispersion
plotDispEsts(cds)

# Step 5
# Perform differential expression with DeSeq and results functions
# perform differential gene expression
cds <- DESeq(cds)
# Obtain differential gene expression results
res <- results(cds)
# Show output
res

# Step 6
# Find p-values less than 0.05 and log2FoldChange > 1 or less than -1
# How many? save in list: diffexpgenes
# Genes with adjust p-value <0.05 and greater than 1
ressigind <- res[which(res$padj < 0.05 & res$log2FoldChange > 1),]
# Genes with adjust p-value <0.05 and less than -1
ressigrep <- res[which(res$padj < 0.05 & res$log2FoldChange < -1),]
# Combine results
diffexpgenes <- rbind(ressigind, ressigrep)
# print list of differentially expressed genes
rownames(diffexpgenes)
nrow(diffexpgenes)
# ~ 1500 genes

# Step 7
# Normalize values of counts data in cds with counts() function and normalized=Tm call this normvalues
normvalues <- counts(cds, normalized=TRUE)
# print normalized values
print(normvalues)
#Show columns names for matching
colnames(normvalues)

# Step 8
# Create new matrix/dataframe containing expression values from normvalues for only diffexpgenes, call this diffexpvalues
diffexpvalues <- normvalues[rownames(normvalues) %in% rownames(diffexpgenes),]
# print the dataframe
print(diffexpvalues)
# see if this matches the list of differentially expressed genes from last step
print(dim(diffexpvalues))

# Step 9
# Create hierarchical cluster
hc <- hclust(dist(diffexpgenes), method="complete")
plot <- plot(hc, main="Dendrogram", horiz=TRUE, cex=0.5)
print(plot)

# Step 10
# Create a heat map using pheatmap package
install.packages("pheatmap")
library(pheatmap)
# Plot the heatmap using designated variables
pheatmap(diffexpvalues, 
         color=colorRampPalette(c("navy","white","firebrick3"))(50), 
         scale="row", 
         cluster_rows=TRUE, 
         cellwidth=5, 
         cellheight=3, 
         fontsize=5)

# Step 11
# Use GOStats package to determine which GO-terms are enriched in diffexpgenes
# install the following packages from Bioconductor:
BiocManager::install("GOstats")
BiocManager::install("GO.db")
BiocManager::install("Category")
BiocManager::install("org.Hs.eg.db")
if (!require("BiocManager", quietly = TRUE))
  install.packages("Category")
# Load libraries
library(GOstats)
library(GO.db)
library(Category)
library(org.Hs.eg.db)
# Create a GOHyperParams object using variables created from steps
params <- new("GOHyperGParams",
              geneIds=rownames(diffexpgenes),
              universeGeneIds=rownames(diffexpvalues),
              annotation="org.Hs.eg",
              ontology="BP",
              pvalueCutoff=0.001,
              conditional=TRUE,
              testDirection="over")
# Run HyperFTest
(overRepresented-hyperGTest(params))
# Print summary of desired columns
summary(overRepresented)[,c(1,2,5,6,7)]
