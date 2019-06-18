#loading library limma needed for normalization
library(limma)

#normalize treated & untreated
untreated_normalized <- normalizeBetweenArrays(NCI_TPW_gep_untreated)
treated_normalized <- normalizeBetweenArrays(NCI_TPW_gep_treated)

#Change names of (not normalized) trated and untreated data sets
treated <- NCI_TPW_gep_treated
untreated <- NCI_TPW_gep_untreated

#-------BIOMARKERS-----------

#Fold-change: values are already log2 transformed, so only substract the two matrices
fcgeneexpression <- treated_normalized - untreated_normalized


#Absolute value of the data in fcgeneexpression, to find the biggest changes, regardless of if they are up- or down-regulated
absvalfcgeneexpression <- abs(fcgeneexpression)

#New matrices with only the columns about gemcitabine:
library("dplyr")

treated_gemcitabine <- select(as_tibble(treated_normalized), contains("gemcitibine"))
treated_gemcitabine <- as.matrix(treated_gemcitabine)
rownames(treated_gemcitabine) <- rownames(treated_normalized)

untreated_gemcitabine <- select(as_tibble(untreated_normalized), contains("gemcitibine"))
untreated_gemcitabine <- as.matrix(untreated_gemcitabine)
rownames(untreated_gemcitabine) <- rownames(untreated_normalized)

absvalfcge_gemcitabine <- select(as_tibble(absvalfcgeneexpression), contains("gemcitibine"))
absvalfcge_gemcitabine <- as.matrix(absvalfcge_gemcitabine)
rownames(absvalfcge_gemcitabine) <- rownames(absvalfcgeneexpression)


#Create new function to get a list of the genes with the biggest difference in expression after treatment for a single column (cell line) of absvalfcge_gembitabine:
orderforbiomarkers <- function(absvalfcge_gemcitabine, n) {
  gfc <- absvalfcge_gemcitabine[ , n]
  orderedgenes <- names(gfc[order(-gfc)])
  return(orderedgenes)
}

#Function to order the components of a vector (descending) and return the names in that order
orderforvectors <- function(x) {
  y <- names(x[order(-x)])
  return(y)
}

#Apply the new funtion orderforvectors to all the columns and create a new matrix with the gene names
matrixofbiomarkers <- apply(absvalfcge_gemcitabine, 2, orderforvectors) #the 2 means it is applied to columns


#Check for most common genes in the first 15 rows of the new matrix -> biomarkers!
topgenes <- sort(table(as.vector(matrixofbiomarkers[1:15,])))

#Convert topgenes into a matrix in descending order of counts
matrixtopgenes <- as.matrix(topgenes[order(-topgenes)])

#Bar plot of topgenes to find the biomarkers
barplot(matrixtopgenes[1:20,], xlab = "Genes", ylab = "Counts", main = "Top differentially expressed genes")
#We decided to take the top 15 genes as our biomarkers

#List of our biomarkers
biomarkers <- names(matrixtopgenes[1:15,1])

