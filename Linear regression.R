#Linear regression Analysis

#Check dimension of untreated_normalized and copynumber to see if they are the same
dim(untreated_normalized)
dim(copynumber)

#Dimensions are completely different:
#1. Copynumber lists more genes than untreated_normalized
# Solution: 
#  - Check if rows (genes) are ordered alphabetically
#  - If not: order both matrixes alphabetically by row name
#  - Delete rows that appear in only one of the matrixes

#2. untreated_normalized has more columns, since there is one column for each drug
# Solution:
# - New matrix with the mean for each cellline in untreated_normalized
# - Order columns alphabetically in both matrixes if not already ordered

# ------- 1. ------- 

#COPY NUMBER:

identical(rownames(copynumber), rownames(copynumber[order(rownames(copynumber)),]))
# False, so order them:
cnforlr <- copynumber[order(rownames(copynumber)),]

commongenes <- intersect(rownames(cnforlr), rownames(utnforlr))
cnforlr <- cnforlr[which(rownames(cnforlr) %in% commongenes),]

#UNTREATED:

identical(rownames(untreated_normalized), rownames(untreated_normalized[order(rownames(untreated_normalized)),]))
# False, so order them:
utnforlr <- untreated_normalized[order(rownames(untreated_normalized)),]

utnforlr <- utnforlr[which(rownames(utnforlr) %in% commongenes),]

# ----------- 2. ----------- 
