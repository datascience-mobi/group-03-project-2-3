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

#Cell lines in matrix utnforlr are untreated, so the part of the colnames stating drug, time and dose (0 nM) are unnecessary
#Delete that part of the name to leave only the cellline

colnames(utnforlr) <- sub("_.*", "", colnames(utnforlr))

#order columns alphabetically 
utnforlr <- utnforlr[, order(colnames(utnforlr))]  
dim(utnforlr)
colnames(utnforlr[,1:15])

#There are  many columns for each cell line, so they need to be merged and have their mean calculated
utnforlr <- sapply(split(seq_len(ncol(utnforlr)),colnames(utnforlr)),function(cis) rowMeans(utnforlr[,cis,drop=F]))
dim(utnforlr)
colnames(utnforlr[,1:15])
#Number of columns has been reduced and now there is only one for each cellline

#Order columns alphabetically in the copy number matrix
cnforlr <- copynumber[, order(colnames(copynumber))]

#Keep in both matrixes only the columns (celllines) that are present in both 
commoncelllines <- intersect(colnames(cnforlr), colnames(utnforlr))
cnforlr <- cnforlr[, which(colnames(cnforlr) %in% commoncelllines)]
utnforlr <- utnforlr[, which(colnames(utnforlr) %in% commoncelllines)]

#Now both matrixes have the same rows (genes) and columns(celllines)

dsforlr <- drugsensitivity[which(rownames(drugsensitivity) %in% commoncelllines),]
utnbmforlr <- utnforlr[which(rownames(utnforlr) %in% biomarkers),]

linearregression <- rbind(utnbmforlr, dsforlr[, "gemcitibine"])
linearregression <- t(linearregression)
namesforlr <- biomarkers[order(biomarkers)]
namesforlr <- namesforlr[-which(namesforlr=="CXCL8")]
namesforlr <- c(namesforlr, "DS_gemcitibine")
colnames(linearregression) <- namesforlr
                      
#Plot each biomarker vs the drug sensitivity

DS_gemcitibine <- linearregression[, "DS_gemcitibine"]
for (i in c(1:14)) {
  GeneX <- linearregression[, i]
  name <- namesforlr[i]
  plot(GeneX, DS_gemcitibine, xlab = namesforlr[i], ylab = "Drug sensitivity")
  abline(lm(DS_gemcitibine~linearregression[,i], data = as.data.frame(linearregression)), col = "red")
}
