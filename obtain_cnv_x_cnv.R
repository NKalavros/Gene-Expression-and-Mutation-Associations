#Load packages
library(TCGAbiolinks)
library(data.table)

#Function variable
project <- "TCGA-BLCA"
#Function variable
queriedgenes <- fread("file:///C:/Users/lab_ema/Desktop/Nikolas/Nikolas Stravopodis/queriedgenes.txt",
                      header = FALSE)
queriedgenes <- queriedgenes$V1

cnv <- getGistic(substring(project,6,9), type = "thresholded")

cnv <- cnv[,c(1,4:ncol(cnv))]

cnv[cnv==-2] <- -1
cnv[cnv==+2] <- +1

#From 50 queried genes, only 47 had CNVs
queriedgenescnv <- subset(cnv,
                          cnv$`Gene Symbol` %in% queriedgenes)

#Getting the rest of the genes
restgenescnv <- subset(cnv,
                       !(cnv$`Gene Symbol` %in% queriedgenes))

#Getting the rownames correct
rownames(queriedgenescnv) <- queriedgenescnv$`Gene Symbol`

#Creating the function to aggregate the names
foo <- function(dfrow){
  gain <- names(dfrow)[dfrow == "1"]
  loss <- names(dfrow)[dfrow == "-1"]
  result <- list(gain,loss)
  names(result) <- c("Gain","Loss")
  return(result)
}

#Aggregating gene names gain loss into a list for queried and the rest
queriedgenesaggregatedcnv <- (apply(queriedgenescnv,1,foo))

restgenesaggregatedcnv <- apply(restgenescnv,1,foo)

###Beginning creating the gain x gain


#Creating the necessary matrices
direction <- matrix(0, nrow = nrow(queriedgenescnv), ncol = nrow(restgenescnv))
pvalues <- rep(0.0,nrow(queriedgenescnv)*ncol(restgenescnv))

count <- 0
allnames <- colnames(queriedgenescnv)

for(i in seq(1:length(queriedgenesaggregatedcnv))){
  gainqueried <- queriedgenesaggregatedcnv[[i]]$Gain
  for(j in seq(1:length(restgenesaggregatedcnv))){
    gainrest <- restgenesaggregatedcnv[[j]]$Gain
    both <- length(intersect(gainqueried,gainrest))
    onlyqueried <- length(setdiff(gainqueried,gainrest))
    onlyrest <- length(setdiff(gainrest,gainqueried))
    noneofthetwo <- length(setdiff(setdiff(allnames,gainqueried),gainrest))
    testtable <- matrix(data = c(both,onlyrest,onlyqueried,noneofthetwo), nrow = 2, ncol =2, byrow=TRUE)
    test <- fisher.test(testtable)
    count <- count + 1
    pvalues[count] <- test$p.value
    if((both >= onlyqueried) & (both >= onlyrest)){
      direction[i,j] <- 1
    }
    print(count)
  }
}

