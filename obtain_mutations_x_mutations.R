#Loading packages
library(TCGAbiolinks)
library(data.table)


#Function variable
project <- "TCGA-BLCA"
#Function variable
queriedgenes <- fread("file:///C:/Users/lab_ema/Desktop/Nikolas/Nikolas Stravopodis/queriedgenes.txt",
                      header = FALSE)
queriedgenes <- queriedgenes$V1

#Query Maf, copied code snippet from previous script

maf <- GDCquery_Maf(substring(project,6,9), pipelines = "muse")

#119159 at first, deleting silents
maf <- subset(maf,
              Variant_Classification != "Silent")
#92532 now

#Keeping only the columns of hugo symbol and sample barcode
maf <- maf[,c(1,16)]

#Aggregating all barcodes by symbol
mutationmatrix <- aggregate(. ~ Hugo_Symbol,
                            data = maf,
                            FUN = paste)
#Keeping only those mutated in at least 10% of samples

#Function variable
#Tried with 5 and 6 7 and 8 and 9, to no avail
mutatednumber <- 10

#Getting the whole list of character vectors with each barcode
length <- mutationmatrix$Tumor_Sample_Barcode

#Finding the length of each character vector in the list using lapply length
length <- lapply(length, length)

#Turning the list into a numeric vector
length <- (unlist(length))

#Keeping only the ones who have more samples mutated than the function variable mutatednumber

#54 genes have been left
mutationmatrix <- mutationmatrix[length > mutatednumber,]

#Substringing to keep only the unique identifiers
mutationmatrix$Tumor_Sample_Barcode <- lapply(mutationmatrix$Tumor_Sample_Barcode,substring,1,15)

#Making names so as to not have problems with hyphens (-)
mutationmatrix$Tumor_Sample_Barcode <- lapply(mutationmatrix$Tumor_Sample_Barcode,make.names)

#Keeping only unique identifiers
mutationmatrix$Tumor_Sample_Barcode <- lapply(mutationmatrix$Tumor_Sample_Barcode,unique)

#Extracting the mutated people of the genes in question, first, assessing which genes are mutated
bool <- queriedgenes %in% mutationmatrix$Hugo_Symbol

#Used mutatednumber = 5
queriedgenes <- queriedgenes[bool]

#Getting mutations in our queried genes
queriedgenesmutations <- subset(mutationmatrix,
                                mutationmatrix$Hugo_Symbol %in% queriedgenes)

#Getting the rest of the mutations
restofmutations <- subset(mutationmatrix,
                          !(mutationmatrix$Hugo_Symbol %in% queriedgenes))

#Creating the count variable

count <- 0
#Getting all the names for the fisher testing
names <- mutationmatrix$Tumor_Sample_Barcode
names <- unique(unlist(names))

#Creating the two empty matrices
direction <- matrix(0,nrow = nrow(queriedgenesmutations),
                    ncol = nrow(restofmutations))
pvalues <- rep(0.0,nrow(queriedgenesmutations)*nrow(restofmutations))

#Looping over all the names and performing the fisher tests
for(i in seq(1:nrow(queriedgenesmutations))){
  mutatedbarcodes <- unlist(queriedgenesmutations$Tumor_Sample_Barcode[i])
  for(j in seq(1:nrow(restofmutations))){
    restmutationsbarcodes <- unlist(restofmutations$Tumor_Sample_Barcode[j])
    both <- length(intersect(restmutationsbarcodes,mutatedbarcodes))
    onlyqueriedgene <- length(setdiff(mutatedbarcodes,restmutationsbarcodes))
    onlyrest <- length(setdiff(restmutationsbarcodes,mutatedbarcodes))
    noneofthetwo <- length(setdiff(setdiff(names,mutatedbarcodes),restmutationsbarcodes))
    testtable <- matrix(data = c(both,onlyrest,onlyqueriedgene,noneofthetwo), nrow = 2, ncol =2, byrow=TRUE)
    test <- fisher.test(testtable)
    count <- count + 1
    pvalues[count] <- test$p.value
    if((both >= onlyqueriedgene) & (both >= onlyrest)){
      direction[i,j] <- 1
    }
    print(count)
  }
}

#Correcting p values and turning them into a matrix for saving
pvaluescorrected <- p.adjust(pvalues,method = "fdr")
sum(pvaluescorrected < 0.05)
pvaluescorrected <- matrix(data = pvaluescorrected, nrow = nrow(queriedgenesmutations),
                           ncol = nrow(restofmutations), byrow = TRUE)

#Creating row and colnames for easier identification
rownames(pvaluescorrected) <- rownames(direction) <- queriedgenesmutations$Hugo_Symbol
colnames(pvaluescorrected) <- colnames(direction) <- restofmutations$Hugo_Symbol

write.table(direction, "mutationxmutation.tsv", sep = "\t")
write.table(pvaluescorrected, "mutationxmutationpvaluescorrect.tsv", sep = "\t")
