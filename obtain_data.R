library(cgdsr)
library(data.table)
library(TCGAbiolinks)
library(biomaRt)
###
#cBioportal part
###

mycgds <- CGDS("http://www.cbioportal.org/")

test(mycgds)  

cancerstudies <- getCancerStudies(mycgds)

blca <- cancerstudies[16,1]

#Add flags here
print(blca)

caselist <- getCaseLists(mycgds, blca)

caselist <- caselist[1,1]

#Add flags here too
print(caselist)

geneticprofiles <- getGeneticProfiles(mycgds,blca)

geneticprofiles <- geneticprofiles$genetic_profile_id

#Function variable
queriedgenes <- fread("file:///C:/Users/lab_ema/Desktop/Nikolas/Nikolas Stravopodis/queriedgenes.txt",
                      header = FALSE)

queriedgenes <- queriedgenes$V1

gistic <- getProfileData(mycgds, genes = c(queriedgenes), geneticProfiles = geneticprofiles[3],
                         caseList = c(caselist))

rnazscore <- getProfileData(mycgds, genes = c(queriedgenes), geneticProfiles = geneticprofiles[5],
                            caseList = c(caselist))

mutations <- getProfileData(mycgds, genes = c(queriedgenes), geneticProfiles = geneticprofiles[8],
                            caseList = c(caselist))

methylation <- getProfileData(mycgds, genes = c(queriedgenes), geneticProfiles = geneticprofiles[7],
                              caseList = c(caselist))

rnahighzscore <- matrix(data = 0, nrow = nrow(rnazscore),
                        ncol = ncol(rnazscore))

rownames(rnahighzscore) <- rownames(rnazscore)
colnames(rnahighzscore) <- colnames(rnazscore)

#Function variable

threshold <- 1.5

for(i in seq(1:nrow(rnazscore))){
  for(j in seq(1:ncol(rnazscore))){
    if((rnazscore[i,j] > threshold)){
      rnahighzscore[i,j] <- 1
    }
    else if(rnazscore[i,j] < -threshold){
      rnahighzscore[i,j] <- -1
    }
  }
}

###
#TCGAbiolinks part
###

#Function variable
project <- "TCGA-BLCA"

#13000
GenomicDataCommons::status()

maf <- GDCquery_Maf(substring(project,6,9), pipelines = "muse")

#90000
maf <- subset(maf,
              Variant_Classification != "Silent")

maf <- maf[,c(1,16)]

mutationmatrix <- aggregate(. ~ Hugo_Symbol,
                            data = maf,
                            FUN = paste)
#Keeping only those mutated in at least 10% of samples

#Function variable
mutatednumber <- 40

length <- mutationmatrix$Tumor_Sample_Barcode
length <- lapply(length, length)
length <- (unlist(length))
mutationmatrix <- mutationmatrix[length > mutatednumber,]

mutationmatrix$Tumor_Sample_Barcode <- lapply(mutationmatrix$Tumor_Sample_Barcode,substring,1,15)
mutationmatrix$Tumor_Sample_Barcode <- lapply(mutationmatrix$Tumor_Sample_Barcode,make.names)
barcodes <- mutationmatrix$Tumor_Sample_Barcode

rna <- fread("file:///C:/Users/lab_ema/Desktop/Nikolas/Nikolas Stravopodis/TCGA-BLCA_mRNA_FPKM_countmatrix.tsv")

rna <- data.frame(rna)

rownames <- rna$V1

rna <- rna[,2:ncol(rna)]

normals <- substring(colnames(rna),14,15) != 11

rna <- rna[,normals]

rownames(rna) <- rownames

rownames(rna) <- substring(rownames(rna),0,15)

colnames(rna) <- substring(colnames(rna),1,15)

colnames(rna) <- make.unique(colnames(rna))

foo <- function(character){
  return(character[character %in% colnames(rna)])
}

mutationmatrix$Tumor_Sample_Barcode <- lapply(mutationmatrix$Tumor_Sample_Barcode,foo)

genes <- rownames(rna)

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

G_list <- getBM(filters= "ensembl_gene_id",
                attributes= c("ensembl_gene_id","hgnc_symbol"),
                values=genes,
                mart= mart)

queriedgenes <- subset(G_list,
                       hgnc_symbol %in% queriedgenes)

mutationxexpression <- matrix(0,nrow = nrow(mutationmatrix),
                              ncol = nrow(queriedgenes))

pvalues <- matrix(0.0,nrow = nrow(mutationmatrix),
                  ncol = nrow(queriedgenes))

pvaluescorrected <- rep(0.0,nrow(mutationmatrix)*nrow(queriedgenes))

count <- 0
###TESTO
for(n in seq(1:length(mutationmatrix$Hugo_Symbol))){
  j <- 1
  gene <- mutationmatrix[n,1]
  samples <- make.names(unlist(mutationmatrix[n,2]))
  duplicatedsamples <- duplicated(samples)
  samples <- samples[!duplicatedsamples]
  for(j in seq(1:length(queriedgenes$hgnc_symbol))){
    queriedgene <- queriedgenes[j,1]
    queriedgene

    data <- rna[queriedgene,]
    mutated <- as.numeric(data[samples])

    nonmutated <- as.numeric(data[, -which(colnames(data) %in% samples)])

    ttest <- t.test(mutated,nonmutated)
    if(ttest$estimate[1] > ttest$estimate[2]){
      mutationxexpression[n,j] <- 1
    }
    count <- count + 1
    pvaluescorrected[count] <- ttest$p.value
    print(count)
  }
  
}
pvaluescorrected <- p.adjust(pvaluescorrected, method = "fdr")
sum(pvaluescorrected < 0.05)

pvalues <- matrix(pvaluescorrected,nrow = nrow(mutationmatrix),
                  ncol = nrow(queriedgenes),byrow = TRUE)

rownames(mutationxexpression) <- rownames(pvalues) <- mutationmatrix$Hugo_Symbol
colnames(mutationxexpression) <- colnames(pvalues) <- queriedgenes$hgnc_symbol

write.table(mutationxexpression, "mutationxexpression.tsv", sep = "\t")
write.table(pvalues, "mutationxexpressionpvaluescorrected.tsv", sep = "\t")
