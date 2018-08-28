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

caselist <- getCaseLists(mycgds, blca)

caselist <- caselist[1,1]

geneticprofiles <- getGeneticProfiles(mycgds,blca)

geneticprofiles <- geneticprofiles$genetic_profile_id

queriedgenes <- fread("/home/nikolas/stravopodis/Bladder Cancer/queriedgenes.txt",
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
colSums(rnahighzscore)

rownames(rnahighzscore)[rnahighzscore[,1] == 1]

###
#TCGAbiolinks part
###

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

length <- mutationmatrix$Tumor_Sample_Barcode
length <- lapply(length, length)
length <- (unlist(length))
mutationmatrix <- mutationmatrix[length > mutatednumber,]

mutationmatrix$Tumor_Sample_Barcode <- lapply(mutationmatrix$Tumor_Sample_Barcode,substring,1,16)

rna <- fread("/home/nikolas/stravopodis/Bladder Cancer/TCGA-BLCA_mRNA_FPKM_countmatrix.tsv")

rna <- data.frame(rna)

rownames <- rna$V1

rna <- rna[,2:ncol(rna)]

normals <- substring(colnames(rna),14,15) != 11

rna <- rna[,normals]

rownames(rna) <- rownames

rownames(rna) <- substring(rownames(rna),0,15)

colnames(rna) <- substring(colnames(rna),1,16)

colnames(rna) <- make.unique(colnames(rna))

genes <- rownames(rna)

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

G_list <- getBM(filters= "ensembl_gene_id",
                attributes= c("ensembl_gene_id","hgnc_symbol"),
                values=genes,
                mart= mart)

queriedgenes <- subset(G_list,
                       hgnc_symbol %in% queriedgenes)

###TESTO
n <- 1
gene <- mutationmatrix[n,1]
samples <- make.names(unlist(mutationmatrix[n,2]))

j <- 1
queriedgene <- queriedgenes[j,1]
queriedgene
duplicated(samples)

data <- rna[queriedgene,]

dim(data)

mutated <- as.numeric(data[samples])

setdiff(colnames(data),samples)

nonmutated <- data[, -which(colnames(data) %in% samples)]

union(samples,colnames(data))

nonmutated <- as.numeric(data[nonmutated])
colnames(data)[duplicated(colnames(data))]
