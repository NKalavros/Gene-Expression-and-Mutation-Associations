library(cgdsr)
library(data.table)
library(TCGAbiolinks)

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

queriedgenes <- fread("file:///C:/Users/lab_ema/Desktop/Nikolas Stravopodis/queriedgenes.txt",
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

###
#TCGAbiolinks part
###
project <- "TCGA-BLCA"

query <- GDCquery(project = project,
                  data.category = "Gene Expression",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq- Counts")


