library(readxl)
library(stringr)
library(tidyverse)
library(VennDiagram)
splitAndCombine <- function(gene, rxn,sep0) { ##one rxn has several genes, this function was used to splite the genes
  
  gene <- str_split(gene, sep0)
  tt<- length(gene)
  gene0 <- list()
  for (i in 1:tt){
    gene0[[i]] <- paste(rxn[i], gene[[i]], sep = "@@@")
    
  }
  
  gene1 <- unique(unlist(gene0))
  gene2 <- str_split(gene1, "@@@" )
  rxnGene <- data.frame(v1=character(length(gene2)),stringsAsFactors = FALSE)
  tt1 <- length(gene2)
  for (j in 1:tt1){
    rxnGene$v1[j] <- gene2[[j]][2]
    rxnGene$v2[j] <- gene2[[j]][1]
  }
  
  return(rxnGene)
}

#panID from TCDB database
newTCDB <- read_excel("Blast results/new transport rxn from TCDB.xlsx", 
                      sheet = "new transport rxn from TCDB")
newTCDB$panID <- str_trim(newTCDB$panID, side = "both")

#panID from bigg blast
newBigg <- read_excel("BiGG/new rxn from pan.xlsx")
newBigg$panID <- str_trim(newBigg$panID, side = "both")

#panID from uniprot blast
newUniprot <- read_excel("Blast results/new ec from pan based on uniprot.xlsx")
newUniprot$panID <- str_trim(newUniprot$panID, side = "both")

#panID from raven blast
newRaven <- read_excel("model files for panGenome and s288c/newRxnCombine.xlsx", 
                       sheet = "newRxnCombine")

newRaven0 <- splitAndCombine (newRaven$`GENE ASSOCIATION`,newRaven$ID,sep0 = "or") 
colnames(newRaven0) <- c('panID','KEGG or Biocyc')
newRaven0$panID <- str_trim(newRaven0$panID, side = "both")


index1 <- which(newRaven0$panID %in% newBigg$panID ==TRUE)
index2 <- which(newRaven0$panID %in% newUniprot$panID ==TRUE)
unique(newRaven0$panID[index2])


#summary
BIGG <- unique(newBigg$panID)
RAVEN <- unique(newRaven0$panID)
Uniprot <- unique(newUniprot$panID)
TCDB <- unique(newTCDB$panID)
#plot the graph
venn.diagram(x= list(BIGG = BIGG, RAVEN = RAVEN, Uniprot = Uniprot, TCDB = TCDB), 
             filename = "new gene in main database.png", height = 1000, width = 1000,resolution =300, imagetype="png", col="transparent",
             fill=c("cornflowerblue","green","darkgreen","darkorchid1"),alpha = 0.50, cex=0.45, cat.cex=0.45)


venn.diagram(x= list(BIGG = BIGG, RAVEN = RAVEN, Uniprot = Uniprot), 
             filename = "new gene in main database.png", height = 1000, width = 1000,resolution =300, imagetype="png", col="transparent",
             fill=c("cornflowerblue","green","darkgreen"),alpha = 0.50, cex=0.45, cat.cex=0.45)









