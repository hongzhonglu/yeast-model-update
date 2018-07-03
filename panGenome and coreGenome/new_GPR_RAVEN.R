library(tidyverse)
library(stringr)
library(readxl)

getSingleReactionFormula <- function(description, reaction_ko, ko) {###description can be any charater of metabolite
  index <- vector()
  result <- vector()
  tt <- vector()
  for (i in 1:length(ko)){
    if(length(match(ko[i],reaction_ko))){
      index <- match(ko[i],reaction_ko)
      tt <- description[index]
      result[i] <- paste0(tt, collapse = ";")
    } else{
      
      result[i] <- NA
    }
  }
  return(result)
}
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

panModel_kegg <- read_excel("model files for panGenome and s288c/panModel_kegg.xlsx")
panModel_metacyc <- read_excel("model files for panGenome and s288c/panModel_metacyc.xlsx")

sceModel_kegg <- read_excel("model files for panGenome and s288c/sceModel_kegg.xlsx")
sceModel_metacyc <- read_excel("model files for panGenome and s288c/sceModel_metacyc.xlsx")


# new rxn

newRxn_kegg <- setdiff(panModel_kegg$ID, sceModel_kegg$ID)
newRxn_metacyc <- setdiff(panModel_metacyc$ID, sceModel_metacyc$ID)

index1 <- which(panModel_kegg$ID %in% newRxn_kegg ==TRUE)
index2 <- which(panModel_metacyc$ID %in% newRxn_metacyc ==TRUE)


newRxn_kegg0 <- panModel_kegg [index1,]
newRxn_kegg0$source <- "KEGG"
newRxn_metacyc0 <- panModel_metacyc[index2,]
newRxn_metacyc0$source <- "metacyc"

newRxnCombine <- rbind.data.frame(newRxn_kegg0,newRxn_metacyc0)

write.table(newRxnCombine,"newRxnCombine.txt", row.names = FALSE, sep = "\t")
