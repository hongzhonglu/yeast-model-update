library(stringr)
library(tidyverse)


getMultipleReactionFormula <- function(description, reaction_ko, ko) {###description can be any charater of metabolite
  index <- vector()
  result <- vector()
  tt <- vector()
  for (i in 1:length(ko)){
    if(length( which (reaction_ko %in%  ko[i]))){
      index <- which (reaction_ko %in%  ko[i])
      tt <- description[index]
      result[i] <- paste0(tt, collapse = ";")
    } else{
      
      result[i] <- NA
    }
  }
  return(result)
}

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

splitAndCombine <- function(gene, rxn) { ##one rxn has several genes, this function was used to splite the genes
  
  gene <- str_split(gene,",")
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


###function estimate whether the gene exist in presnet model

geneExist <- function (original, newgenelist){  ##input: original is from yeast model v7.7; 
  index <- vector()
  for (i in 1: length(newgenelist)){
    if (length(which(newgenelist[i] %in% original))){
      index[i] <- "YES"
    } else {
      index[i] <- "NO"
    }
  }
  return(index)
}



#small task-reaction summary from draft yeast map
library(readxl)
rxn_cellDesigner <- read_excel("rxn_cellDesigner.xlsx", 
                               sheet = "rxn_cellDesigner")

met_cellDesigner <- read_excel("rxn_cellDesigner.xlsx", 
                               sheet = "met_cellDesigner")

#reactants split, replace and combine
rxn_reactant <- splitAndCombine(rxn_cellDesigner$reactants,rxn_cellDesigner$id)
rxn_reactant$v1_name <- getSingleReactionFormula(met_cellDesigner$name,met_cellDesigner$id,rxn_reactant$v1)
rxn_cellDesigner$reactants_name <- getMultipleReactionFormula(rxn_reactant$v1_name,rxn_reactant$v2,rxn_cellDesigner$id)

#products split, replace and combine
rxn_product <- splitAndCombine(rxn_cellDesigner$products,rxn_cellDesigner$id)
rxn_product$v1_name <- getSingleReactionFormula(met_cellDesigner$name,met_cellDesigner$id,rxn_product$v1)
rxn_cellDesigner$products_name <- getMultipleReactionFormula(rxn_product$v1_name,rxn_product$v2,rxn_cellDesigner$id)

#connect the reactants with products using "="
rxn_cellDesigner$reaction_detail <- paste(rxn_cellDesigner$reactants_name, rxn_cellDesigner$products_name, sep=" = ") 
for (i in 1:length(rxn_cellDesigner$id)){
  if (str_detect(rxn_cellDesigner$reaction_detail[i],";")){
    rxn_cellDesigner$reaction_detail[i] = str_replace_all(rxn_cellDesigner$reaction_detail[i],";"," + ")
  } else{
    rxn_cellDesigner$reaction_detail[i] = rxn_cellDesigner$reaction_detail[i]
  }
  
}
write.table(rxn_cellDesigner,"rxn_cellDesigner_detail.txt", row.names = FALSE, sep="\t")





#read the model
library(readxl)
yeastModel <- read_excel("yeastGEM for yeast map.xls", 
                         sheet = "Reaction List")
metabolite_yeastModel <- read_excel("yeastGEM for yeast map.xls", 
                                    sheet = "Metabolite List")## a new version of function-splitAndCombine

yeastModel$Reaction_without_compartment <- str_replace_all(yeastModel$Description_new,"\\[.*?\\]", "") %>%
  str_replace_all(.,"\\<\\=\\>", "=") %>%
  str_replace_all(.,"\\-\\>", "=")

write.table(yeastModel,"yeastModel_detail.txt", row.names = FALSE, sep="\t")





## give the kegg id for metabolites draft yeast map
splitAndCombine0 <- function(gene, rxn, sep0) { ##one rxn has several genes, this function was used to splite the genes
  
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
} ## new version of split and then combination
kegg_ALL_compound <- read.delim2("kegg_ALL_compound.txt", header = FALSE,sep="\t", stringsAsFactors = FALSE)
kegg_ALL_compound$V1 <- str_replace_all(kegg_ALL_compound$V1, "cpd:","")

kegg_compound_name <- splitAndCombine0(kegg_ALL_compound$V2, kegg_ALL_compound$V1,sep0=";")
colnames(kegg_compound_name) <- c('name','keggid')
kegg_compound_name_refine <- kegg_compound_name[which(duplicated(kegg_compound_name$keggid)==FALSE),]


## find the name from kegg database name based on keggID in yeastGEM model
metabolite_yeastModel$name_from_kegg <- getMultipleReactionFormula(kegg_compound_name_refine$name,kegg_compound_name_refine$keggid,metabolite_yeastModel$KEGGID)

for (i in 1:length(metabolite_yeastModel$`Metabolite name`)){
  if(!is.na(metabolite_yeastModel$name_from_kegg[i])){
    metabolite_yeastModel$name_based_kegg[i] <- metabolite_yeastModel$name_from_kegg[i]
  } else{
    metabolite_yeastModel$name_based_kegg[i] <- metabolite_yeastModel$`Metabolite description`[i]
  }
  
}


metabolite_yeastModel$`Metabolite description` <- str_replace_all(metabolite_yeastModel$`Metabolite description`,"\\[.*?\\]", "")
metabolite_yeastModel$name_based_kegg <- str_replace_all(metabolite_yeastModel$name_based_kegg,"\\[.*?\\]", "")


write.table(metabolite_yeastModel, "metabolite_yeastModel for mapping.txt", row.names = FALSE, sep = "\t")
