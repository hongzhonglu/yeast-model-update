library(stringr)
library(tidyverse)
library(readxl)
library(igraph)
library(networkD3)
library(hongR)

########small task###update the metabolites formula
#read the model

rxn <- read_excel("data/yeastGEM_latest version.xls", 
                                             sheet = "Reaction List")
metabolite <-  read_excel("data/yeastGEM_latest version.xls", 
                          sheet = "Metabolite List")

splitAndCombine0 <- function(gene, rxn, sep0) { ##one rxn has several genes, this function was used to splite the genes
  
  gene <- str_split(gene, sep0)
  tt<- length(gene)
  gene0 <- list()
  for (i in 1:tt){
    gene0[[i]] <- paste(rxn[i], gene[[i]], sep = "@@@")
    
  }
  
  gene1 <- unlist(gene0)
  gene2 <- str_split(gene1, "@@@" )
  rxnGene <- data.frame(v1=character(length(gene2)),stringsAsFactors = FALSE)
  tt1 <- length(gene2)
  for (j in 1:tt1){
    rxnGene$v1[j] <- gene2[[j]][2]
    rxnGene$v2[j] <- gene2[[j]][1]
  }
  
  return(rxnGene)
}

rxn_split <- splitAndCombine0(rxn$Reaction,rxn$Abbreviation,sep=" ")
rxn_split$v3 <- getSingleReactionFormula(metabolite$`Metabolite description`,metabolite$`Metabolite name`,rxn_split$v1)

for (i in 1:length(rxn_split$v2)){
  if(rxn_split$v3[i]=="NA"){
    rxn_split$v3[i] <- rxn_split$v1[i]
  } else{
    rxn_split$v3[i] <- rxn_split$v3[i]
  }
  

}
rxn$Description_new <- getMultipleReactionFormula(rxn_split$v3,rxn_split$v2,rxn$Abbreviation)
rxn$Description_new <- str_replace_all(rxn$Description_new,";"," ")
write.table(rxn, "yeastGEM with update reaction formula.txt", row.names = FALSE, sep = "\t")
########small task###update the metabolites formula


rxn_metabolite_mapping <- select(rxn,Abbreviation, Reaction)
rxn_metabolite_mapping$Reaction <- str_replace_all(rxn_metabolite_mapping$Reaction, "\\-\\>","\\<\\=\\>" )


getMappingReaction <- function(reaction,id) {
  # input: reaction equation
  # input: reaction id
  # this function is used to establish the relation between metabolite according to certain reaction
  
rxn_metabolite_split0 <- splitAndCombine0(reaction,id,sep=" \\<\\=\\> ")
reactants <- str_split(rxn_metabolite_split0$v1[1], " ")
reactants <- unlist(reactants)
reactants_index <- which(str_detect(reactants,"\\[")==TRUE)
reactants <- reactants[reactants_index]

products <- str_split(rxn_metabolite_split0$v1[2]," ")
products <- unlist(products)
products_index <- which(str_detect(products,"\\[")==TRUE)
products <- products[products_index]


m <- length(reactants)
n <- length(products)
mn <- c()
for (i in seq(m)){
  for(j in seq(n)){
    mn <- c(paste(reactants[i],products[j],sep = ";"),mn)
  }
}

return(mn)
}


## process to establish the metabolite mapping relation
ss <- list()
ss0 <- list()
for (i in 1:3496){
ss[[i]] <- getMappingReaction(rxn_metabolite_mapping$Reaction[i],rxn_metabolite_mapping$Abbreviation[i])
ss0[[i]] <- paste(rxn_metabolite_mapping$Abbreviation[i],ss[[i]],sep = ";" )
}

mappedMetabolites <- unlist(ss0)
head(mappedMetabolites)

ss1 <- str_split(mappedMetabolites,";")
ss2 <- data.frame(id = character(length(ss1)), stringsAsFactors = FALSE)
for (i in seq(length(ss1))){
  ss2$id[i] <- ss1[[i]][1]
  ss2$from[i] <- ss1[[i]][2]
  ss2$to[i] <- ss1[[i]][3]
  ss2$combine[i] <- paste(ss1[[i]][2],ss1[[i]][3])
}

ss2 <- ss2[which(ss2$to!=""),]

ss3 <- ss2[!duplicated(ss2$combine),] # remove the duplicated relation of metabolites ?

ss3 <- select(ss3, from, to)

simpleNetwork(ss3[1:100,], fontSize = 18, opacity = 1, zoom = TRUE,fontFamily = "sans-serif")


