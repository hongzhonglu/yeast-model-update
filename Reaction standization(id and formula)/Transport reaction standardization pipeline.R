library(readxl)
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




#summarize the reaction from kegg database
reaction_kegg  <- read.delim2("data/reaction from kegg.txt", sep = "\t", stringsAsFactors = FALSE)
for (i in seq(length(reaction_kegg$rxnID))){
  if(reaction_kegg$X.1[i] == "" & reaction_kegg$X.2[i] == "" & reaction_kegg$X.3[i] == ""){
    reaction_kegg$reaction[i] <- reaction_kegg$X[i]
    reaction_kegg$name[i] <- "NONE"
  } else if(reaction_kegg$X.2[i] == "" & reaction_kegg$X.3[i] == ""){
    reaction_kegg$reaction[i] <- reaction_kegg$X.1[i]
    reaction_kegg$name[i] <- reaction_kegg$X[i]
  } else if(reaction_kegg$X.3[i] == ""){
    reaction_kegg$reaction[i] <- reaction_kegg$X.2[i]
    reaction_kegg$name[i] <- paste(reaction_kegg$X[i],reaction_kegg$X.1[i],sep = ";")
  } else {
    reaction_kegg$reaction[i] <- reaction_kegg$X.3[i]
    reaction_kegg$name[i] <- paste(reaction_kegg$X[i],reaction_kegg$X.1[i], reaction_kegg$X.2[i],sep = ";")
  }
}

reaction_kegg0  <- select(reaction_kegg, rxnID, name, reaction)
reaction_kegg0$reaction <- str_trim(reaction_kegg0$reaction) %>%
  str_replace_all(.,"NAD\\+", "NAD\\(\\+\\)") %>%
  str_replace_all(.,"H\\+", "H\\(\\+\\)") %>%
  str_replace_all(.,"NADP\\+", "NADP\\(\\+\\)") %>%
  str_to_lower(.)

write.table(reaction_kegg0, "result/reaction_kegg summary.txt", row.names = FALSE, sep = "\t")



##find the keggid for the reactions automatically
newGPR_with_subsystem <- read_excel("data/newGPR with subsystem, reversibility and compartment_improved.xlsx")
newGPR_with_subsystem <- select(newGPR_with_subsystem, ID0, reaction, subsystem)
newGPR_with_subsystem$reaction <- str_trim(newGPR_with_subsystem$reaction) %>%
  str_replace_all(.,"  \\+ ", " \\+ ") %>%
  str_replace_all(.,"  <", " <") %>%
  str_replace_all(.,"\\)<", "\\) <") 

newGPR_with_subsystem$reaction0 <- str_to_lower(newGPR_with_subsystem$reaction)
newGPR_with_subsystem$rxnID <- getSingleReactionFormula(reaction_kegg0$rxnID,reaction_kegg0$reaction,newGPR_with_subsystem$reaction0)
write.table(newGPR_with_subsystem, "result/newGPR_with_subsystem.txt", row.names = FALSE, sep = "\t")






##transport reaction standization
##find the transport reaction id from metanet dabase
##it should be noted that for one metabolite transport, there exist different reaction formula
##thus the final transport reaction need manual check by compare the related gene annotation


newGPR2 <- read_excel("data/newGPR with subsystem, reversibility and compartment_improved.xlsx", 
                                                                           sheet = "rxn split for matlab")
newGPR20 <- select(newGPR2, rxnID, Standard_metName)
met <- read_excel("data/newGPR with subsystem, reversibility and compartment_improved.xlsx", 
                                                                           sheet = "Metabolite list")
newGPR20$keggid <- getSingleReactionFormula(met$keggID,met$Standard_metName,newGPR20$Standard_metName)
newGPR20$chebiid <- getSingleReactionFormula(met$ChEBI_ID,met$Standard_metName,newGPR20$Standard_metName)

#find chebiID based on keggid
kegg_chebi <- read.table('data/kegg_chebi.txt', sep ="\t", stringsAsFactors = FALSE)
kegg_chebi$V1 <- str_replace_all(kegg_chebi$V1, "cpd:", "")
kegg_chebi$V2 <- str_replace_all(kegg_chebi$V2, "chebi:","CHEBI:")
newGPR20$chebiid2 <- getSingleReactionFormula(kegg_chebi$V2,kegg_chebi$V1,newGPR20$keggid)
for (i in seq(length(newGPR20$chebiid2))){
  if(newGPR20$chebiid2[i]=="NA"){
    newGPR20$chebiid2[i] <- newGPR20$chebiid[i]
  } else{
    newGPR20$chebiid2[i] <- newGPR20$chebiid2[i]
  }
}

#find the mnxid based on chebiid
metanetx_metref <- read_tsv('data/chem_xref.tsv')
metanetx_metref$XREF <- str_replace_all(metanetx_metref$XREF, "chebi:","CHEBI:")
newGPR20$MNX_ID <- getMultipleReactionFormula(metanetx_metref$MNX_ID,metanetx_metref$XREF,newGPR20$chebiid2)

# choose transport
newGPR_with_subsystem <- read_excel("data/newGPR_with_subsystem.xlsx")
newGPR_transport <- filter(newGPR_with_subsystem, str_detect(newGPR_with_subsystem$subsystem,"transport"))
index0 <- which(newGPR20$rxnID %in% newGPR_transport$ID0 ==TRUE)
newGPR20_transport <- newGPR20[index0,]

#the currency metabolite is connected with so many reactions and should be removed.
#find the currency met and remove them
met0 <- newGPR20_transport %>% 
  count(Standard_metName) %>%
  arrange(.,desc(n))

currency_met <- c(met0$Standard_metName[1:4],'phosphate')

index1 <- which(newGPR20_transport$Standard_metName %in% currency_met ==FALSE)
newGPR20_transport <- newGPR20_transport[index1,]
newGPR20_transport <- newGPR20_transport[!duplicated(newGPR20_transport$rxnID),]


#find reaction connect with each metabolite
metanetx_rxn <- read_tsv("data/reac_prop_metanetx.tsv")
rxn_id_list <- list()
rxn_list <- list()
rxn_list0 <- list() # strore the reactions using the MNX_ID
for (i in seq(length(newGPR20_transport$MNX_ID))){
  rxn_id_list[[i]] <- which(str_detect(metanetx_rxn$Equation,paste(newGPR20_transport$MNX_ID[i],"@",sep = "")) ==TRUE)
}


names(rxn_id_list) <- newGPR20_transport$rxnID


for (i in seq(length(newGPR20_transport$MNX_ID))){
  rxn_list[[i]] <- metanetx_rxn$Description[rxn_id_list[[i]]]
  rxn_list0[[i]] <- metanetx_rxn$Equation[rxn_id_list[[i]]]
  rxn_list0[[i]] <- paste(newGPR20_transport$rxnID[i],rxn_list0[[i]],sep = "@@")
}

names(rxn_list0) <- newGPR20_transport$rxnID


rxn_withMet <- unlist(rxn_list0)
rxn_withMet0 <- data.frame(met = rxn_withMet)
rxn_withMet0 <- rxn_withMet0 %>% separate(met,into = c("met","reaction"), sep = "@@")

index3 <- which(str_detect(rxn_withMet0$reaction, "BIOMASS")==FALSE)
rxn_withMet0 <- rxn_withMet0[index3,]

rxn_withMet0 <- rxn_withMet0[!duplicated(rxn_withMet0$reaction),]
rxn_withMet0$formula <- getSingleReactionFormula(metanetx_rxn$Description,metanetx_rxn$Equation,rxn_withMet0$reaction)
rxn_withMet0$formula0  <- rxn_withMet0$formula

rxn_withMet0$MNX_ID <- getSingleReactionFormula(metanetx_rxn$MNX_ID,metanetx_rxn$Equation,rxn_withMet0$reaction)


rxn_withMet0 <- rxn_withMet0 %>% separate(formula0,into = c("reactant","product"), sep = "=")

reactant <- str_split(rxn_withMet0$reactant," \\+ ")
for (i in seq(length(reactant))){
  reactant[[i]] <- str_trim(reactant[[i]], side = "both")
}

common_met <- list()
len1 <- vector()
product <- str_split(rxn_withMet0$product," \\+ ")
for (i in seq(length(product))){
  product[[i]] <- str_trim(product[[i]], side = "both")
  common_met[[i]] <- intersect(reactant[[i]], product[[i]])
  len1[i] <- length(common_met[[i]])
  }

rxn_withMet0$common_met_num <- len1

rxn_withMet0 <- filter(rxn_withMet0, len1 > 0)
write.table(rxn_withMet0,"result/rxn_withMet0.txt", row.names = FALSE, sep = "\t")





