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


newGPR <- read_excel("data/newGPR with subsystem, reversibility and compartment_improved.xlsx")
rxn <- filter(newGPR, str_detect(newGPR$subsystem, "transport")==FALSE) 
rxn <- filter(rxn, str_detect(rxn$rxnID, "Exist")==FALSE) 
rxn$rxnID <- str_trim(rxn$rxnID, side = "both")

#find the rhea id based on keggid
#input the rhea id and keggid
rhea2kegg <- read_tsv('data/rhea2kegg_reaction.tsv')
rhea2kegg$RHEA_ID <- paste("RHEA:",rhea2kegg$RHEA_ID, sep = "")
rhea2kegg$MASTER_ID <- paste("RHEA:",rhea2kegg$MASTER_ID, sep = "")
rxn$RHEA_ID <- getSingleReactionFormula(rhea2kegg$MASTER_ID,rhea2kegg$ID,rxn$rxnID)


#find the rhea id based on biocycid
rhea2biocyc <- read_tsv('data/rhea2metacyc.tsv')
rhea2biocyc$RHEA_ID <- paste("RHEA:",rhea2biocyc$RHEA_ID, sep = "")
rhea2biocyc$MASTER_ID <- paste("RHEA:",rhea2biocyc$MASTER_ID, sep = "")
rxn$RHEA_ID2 <- getSingleReactionFormula(rhea2biocyc$MASTER_ID,rhea2biocyc$ID,rxn$rxnID)


#merge the two rhea ID

for (i in seq(length(rxn$ID0))){
  if(rxn$RHEA_ID[i]=='NA'){
    rxn$RHEA_ID[i] <- rxn$RHEA_ID2[i]
  } else{
    rxn$RHEA_ID[i] <- rxn$RHEA_ID[i]
  }
  
}




#find the metnetid based on keggid
reac_metnet <- read_tsv('data/reac_prop_metanetx.tsv')
reac_metnet$Source <- str_replace_all(reac_metnet$Source, "kegg:","")
rxn$metnetID <- getMultipleReactionFormula(reac_metnet$MNX_ID,reac_metnet$Source,rxn$rxnID)


#find the metnetid based on rheaid
rxn$RHEA_ID <-str_replace_all(rxn$RHEA_ID, "RHEA","rhea")
rxn$metnetID2 <- getMultipleReactionFormula(reac_metnet$MNX_ID,reac_metnet$Source,rxn$RHEA_ID)



#merge the metnetID
for (i in seq(length(rxn$metnetID2))){
  if(is.na(rxn$metnetID2[i])){
    rxn$metnetID2[i] <- rxn$metnetID[i]
  } else{
    rxn$metnetID2[i] <- rxn$metnetID2[i]
  }

  if(str_detect(rxn$rxnID[i],'MNXR' )){
    rxn$metnetID2[i] <- rxn$rxnID[i]
  } else{
    rxn$metnetID2[i] <- rxn$metnetID2[i]
  }

}


rxn_final <- filter(rxn, !is.na(rxn$metnetID2) | rxn$RHEA_ID != 'NA' )


# find reaction based on metnet id and rhea id
# input the metnet database

rxn_final$formula_metnet <- getSingleReactionFormula(reac_metnet$Description,reac_metnet$MNX_ID,rxn_final$metnetID2)
rxn_final$formula_metnet_balance <- getSingleReactionFormula(reac_metnet$Balance,reac_metnet$MNX_ID,rxn_final$metnetID2)
rxn_final$Equation <- getMultipleReactionFormula(reac_metnet$Equation,reac_metnet$MNX_ID,rxn_final$metnetID2)

# input the rhea database
rhea_reaction_summary <- read_csv("data/rhea reaction summary.csv")
rhea_reaction_summary$masterID <- paste("rhea:",rhea_reaction_summary$masterID, sep = "")
rxn_final$formula_rhea <- getMultipleReactionFormula(rhea_reaction_summary$formula,rhea_reaction_summary$masterID,rxn_final$RHEA_ID)

for (i in seq(length(rxn_final$ID0))){
  if(rxn_final$formula_metnet[i] =='NA'){
    rxn_final$formula_metnet[i] <- rxn_final$formula_rhea[i]
    rxn_final$formula_metnet_balance[i] <- 'true'
    } else{
    rxn_final$formula_metnet[i] <- rxn_final$formula_metnet[i] 
    rxn_final$formula_metnet_balance[i] <- rxn_final$formula_metnet_balance[i]
    
    }
  
}

write.table(rxn_final, "result/rxn_final.txt", row.names = FALSE, sep = "\t")






##find the standard formula of metabolite based on reaction compostions from metnetx database
rxn_final0 <- select(rxn_final, ID0, Equation)
rxn_list <- str_split(rxn_final0$Equation, " = ")
for (i in seq(length(rxn_list))){
  rxn_list[[i]][1] <- paste("reactant",rxn_list[[i]][1],sep = "@@" )
  rxn_list[[i]][1] <- paste(rxn_final0$ID0[i],rxn_list[[i]][1],sep = "@" )
  
  rxn_list[[i]][2] <- paste("product",rxn_list[[i]][2],sep = "@@" )
  rxn_list[[i]][2] <- paste(rxn_final0$ID0[i],rxn_list[[i]][2],sep = "@" )
  }

rxn_unlist <- unlist(rxn_list)

rxn0_list <- str_split(rxn_unlist, "@@")

ss1 <- vector()
ss2 <- vector()
for (i in seq(length(rxn0_list))){
  ss1[i] <- rxn0_list[[i]][1]
  ss2[i] <- rxn0_list[[i]][2]
  }
ss2_list <- str_split(ss2, " \\+ ")

for (i in seq(length(rxn0_list))){
ss2_list[[i]] <- paste(ss1[[i]], ss2_list[[i]], sep = "@@")
}

ss2_unlist <- unlist(ss2_list)

ss3_list <- str_split(ss2_unlist, "@@")


ss4 <- vector()
ss5 <- vector()
for (i in seq(length(ss3_list))){
  ss4[i] <- ss3_list[[i]][1]
  ss5[i] <- ss3_list[[i]][2]
}

mnx_id <- str_replace_all(ss5, "@MNXD[:alnum:]","") %>%
  str_replace_all(.,"[:digit:] ","") %>%
  str_replace_all(.,"\\(n\\)","") %>%
  str_replace_all(.,"\\(n\\+1\\)","")


rxn_met <- data.frame(reaction = ss4, mnxID = mnx_id, stringsAsFactors = FALSE)
rxn_met <- rxn_met %>% separate(.,reaction, into = c('ID','compostion'), sep = "@")
chem_prop <- read_tsv('data/chem_prop_metanetx.tsv')
rxn_met$Formula <- getMultipleReactionFormula(chem_prop$Formula,chem_prop$MNX_ID,rxn_met$mnxID)
rxn_met$Charge <- getMultipleReactionFormula(chem_prop$Charge,chem_prop$MNX_ID,rxn_met$mnxID)
rxn_met$Description <- getMultipleReactionFormula(chem_prop$Description,chem_prop$MNX_ID,rxn_met$mnxID)

write.table(rxn_met, "rxn_obtain the metabolite based on reaction id.txt", row.names = FALSE, sep = "\t")








##find the standard formula of metabolite based on chebiID, keggID
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

#find the metabolite formula based on mnx_id
chem_prop <- read_tsv('data/chem_prop_metanetx.tsv')
newGPR20$Formula <- getMultipleReactionFormula(chem_prop$Formula,chem_prop$MNX_ID,newGPR20$MNX_ID)
newGPR20$Charge <- getMultipleReactionFormula(chem_prop$Charge,chem_prop$MNX_ID,newGPR20$MNX_ID)

write.table(newGPR20, "result/met_with metnetxID and formula.txt", row.names = FALSE, sep = "\t")

