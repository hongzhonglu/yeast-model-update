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

#an example
#newGPR <- read_excel("data/newGPR with subsystem, reversibility and compartment_improved.xlsx")
#rxn <- filter(newGPR, str_detect(newGPR$subsystem, "transport")==FALSE) 
#rxn <- filter(rxn, str_detect(rxn$rxnID, "Exist")==FALSE) 
#rxn$rxnID <- str_trim(rxn$rxnID, side = "both")



#input the rection information with RxnID from kegg databas
#input data should have a column named 'rxnID', in this column, there are reactionID from kegg database,like R06513
rxn <- read.delim2("data/newRxn_panGenome_kegg_web.txt", stringsAsFactors = FALSE)


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

for (i in seq(length(rxn$rxnID))){
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


# now some rxnID from kegg could be mapped onto Rhea or metnetx database, while some can't
# in the above first case, the metabolite formula could be obtained based on Rhea or metanetx
# while for the second case, keggid can be obtained for each reactions


#rxn_final <- filter(rxn, !is.na(rxn$metnetID2) | rxn$RHEA_ID != 'NA' )
rxn_final <- rxn #sometimes need some filters

# find reaction based on metnet id and rhea id
# input the metnet database
rxn_final$formula_metnet <- getSingleReactionFormula(reac_metnet$Description,reac_metnet$MNX_ID,rxn_final$metnetID2)
rxn_final$formula_metnet_balance <- getSingleReactionFormula(reac_metnet$Balance,reac_metnet$MNX_ID,rxn_final$metnetID2)
rxn_final$Equation <- getMultipleReactionFormula(reac_metnet$Equation,reac_metnet$MNX_ID,rxn_final$metnetID2)

# input the rhea database
rhea_reaction_summary <- read_csv("data/rhea reaction summary.csv")
rhea_reaction_summary$masterID <- paste("rhea:",rhea_reaction_summary$masterID, sep = "")
rxn_final$formula_rhea <- getMultipleReactionFormula(rhea_reaction_summary$formula,rhea_reaction_summary$masterID,rxn_final$RHEA_ID)

for (i in seq(length(rxn_final$rxnID))){
  if(rxn_final$formula_metnet[i] =='NA'){
    rxn_final$formula_metnet[i] <- rxn_final$formula_rhea[i]
    #rxn_final$formula_metnet_balance[i] <- 'true'
    } else{
    rxn_final$formula_metnet[i] <- rxn_final$formula_metnet[i] 
    rxn_final$formula_metnet_balance[i] <- rxn_final$formula_metnet_balance[i]
    
    }
  
}


# increase the reveserbility check based on BiGG or Rhea database
# for example based on rhea database
head(rxn_final$formula_rhea)
reverse <- str_count(rxn_final$formula_rhea, ";")
for (i in seq(length(rxn_final$rxnID))){
  if(!is.na(rxn_final$formula_rhea[i]) & reverse[i] >= 1 ){
    rxn_final$reverse[i] <- 'yes'
  } 
  else if (!is.na(rxn_final$formula_rhea[i]) & reverse[i] != 1) {
    rxn_final$reverse[i] <- 'no'
    
  } else{
    rxn_final$reverse[i] <- 'not sure'
  }
  
}



write.table(rxn_final, "result/new reactions from kegg based on panGenome.txt", row.names = FALSE, sep = "\t")





##find the standard formula of metabolite based on reaction compostions from metnetx database
rxn_final0 <- select(rxn_final, rxnID, Equation)
colnames(rxn_final0) <- c('ID0', 'Equation')

splitRxnToMetabolite <- function(reationFrame, sep0){
   #reationFrame <- rxn_final0
   #input:
   #dataframe: column 1--ID0; column 2--Equation, which is the reaction formula
   rxn_list <- str_split(reationFrame$Equation, sep0)
   for (i in seq(length(rxn_list))){
      rxn_list[[i]][1] <- paste("reactant",rxn_list[[i]][1],sep = "@@" )
      rxn_list[[i]][1] <- paste(reationFrame$ID0[i],rxn_list[[i]][1],sep = "@" )
  
      rxn_list[[i]][2] <- paste("product",rxn_list[[i]][2],sep = "@@" )
      rxn_list[[i]][2] <- paste(reationFrame$ID0[i],rxn_list[[i]][2],sep = "@" )
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
   rxn_met <- data.frame(reaction = ss4, MetID = ss5, stringsAsFactors = FALSE)
   rxn_met <- rxn_met %>% separate(.,reaction, into = c('ID','compostion'), sep = "@")
   #remove the coefficient for each metabolites, which will be easy to standardize this metabolite
   rxn_met$MetID <- str_replace_all(rxn_met$MetID, "[:digit:] ","") %>%
     str_replace_all(.,"\\(n\\)","") %>%
     str_replace_all(.,"\\(n\\+1\\)","") %>%
     str_replace_all(.,"\\(n\\-2\\)","")
 return(rxn_met)
}


#remove the compartment information from metnet database
rxn_met <- splitRxnToMetabolite(rxn_final0, sep0 = " = ")
rxn_met$MetID <- str_replace_all(rxn_met$MetID, "@MNXD[:alnum:]","")
chem_prop <- read_tsv('data/chem_prop_metanetx.tsv')
rxn_met$Formula <- getMultipleReactionFormula(chem_prop$Formula,chem_prop$MNX_ID,rxn_met$MetID)
rxn_met$Charge <- getMultipleReactionFormula(chem_prop$Charge,chem_prop$MNX_ID,rxn_met$MetID)
rxn_met$Description <- getMultipleReactionFormula(chem_prop$Description,chem_prop$MNX_ID,rxn_met$MetID)
rxn_met0 <- filter(rxn_met, rxn_met$MetID != 'NA')






##find the standard formula of metabolite based on chebiID, keggID
rxn_needMet <- filter(rxn_met, rxn_met$MetID =='NA')
rxn_needMet_list <- unique(rxn_needMet$ID)
index0 <- which(rxn$rxnID %in% rxn_needMet_list)
rxn2 <- rxn[index0,]
rxn2_for_met <- select(rxn2, rxnID, formula)
colnames(rxn2_for_met) <- c('ID0', 'Equation')
rxn_met2 <- splitRxnToMetabolite(rxn2_for_met, sep0 = " <=> ")

#based on name find the keggid for each metabolites
#prepare the name from kegg database
name_kegg <- read.csv2("data/compound_kegg.txt", sep="\t", stringsAsFactors = FALSE, header = FALSE)
name_kegg$V1 <- str_replace_all(name_kegg$V1,"cpd:","")
name_kegg0 <- splitAndCombine(name_kegg$V2,name_kegg$V1,sep0 = ";")
colnames(name_kegg0) <- c("metabolite_name","KEGGID")
name_kegg0$name_unify <-  str_to_lower(name_kegg0$metabolite_name) %>%
  str_trim(., side = "both")
rxn_met2$name_unify <- str_to_lower(rxn_met2$MetID) %>%
  str_trim(., side = "both")
rxn_met2$keggid <- getMultipleReactionFormula(name_kegg0$KEGGID, name_kegg0$name_unify, rxn_met2$name_unify)



#find chebiID based on keggid
kegg_chebi <- read.table('data/kegg_chebi.txt', sep ="\t", stringsAsFactors = FALSE)
kegg_chebi$V1 <- str_replace_all(kegg_chebi$V1, "cpd:", "")
kegg_chebi$V2 <- str_replace_all(kegg_chebi$V2, "chebi:","CHEBI:")
rxn_met2$chebiid <- getSingleReactionFormula(kegg_chebi$V2,kegg_chebi$V1,rxn_met2$keggid)

#find the mnxid based on chebiid
metanetx_metref <- read_tsv('data/chem_xref.tsv')
metanetx_metref$XREF <- str_replace_all(metanetx_metref$XREF, "chebi:","CHEBI:")
rxn_met2$MNX_ID <- getMultipleReactionFormula(metanetx_metref$MNX_ID,metanetx_metref$XREF,rxn_met2$chebiid)

#find the metabolite formula based on mnx_id
rxn_met2$Formula <- getMultipleReactionFormula(chem_prop$Formula,chem_prop$MNX_ID,rxn_met2$MNX_ID)
rxn_met2$Charge <- getMultipleReactionFormula(chem_prop$Charge,chem_prop$MNX_ID,rxn_met2$MNX_ID)


#combine the metabolite formula found from different ways
#find the chebiid and keggid based on metnetxID information
metanetx_metref0 <- filter(metanetx_metref, str_detect(metanetx_metref$XREF,"CHEBI"))
metanetx_metref1 <- filter(metanetx_metref, str_detect(metanetx_metref$XREF,"kegg"))
rxn_met0$chebiid <- getMultipleReactionFormula(metanetx_metref0$XREF,metanetx_metref0$MNX_ID,rxn_met0$MetID)
rxn_met0$keggid <- getSingleReactionFormula(metanetx_metref1$XREF,metanetx_metref0$MNX_ID,rxn_met0$MetID)
ID1 <- paste('M', 1:123, sep = "")
mm <- splitAndCombine(rxn_met0$chebiid,ID1, sep0 = ";")
mm$keggid <- getSingleReactionFormula(kegg_chebi$V1,kegg_chebi$V2,mm$v1)
mm <- filter(mm, str_detect(mm$keggid, "C"))
rxn_met0$ID1 <- ID1
rxn_met0$keggid <- getSingleReactionFormula(mm$keggid,mm$v2,rxn_met0$ID1)
#replace the old chebiid WITH new from kegg database
rxn_met0$chebiid <- getSingleReactionFormula(mm$v1,mm$v2,rxn_met0$ID1)


ss <- colnames(rxn_met0)
ss[3] <- "MNX_ID"
colnames(rxn_met0) <- ss
ss <- ss[-9]
ss2 <- colnames(rxn_met2)
ss2[3] <- "Description"
colnames(rxn_met2) <- ss2

rxn_met2 <- select(rxn_met2, ss)
rxn_met0 <- select(rxn_met0,ss)
rxn_met_final <- rbind.data.frame(rxn_met0, rxn_met2)

write.table(rxn_met_final, "result/metabolite standardization for new reactions from kegg based on panGenome.txt", row.names = FALSE, sep = "\t")











