library(stringr)
library(tidyverse)
library(readxl)
library(igraph)
library(networkD3)

###match
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




#part one
#prepare the reaction format
########small task###update the metabolites formula
#read the model

rxn <- read_excel("yeastGEM_latest version.xls", 
                                             sheet = "Reaction List")
metabolite <-  read_excel("yeastGEM_latest version.xls", 
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


# give the type: reactant and product for each metabolite
#for r_0001
 rxn_split$v3 <- str_replace_all(rxn_split$v3,"->", "<=>")
 
 getMetType <- function(rxnID, rxn_split){
 ss <- which(rxn_split$v2 %in% rxnID==TRUE)
 ss_split <- which(rxn_split$v2 %in% rxnID ==TRUE & rxn_split$v3 %in% "<=>" == TRUE)
 met_type1 <- vector()
 met_type2 <- vector()
 for (i in ss[1]:(ss_split-1)){
   met_type1[i-ss[1]+1] <- "reactant"
 }
 
 for (i in (ss_split):ss[length(ss)]){
   met_type2[i-ss_split+1] <- "product"
 }
 
 met_type <- c(met_type1,met_type2)
 return(met_type)
}


rxn_ID <- unique(rxn_split$v2)
met_type0 <- vector() 
for (i in 1:length(rxn_ID)){
  met_type0 <- c(met_type0, getMetType(rxn_ID[i],rxn_split))
  
}
rxn_split$type <- met_type0
rxn_split_refine <-rxn_split[which(str_detect(rxn_split$v3,"]")==TRUE),]

rxn_split_refine$metabolit_formula <- getSingleReactionFormula(metabolite$`Metabolite formula`,metabolite$`Metabolite description`,rxn_split_refine$v3)
rxn_split_refine$v3 <- str_replace_all(rxn_split_refine$v3,"<=>","->")




###get carbon number for each metabolite
getMetaboliteComposition <- function(MET,type="C"){
if (str_detect(MET, paste(type,"[0-9]+", sep=""))){
 carbon_number<- str_extract_all(MET, paste(type,"[0-9]+", sep=""))
} else if (str_detect(MET, "C") & str_detect(MET, paste(type,"[0-9]+", sep=""))==FALSE ){
carbon_number <- str_extract_all(MET, "C")  
} else{
 carbon_number <- ""
}
  return(carbon_number[[1]][1])
}

getCompositionNum <- function(MET,type="C"){
  if (str_detect(MET, paste(type,"[0-9]+", sep=""))){
    carbon_number<- str_extract_all(MET, "[0-9]+")
  } else if (str_detect(MET, "C") & str_detect(MET, paste(type,"[0-9]+", sep=""))==FALSE ){
    carbon_number <- "1"
  } else{
    carbon_number <- ""
  }
  return(carbon_number[[1]][1])
}

for (i in 1:length(rxn_split_refine$v2)){
  rxn_split_refine$carbonCompostion[i] <- getMetaboliteComposition(rxn_split_refine$metabolit_formula[i])

}

for (i in 1:length(rxn_split_refine$v2)){
  rxn_split_refine$carbonNumber[i] <- getCompositionNum(rxn_split_refine$carbonCompostion[i])
  
}


##remove currency metabolites
currency_metabolites <- c("H+", "H2O", "ATP", "phosphate", "coenzyme A", "ADP", "diphosphate", "NADP(+)", "NADPH", "carbon dioxide", "AMP","NAD", "CMP", "NADH")     
rxn_split_refine$simple_name <- str_replace_all(rxn_split_refine$v3,"\\[.*?\\]", "")
which(rxn_split_refine$simple_name %in% currency_metabolites ==FALSE)
rxn_split_refine <- rxn_split_refine[which(rxn_split_refine$simple_name %in% currency_metabolites ==FALSE),]

#remove the reactions with only one product or only one reactant
analysis_rxn <- rxn_split_refine %>% 
  count(v2) %>% ## calculate the number of each metabolite
  arrange(., desc(n)) ## order the metabolites based on the number

rxn_with_one_met <- which(analysis_rxn$n ==1)
rxn_remove <- analysis_rxn$v2[rxn_with_one_met]
which(rxn_split_refine$v2 %in% rxn_remove ==FALSE)
rxn_split_refine <- rxn_split_refine[which(rxn_split_refine$v2 %in% rxn_remove ==FALSE),]

#further remove the reactions with only one product or only one reactant
rxn_index <- list()
rxn_unique <- unique(rxn_split_refine$v2)
for (i in 1:length(rxn_unique)){
  rxn_index[[i]] <- which(rxn_split_refine$v2 %in% rxn_unique[i])
  
}

rxn_metabotite_type <- list()
for (i in 1:length(rxn_unique)){
  rxn_metabotite_type[[i]] <- rxn_split_refine$type[rxn_index[[i]]]
}

ss<- vector()
for (i in 1:length(rxn_unique)){
  if("reactant" %in% rxn_metabotite_type[[i]] & "product" %in% rxn_metabotite_type[[i]]){
    
    ss[i] <- i
  } else{
    ss[i] <- NA
 }
  
}

rxn_final <- rxn_unique[!is.na(ss)]
rxn_split_refine <- rxn_split_refine[which(rxn_split_refine$v2 %in% rxn_final == TRUE),]




## define the base reactant and product for cellDesigner
## it should be noted that the currency metabolites is not considered in this phase
rxn_split_refine$note <- ""
rownames(rxn_split_refine) <- 1:length(rxn_split_refine$v1)

## define the base reactant and product for cellDesigner
DefineBaseMetabolite <- function(rxnID, rxn_split_refine){
  ## input: rxnID: a reaction ID
  ## input: rxn_split_refine: datafram contains the relation of rxnID and metabolite, as well as metabolite formula and 
  ## and carbon number
  
  ## find the base reactant (with largest Carbon number)
  if(length(which(rxn_split_refine$v2 %in% rxnID ==TRUE & rxn_split_refine$type =="reactant")) !=0){
    rxn_reactant <- which(rxn_split_refine$v2 %in% rxnID ==TRUE & rxn_split_refine$type =="reactant")
    ss <- max(rxn_split_refine$carbonNumber[rxn_reactant],na.rm = TRUE)
    tt <- which(rxn_split_refine$carbonNumber[rxn_reactant] == ss)
    nn <- vector()
    nn[1] <- rxn_reactant[1]-1+tt
    
  } else{
    nn[1] <- ""
  }
  ## find the base product (with largest carbon nunber)
  if (length(which( rxn_split_refine$v2 %in% rxnID ==TRUE & rxn_split_refine$type =="product")) !=0){
    rxn_product <- which( rxn_split_refine$v2 %in% rxnID ==TRUE & rxn_split_refine$type =="product")
    ss0 <- max(rxn_split_refine$carbonNumber[rxn_product],na.rm = TRUE)
    tt0 <- which(rxn_split_refine$carbonNumber[rxn_product] == ss0)
    nn[2] <- rxn_product[1]-1+tt0
    
  } else {
    nn[2] <-""
  }
  
  return(nn)
}

metabolite_type <- list()
rxn_ID <- unique(rxn_split_refine$v2)
for (i in 1:length(rxn_ID)){
  metabolite_type[[i]] <- DefineBaseMetabolite(rxn_ID[i],rxn_split_refine)
}
metabolite_type <- unlist(metabolite_type)

for (i in 1:length(rxn_split_refine$v2)){
 if (i %in% metabolite_type){
  rxn_split_refine$note[i] <- "base"
} else {
  rxn_split_refine$note[i] <- ""
}
}


#classify the reactions based on subsystems
rxn_system <- select(rxn,Abbreviation,Subsystem_new)
index00 <- which(str_detect(rxn_system$Subsystem_new,"transport")==TRUE)
rxn_system$Subsystem_new[index00] <-"transport"
subsystem <- rxn_system %>% 
  count(Subsystem_new)

# give subsystem for each reaction
rxn_split_refine$subsystem <- getSingleReactionFormula(rxn_system$Subsystem_new,rxn_system$Abbreviation,rxn_split_refine$v2)

index0 <- which(str_detect(rxn_split_refine$subsystem, "glycolysis")==TRUE)
index1 <- which(str_detect(rxn_split_refine$subsystem, "citrate cycle")==TRUE)
index2 <- which(str_detect(rxn_split_refine$subsystem, "pentose phosphate pathway")==TRUE)
index_combine <- c(index0, index1, index2)
rxn_core_carbon <- rxn_split_refine[index_combine,]
met_core_carbon <- unique(rxn_core_carbon$v3)

# find the transport reactions to connect the gap in the above systerm
index_transport <- which(str_detect(rxn_split_refine$subsystem, "transport")==TRUE)
rxn_transport <- rxn_split_refine[index_transport,]
rxn_transport_id <- unique(rxn_transport$v2)

getConnectedTransport <- function (id,rxn_transport_id, met_core_carbon){
  met_transport_index <- which(rxn_transport$v2 %in% rxn_transport_id[id] )
  met_transport <- rxn_transport$v3[met_transport_index]
  ss <- vector()
  ss <- met_transport %in% met_core_carbon
  tt <- sum(ss)
  if (tt == length(met_transport)){
     mm <- id
  } else{
  mm <- NA
  
 }
return(mm)
}

connected_rxn <- vector()
for (i in 1:length(rxn_transport_id)){
  connected_rxn[i] <- getConnectedTransport(id=i, rxn_transport_id, met_core_carbon)

}

connect_rxn0 <- connected_rxn[!is.na(connected_rxn)]



##add the connected transport reactions
connect_rxn0
trasport_choosed <- rxn_transport_id[connect_rxn0]
trasport_rxn_choosed <- rxn_transport[which(rxn_transport$v2 %in% trasport_choosed==TRUE),]
rxn_core_carbon <- rbind.data.frame(rxn_core_carbon,trasport_rxn_choosed)



#part two
#prepare the metabolites formula

met_core_carbon <- select(rxn_core_carbon,v2,v3,type)
colnames(met_core_carbon) <- c("rxnID","name","type")
unique_met_core_carbon <- unique(met_core_carbon$name)
unique_met_annotation <- data.frame(species=character(length(unique_met_core_carbon)),stringsAsFactors = FALSE)
unique_met_annotation$name <- unique_met_core_carbon 
unique_met_annotation$species <- paste("s", 1:length(unique_met_core_carbon), sep="")
unique_met_annotation$id <-paste("sa", 1:length(unique_met_core_carbon), sep="") 
unique_met_annotation$MetaID <- paste("CDMT0000", 1:length(unique_met_core_carbon), sep="") 
unique_met_annotation$x <-""
unique_met_annotation$y <-""
unique_met_annotation$rxnID <- getSingleReactionFormula(met_core_carbon$rxnID,met_core_carbon$name,unique_met_annotation$name)
unique_met_annotation <- select(unique_met_annotation,rxnID,id,species,name,x,y,MetaID)
write.table(unique_met_annotation, "met_core_carbon.txt", row.names = FALSE, sep="\t")



#prepare the rxn formula
rxn_core_carbon_cellD <- select(rxn_core_carbon, v2,v3,type,note)
rxn_core_carbon_cellD$species <- getSingleReactionFormula(unique_met_annotation$species,unique_met_annotation$name,rxn_core_carbon_cellD$v3)
rxn_core_carbon_cellD$id <- getSingleReactionFormula(unique_met_annotation$id,unique_met_annotation$name,rxn_core_carbon_cellD$v3)
rxn_core_carbon_cellD$MetaID <- getSingleReactionFormula(unique_met_annotation$MetaID,unique_met_annotation$name,rxn_core_carbon_cellD$v3)


rxn_core_carbon_cellD0 <- select(rxn_core_carbon_cellD, v2,v3,species,id, type,MetaID,note)
colnames(rxn_core_carbon_cellD0) <-c("rxnID","name","specie","id","type","MetaID","note")
write.table(rxn_core_carbon_cellD0, "reaction_core_carbon.txt", row.names = FALSE, sep="\t")




