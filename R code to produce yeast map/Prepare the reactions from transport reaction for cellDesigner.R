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


which(str_detect(rxn_split$v3,"]")==TRUE)
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


## define the base reactant and product for cellDesigner
## it should be noted that the currency metabolites is not considered in this phase
rxn_split_refine$note <- ""
rownames(rxn_split_refine) <- 1:12782

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



##when run the above code, get the followed error information
##Error in nn[2] <- rxn_product[1] - 1 + tt0 : replacement has length zero
##In addition: There were 50 or more warnings (use warnings() to see the first 50)



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

rxn_transport <- filter(rxn_split_refine,subsystem =="transport")





#part two
#prepare the metabolites formula

met_transport <- select(rxn_transport,v2,v3,type)
unique_met_transport <- unique(met_transport$v3)

unique_met_annotation <- data.frame(species=character(length(unique_met_transport)),stringsAsFactors = FALSE)
unique_met_annotation$name <- unique_met_transport 
unique_met_annotation$species <- paste("s", 1:1311, sep="")

colnames(met_transport) <- c("rxnID","name","type")
met_transport$species <- getSingleReactionFormula(unique_met_annotation$species,unique_met_annotation$name,met_transport$name)
met_transport$id <-paste("sa", 1:2100, sep="") 
met_transport$MetaID <- paste("CDMT0000", 1:2100, sep="") 
met_transport$x <-""
met_transport$y <-""
met_transport <- select(met_transport,id,species,name,x,y,rxnID,type,MetaID)
write.table(met_transport, "met_transport.txt", row.names = FALSE, sep="\t")



#prepare the rxn formula
rxn_transport_cellD <- select(rxn_transport, v2,v3,type,note)
rxn_transport_cellD$species <- getSingleReactionFormula(unique_met_annotation$species,unique_met_annotation$name,rxn_transport_cellD$v3)
rxn_transport_cellD$id <- met_transport$id #should be careful
rxn_transport_cellD$MetaID <- met_transport$MetaID #should be careful


rxn_transport_cellD0 <- select(rxn_transport_cellD, v2,v3,species,id, type,MetaID,note)
colnames(rxn_transport_cellD0) <-c("rxnID","name","specie","id","type","MetaID","note")
write.table(rxn_transport_cellD0, "reaction_transport.txt", row.names = FALSE, sep="\t")




