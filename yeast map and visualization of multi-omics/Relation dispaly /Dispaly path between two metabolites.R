library(stringr)
library(tidyverse)
library(readxl)
library(igraph)
library(networkD3)
library(hongR)


#part one
#prepare the reaction format
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



## Define the currency metabolite in each subsystem
rxn_system <- select(rxn,Abbreviation,Subsystem_new)
metabolite_withoutCompartment <- rxn_split_refine
metabolite_withoutCompartment$name_simple <-  str_replace_all(metabolite_withoutCompartment$v3,"\\[.*?\\]","")
metabolite_withoutCompartment$subsystem <- getSingleReactionFormula(rxn_system$Subsystem_new,rxn_system$Abbreviation,metabolite_withoutCompartment$v2)


## define the general currency in whole model
analysis_metabolites <- metabolite_withoutCompartment %>% 
  count(name_simple) %>% ## calculate the number of each metabolite
  arrange(., desc(n)) ## order the metabolites based on the number
currency_metabolites_general <- analysis_metabolites$name_simple[1:14] ## find the currency metabolites



## classify the reactions based on subsystems
rxn_system <- select(rxn,Abbreviation,Subsystem_new)
index00 <- which(str_detect(rxn_system$Subsystem_new,"transport")==TRUE)
rxn_system$Subsystem_new[index00] <-"transport"
subsystem <- rxn_system %>% 
  count(Subsystem_new) %>%
  arrange(.,desc(n))

rxn_split_refine$subsystem <- getSingleReactionFormula(rxn_system$Subsystem_new,rxn_system$Abbreviation,rxn_split_refine$v2)




## define specific metabolites
## choosing the subsystem, the specific currency metabolite and reaction input will choosed on base on it
## choose several subsytem
#subsystem1 <- "glycolysis / gluconeogenesis \\( sce00010 \\)"
#subsystem2 <- "citrate cycle \\(tca cycle\\) \\( sce00020 \\)"
#subsystem3 <- "pentose phosphate pathway \\( sce00030 \\)"

#index1 <- which(str_detect(metabolite_withoutCompartment$subsystem, subsystem1)==TRUE)
#index2 <- which(str_detect(metabolite_withoutCompartment$subsystem, subsystem2)==TRUE)
#index3 <- which(str_detect(metabolite_withoutCompartment$subsystem, subsystem3)==TRUE)
#index_combine <- c(index1,index2,index3)


# choose one subsystem
subsystem1 <- "pyruvate metabolism \\( sce00620 \\)"
index_combine <- which(str_detect(metabolite_withoutCompartment$subsystem, subsystem1)==TRUE)

## define the currency in specific subsystem
metabolite_subsystem <- metabolite_withoutCompartment[index_combine,]
metabote_analysis_subsystem <- metabolite_subsystem %>%
  count(name_simple) %>%
  arrange(.,desc(n))
currency_metabolites_from_subsystem <- metabote_analysis_subsystem$name_simple[1]



#combine the general currency and specific currency
currency_metabolites <- unique(c(currency_metabolites_general,currency_metabolites_from_subsystem))

## remove currency metabolites
## based on Zhengming's suggestion, we will remain these currency metabolites in each subsystem
## but there will different methods to format for the currecy metabolite and general metabolite
## for these currency metabolite, there carbon number will be set as "" to avoid the wrong defination of base reactant or product

rxn_split_refine$simple_name <- str_replace_all(rxn_split_refine$v3,"\\[.*?\\]", "")
index_currency <- which(rxn_split_refine$simple_name %in% currency_metabolites == TRUE) ### it should be noted: currency_metabolites should be small range
rxn_split_refine$carbonNumber[index_currency] <- ""



#remove the reactions with only one metabolite
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
  ## test rxnID <- "r_0520"
  ## find the base reactant (with largest Carbon number)
  if(length(which(rxn_split_refine$v2 %in% rxnID ==TRUE & rxn_split_refine$type =="reactant")) !=0){
    rxn_reactant <- which(rxn_split_refine$v2 %in% rxnID ==TRUE & rxn_split_refine$type =="reactant")
    nn <- vector()
    if(all(rxn_split_refine$carbonNumber[rxn_reactant]=="")){ ##estimate whether it contains the carbon number
      nn[1] <- rxn_reactant[1]} else{
    ss <- max(as.numeric(rxn_split_refine$carbonNumber[rxn_reactant]),na.rm = TRUE)
    tt <- which(as.numeric(rxn_split_refine$carbonNumber[rxn_reactant]) == ss)
    
    nn[1] <- rxn_reactant[1]-1+tt
     } 
  } else{
    nn[1] <- ""
  }
  ## find the base product (with largest carbon nunber)
  if (length(which( rxn_split_refine$v2 %in% rxnID ==TRUE & rxn_split_refine$type =="product")) !=0){
    rxn_product <- which( rxn_split_refine$v2 %in% rxnID ==TRUE & rxn_split_refine$type =="product")
    
    if(all(rxn_split_refine$carbonNumber[rxn_product] =="")){
      nn[2] <-rxn_product[1]
    } else{
    ss0 <- max(as.numeric(rxn_split_refine$carbonNumber[rxn_product]),na.rm = TRUE)
    tt0 <- which(as.numeric(rxn_split_refine$carbonNumber[rxn_product]) == ss0)
    nn[2] <- rxn_product[1]-1+tt0
    
    }
    } else {
    nn[2] <-""
  }
  
  return(nn)
}

metabolite_type <- list()
rxn_ID <- unique(rxn_split_refine$v2)



## give the base definition for each reaction
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





## choose several subsytem
#subsystem1 <- "glycolysis / gluconeogenesis \\( sce00010 \\)"
#subsystem2 <- "citrate cycle \\(tca cycle\\) \\( sce00020 \\)"
#subsystem3 <- "pentose phosphate pathway \\( sce00030 \\)"

#index1 <- which(str_detect(rxn_split_refine$subsystem, subsystem1)==TRUE)
#index2 <- which(str_detect(rxn_split_refine$subsystem, subsystem2)==TRUE)
#index3 <- which(str_detect(rxn_split_refine$subsystem, subsystem3)==TRUE)
#index_combine0 <- c(index1,index2,index3)



# choose one subsystem
subsystem1 <- "pyruvate metabolism \\( sce00620 \\)"
index_combine0 <- which(str_detect(rxn_split_refine$subsystem, subsystem1)==TRUE)


###########choose reaction
rxn_core_carbon <- rxn_split_refine[index_combine0,]
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


## remove connected reactions which is composed of currency metabolites
rxn_id_subsytem <- unique(rxn_core_carbon$v2)
met_inRXNsubsytem <- list()
for (i in seq(length(rxn_id_subsytem))){
met_inRXNsubsytem[[i]] <- rxn_core_carbon$simple_name[which(rxn_core_carbon$v2 %in% rxn_id_subsytem[i] ==TRUE)]
}

rxn_choose <- vector()
for (i in seq(length(rxn_id_subsytem))){
  if(all(met_inRXNsubsytem[[i]] %in% currency_metabolites_general ==TRUE)){
    rxn_choose[i] <- FALSE
  } else{
    rxn_choose[i] <- TRUE
  }
  }

rxnID_choose <- rxn_id_subsytem[rxn_choose]
rxn_core_carbon <- rxn_core_carbon[which(rxn_core_carbon$v2 %in% rxnID_choose ==TRUE),]




#part two
#prepare the metabolites formula
met_core_carbon <- select(rxn_core_carbon,v2,v3,type,simple_name)
colnames(met_core_carbon) <- c("rxnID","name","type", "simple_name")
##replace the type
met_core_carbon$type <- "SIMPLE_MOLECULE"

#divide the met into two types
met_no_currency <- met_core_carbon[which(met_core_carbon$simple_name %in% currency_metabolites ==FALSE),]

# remove the duplicated metabolite in met which is not currency
met_no_currency0 <- met_no_currency[!duplicated(met_no_currency$name),]
met_currency <- met_core_carbon[which(met_core_carbon$simple_name %in% currency_metabolites == TRUE),]

protein_annotation <- data.frame(rxnID=character(length(rxnID_choose)),stringsAsFactors = FALSE)
gene_annotation <- data.frame(rxnID=character(length(rxnID_choose)),stringsAsFactors = FALSE)
protein_annotation$rxnID <- rxnID_choose
protein_annotation$name <- str_replace_all(protein_annotation$rxnID, "r", "p")
protein_annotation$type <- "PROTEIN"

gene_annotation$rxnID <- rxnID_choose
gene_annotation$name <- str_replace_all(gene_annotation$rxnID, "r", "g")
gene_annotation$type <- "GENE"

met_final <- rbind.data.frame(select(met_no_currency0,"rxnID","name","type"), select(met_currency,"rxnID","name","type"), protein_annotation, gene_annotation)

## sort by rxnID
met_final <- met_final[order(met_final$rxnID),]
met_final0 <- met_final$name

## define the unique species
unique_species <- data.frame(species=character(length(unique(met_final0))),stringsAsFactors = FALSE)
unique_species$name <- unique(met_final0)
unique_species$species <- paste("s", 1:length(unique(met_final0)), sep="") 
met_annotation <- data.frame(species=character(length(met_final0)),stringsAsFactors = FALSE)
met_annotation$name <- met_final0 
met_annotation$species <- getSingleReactionFormula(unique_species$species,unique_species$name,met_annotation$name)
met_annotation$id <-paste("sa", 1:length(met_final0), sep="") 
met_annotation$MetaID <- paste("CDMT0000", 1:length(met_final0), sep="") 

#define the size of whole graph
x_vector <- seq(100, 12000, by=100)
y_vector <- seq(100, 15000, by=50)
met_annotation$x <-rep(x_vector, each = 30)[1:length(met_annotation$id)]
met_annotation$y <-rep(y_vector, times = 12)[1:length(met_annotation$id)]

met_annotation$type <- met_final$type
met_annotation$rxnID <- met_final$rxnID
met_annotation$metID <- paste("m", 1:length(met_final0), sep="")
met_annotation <- select(met_annotation,metID,id,species,name,x,y,type,MetaID,rxnID) # will be the data source for import metabolites into graph



#prepare the rxn formula
rxn_core_carbon_cellD <- select(rxn_core_carbon, v2,v3,type,note,simple_name)
colnames(rxn_core_carbon_cellD) <- c("rxnID","name","type","note","simple_name")
# this metabolite in rxn should be divided into two types
rxn_core_with_currency <- rxn_core_carbon_cellD[which(rxn_core_carbon_cellD$simple_name %in% currency_metabolites ==TRUE),]
rxn_core_with_currency$combine_rxnID_name <- paste(rxn_core_with_currency$rxnID,rxn_core_with_currency$name)

# give the information of currency metabolites in rxn by mapping rhe "combine_rxnID_name"
met_annotation2 <- met_annotation
met_annotation2$combine_rxnID_name <- paste(met_annotation2$rxnID,met_annotation2$name)
rxn_core_with_currency$specie <- getSingleReactionFormula(met_annotation2$species,met_annotation2$combine_rxnID_name,rxn_core_with_currency$combine_rxnID_name)
rxn_core_with_currency$id <- getSingleReactionFormula(met_annotation2$id,met_annotation2$combine_rxnID_name,rxn_core_with_currency$combine_rxnID_name)
rxn_core_with_currency$MetaID <- getSingleReactionFormula(met_annotation2$MetaID,met_annotation2$combine_rxnID_name,rxn_core_with_currency$combine_rxnID_name)

# give the information of general metabolites in rxn by mapping the name of metabolite
rxn_core_without_currency <- rxn_core_carbon_cellD[which(rxn_core_carbon_cellD$simple_name %in% currency_metabolites ==FALSE),]
rxn_core_without_currency$specie <- getSingleReactionFormula(met_annotation$species,met_annotation$name,rxn_core_without_currency$name)
rxn_core_without_currency$id <- getSingleReactionFormula(met_annotation$id,met_annotation$name,rxn_core_without_currency$name)
rxn_core_without_currency$MetaID <- getSingleReactionFormula(met_annotation$MetaID,met_annotation$name,rxn_core_without_currency$name)

names_column <- c("rxnID","name","specie","id","type","MetaID","note")
rxn_core_carbon_cellD0 <- rbind.data.frame(rxn_core_without_currency[,names_column], rxn_core_with_currency[,names_column])
rxn_core_carbon_cellD0 <- rxn_core_carbon_cellD0[order(rxn_core_carbon_cellD0$rxnID),] # will be data source to import the reactions

#prepare the protein and gene
#define the proteinID information
protein_annotation0 <- filter(met_annotation, type=="PROTEIN") %>%
  select(.,rxnID,species,id,MetaID)
colnames(protein_annotation0) <- c("rxnID","protein_specie","protein_id","MetaID_p")

#define the the gene information
gene_annotation0 <- filter(met_annotation, type=="GENE") %>%
  select(.,rxnID,species,id,MetaID)
colnames(gene_annotation0) <- c("rxnID","gene_specie","gene_id","MetaID_g")

#define the final gene-protein-reaction information
gpr <- merge.data.frame(protein_annotation0,gene_annotation0) # will be the data source to import the reactions





