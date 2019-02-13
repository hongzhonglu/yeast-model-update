library(stringr)
library(tidyverse)
library(readxl)
library(igraph)
library(networkD3)
library(hongR)

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


#this function is used to split the reaction
splitRxnToMetabolite <- function(reationFrame, sep0){
  #reationFrame <- rxn_final0
  #input:
  #dataframe: column 1--ID0; column 2--Equation, which is the reaction formula
  reationFrame$Equation <- str_replace_all(reationFrame$Equation, "->", "<=>")
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
  rxn_met$MetID <- str_trim(rxn_met$MetID, side = "both")
  for (i in seq_along(rxn_met$ID)){
    if(str_detect(rxn_met$MetID[i], "^M_")){
      rxn_met$MetID[i] <- paste(1, rxn_met$MetID[i], " ")
    }
  }
  
  #remove the coefficient for each metabolites, which will be easy to standardize this metabolite
  rxn_met <- rxn_met %>% separate(.,MetID, into = c('coefficient','MetID'), sep = " ")
  #rxn_met$MetID <- str_replace_all(rxn_met$MetID, "[:digit:] ","") %>%
  #   str_replace_all(.,"\\(n\\)","") %>%
  #   str_replace_all(.,"\\(n\\+1\\)","") %>%
  #   str_replace_all(.,"\\(n\\-2\\)","")
  # the followed two step is used to process data from mnx database
  rxn_met$MetID <- str_replace_all(rxn_met$MetID, "@MNXD[:alnum:]","")
  rxn_met$MetID <- str_replace_all(rxn_met$MetID, "@BOUNDARY","")
  return(rxn_met)
}

#this function is used to split the reaction
splitRxnToMetabolite.Ecoli <- function(reationFrame, sep0){
  #reationFrame <- rxn_final0
  #input:
  #dataframe: column 1--ID0; column 2--Equation, which is the reaction formula
  reationFrame$Equation <- str_replace_all(reationFrame$Equation, "->", "<=>")
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
  rxn_met$MetID <- str_trim(rxn_met$MetID, side = "both")
  for (i in seq_along(rxn_met$ID)){
    if(!str_detect(rxn_met$MetID[i], "^\\d+\\.*\\d* ")){
      rxn_met$MetID[i] <- paste(1, rxn_met$MetID[i], " ")
    }
  }
  #remove the coefficient for each metabolites, which will be easy to standardize this metabolite
  rxn_met <- rxn_met %>% separate(.,MetID, into = c('coefficient','MetID'), sep = " ")
  #rxn_met$MetID <- str_replace_all(rxn_met$MetID, "[:digit:] ","") %>%
  #   str_replace_all(.,"\\(n\\)","") %>%
  #   str_replace_all(.,"\\(n\\+1\\)","") %>%
  #   str_replace_all(.,"\\(n\\-2\\)","")
  # the followed two step is used to process data from mnx database
  rxn_met$MetID <- str_replace_all(rxn_met$MetID, "@MNXD[:alnum:]","")
  rxn_met$MetID <- str_replace_all(rxn_met$MetID, "@BOUNDARY","")
  return(rxn_met)
}

# this function is designed for yeast GEM to split the reaction
splitRxnToMetabolite.Yeast <- function(rxn_inf, metabolite_inf) {
  #input
  # rxn_inf, rxn annotation of yeast GEM
  # metabolite_inf, metabolite annotation of yeast GEM
  #output
  # rxn_split_refine a dataframe contained the split format of yeast GEM
  
  rxn_split <- splitAndCombine0(rxn_inf$Reaction, rxn_inf$Abbreviation, sep = " ")
  rxn_split$v3 <- getSingleReactionFormula(metabolite_inf$`Metabolite description`, metabolite_inf$`Metabolite name`, rxn_split$v1)
  
  for (i in 1:length(rxn_split$v2)) {
    if (rxn_split$v3[i] == "NA") {
      rxn_split$v3[i] <- rxn_split$v1[i]
    } else {
      rxn_split$v3[i] <- rxn_split$v3[i]
    }
  }
  # give the type: reactant and product for each metabolite
  # for r_0001
  rxn_split$v3 <- str_replace_all(rxn_split$v3, "->", "<=>")
  rxn_ID <- unique(rxn_split$v2)
  met_type0 <- vector()
  for (i in 1:length(rxn_ID)) {
    met_type0 <- c(met_type0, getMetType(rxn_ID[i], rxn_split))
  }
  
  rxn_split$type <- met_type0
  rxn_split_refine <- rxn_split[which(str_detect(rxn_split$v3, "]") == TRUE), ]
  
  ## classify the reactions based on subsystems
  rxn_system <- select(rxn_inf, Abbreviation, Subsystem_new)
  index00 <- which(str_detect(rxn_system$Subsystem_new, "transport") == TRUE)
  rxn_system$Subsystem_new[index00] <- "transport"
  subsystem <- rxn_system %>%
    count(Subsystem_new) %>%
    arrange(., desc(n))
  
  rxn_split_refine$subsystem <- getSingleReactionFormula(rxn_system$Subsystem_new, rxn_system$Abbreviation, rxn_split_refine$v2)
  return(rxn_split_refine)
  
}

#when removing the currency metabolites. some reaction could cotain only metabolite
#this function is used to remove the exchange reaction and reactions with only one reactant or product
removeRxnWithSingleMet <- function(rxn_split=rxn_split_refine) {
  #input
  #rxn_split, a split rxn format
  #output
  #rxn_split, a split rxn in which a reaction with single metabolite or only with reactant (product) was removed
  
  analysis_rxn <- rxn_split %>%
    count(v2) %>% ## calculate the number of each metabolite
    arrange(., desc(n)) ## order the metabolites based on the number
  
  rxn_with_one_met <- which(analysis_rxn$n == 1)
  rxn_remove <- analysis_rxn$v2[rxn_with_one_met]
  which(rxn_split$v2 %in% rxn_remove == FALSE)
  rxn_split <- rxn_split[which(rxn_split$v2 %in% rxn_remove == FALSE), ]
  
  # further remove the reactions with only one product or only one reactant
  rxn_index <- list()
  rxn_unique <- unique(rxn_split$v2)
  for (i in 1:length(rxn_unique)) {
    rxn_index[[i]] <- which(rxn_split$v2 %in% rxn_unique[i])
  }
  
  rxn_metabotite_type <- list()
  for (i in 1:length(rxn_unique)) {
    rxn_metabotite_type[[i]] <- rxn_split$type[rxn_index[[i]]]
  }
  
  ss <- vector()
  for (i in 1:length(rxn_unique)) {
    if ("reactant" %in% rxn_metabotite_type[[i]] & "product" %in% rxn_metabotite_type[[i]]) {
      ss[i] <- i
    } else {
      ss[i] <- NA
    }
  }
  rxn_final <- rxn_unique[!is.na(ss)]
  rxn_split <- rxn_split[which(rxn_split$v2 %in% rxn_final == TRUE), ]
  return(rxn_split)
}


# this function is used to estimate the metabolite type in a reaction.
getMetType <- function(rxnID_inf, rxn_split_inf){
  # input  rxnID_inf a list of reaction id
  # input  rxn_split_inf a dataframe with three columns
  # output met_type, a list contains the type for all the metabolites in a reaction
  ss <- which(rxn_split_inf$v2 %in% rxnID_inf==TRUE)
  ss_split <- which(rxn_split_inf$v2 %in% rxnID_inf ==TRUE & rxn_split_inf$v3 %in% "<=>" == TRUE)
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

# the fuction is used to extract the carbon part for a metabolite
getMetaboliteComposition <- function(MET, type = "C") {
  #input a metabolite formula, like C5H6
  #output the carbon part of a metabolite, like C5
  if (str_detect(MET, paste(type, "[0-9]+", sep = ""))) {
    carbon_number <- str_extract_all(MET, paste(type, "[0-9]+", sep = ""))
  } else if (str_detect(MET, "C") & str_detect(MET, paste(type, "[0-9]+", sep = "")) == FALSE) {
    carbon_number <- str_extract_all(MET, "C")
  } else {
    carbon_number <- ""
  }
  return(carbon_number[[1]][1])
}

# the function is used to extract the number of carbon based on the above function
getCompositionNum <- function(MET,type="C"){
  # input a carbon part of a molecular, like C5
  # output the number of carbon, like 5
  if (str_detect(MET, paste(type,"[0-9]+", sep=""))){
    carbon_number<- str_extract_all(MET, "[0-9]+")
  } else if (str_detect(MET, "C") & str_detect(MET, paste(type,"[0-9]+", sep=""))==FALSE ){
    carbon_number <- "1"
  } else{
    carbon_number <- ""
  }
  return(carbon_number[[1]][1])
}

# This function is used to define the currency of metabolites in the model
DefineCurrencyMet <- function(rxn_split_refine0, subsystem0,numberGEM, numberSubsystem) {
  #input
  #rxn_split_refine0, a split rxn format with detailed annotation with 10 columns
  #subsystem0, chose subsytem
  #numberGEM,a number used to limit the number of currency metabolite from the whole model
  #numberSubsystem, a number used to limit the number of currency metabolite from each chose subsystems
  #output
  #currency_metabolites, a vector of currency metabolite
  
  rxn_system0 <- select(rxn_split_refine0, v2, subsystem)
  colnames(rxn_system0) <-  c('Abbreviation','subsystem')
  rxn_system0 <- rxn_system0[!duplicated(rxn_system0$Abbreviation),]
  metabolite_withoutCompartment <- rxn_split_refine0
  metabolite_withoutCompartment$name_simple <- str_replace_all(metabolite_withoutCompartment$v3, "\\[.*?\\]", "")
  metabolite_withoutCompartment$subsystem <- getSingleReactionFormula(rxn_system0$subsystem, rxn_system0$Abbreviation, metabolite_withoutCompartment$v2)
  
  ## define the general currency in whole model
  analysis_metabolites <- metabolite_withoutCompartment %>%
    count(name_simple) %>% ## calculate the number of each metabolite
    arrange(., desc(n)) ## order the metabolites based on the number
  currency_metabolites_general <- analysis_metabolites$name_simple[1:numberGEM] ## find the currency metabolites
  
  cat('statistical analysis of metabolites in model')
  print(analysis_metabolites)
  
  # choose one subsystem
  # subsystem0 <- "pyruvate metabolism \\( sce00620 \\)"
  index_combine <- which(str_detect(metabolite_withoutCompartment$subsystem, subsystem0) == TRUE)
  
  ## define the currency in specific subsystem
  metabolite_subsystem <- metabolite_withoutCompartment[index_combine, ]
  metabote_analysis_subsystem <- metabolite_subsystem %>%
    count(name_simple) %>%
    arrange(., desc(n))
  cat('statistical analysis of metabolites in subsystem')
  print(metabote_analysis_subsystem)
  currency_metabolites_from_subsystem <- metabote_analysis_subsystem$name_simple[1:numberSubsystem]
  # combine the general currency and specific currency
  currency_metabolites <- unique(c(currency_metabolites_general, currency_metabolites_from_subsystem))
  cat('Choose currency metabolite for map of subsystem \n')
  print(currency_metabolites)
  return(currency_metabolites[!is.na(currency_metabolites)])
}

## define the base reactant and product for cellDesigner
DefineBaseMetabolite <- function(rxnID0, rxn_split_refine_inf){
  ## input: rxnID0: a reaction ID
  ## input: rxn_split_refine_inf: datafram contains the relation of rxnID0 and metabolite, as well as metabolite formula and 
  ## and carbon number
  ## test rxnID0 <- "r_0520"
  ## find the base reactant (with largest Carbon number)
  if(length(which(rxn_split_refine_inf$v2 %in% rxnID0 ==TRUE & rxn_split_refine_inf$type =="reactant")) !=0){
    rxn_reactant <- which(rxn_split_refine_inf$v2 %in% rxnID0 ==TRUE & rxn_split_refine_inf$type =="reactant")
    nn <- vector()
    if(all(rxn_split_refine_inf$carbonNumber[rxn_reactant]=="")){ ##estimate whether it contains the carbon number
      nn[1] <- rxn_reactant[1]} else{
        ss <- max(as.numeric(rxn_split_refine_inf$carbonNumber[rxn_reactant]),na.rm = TRUE)
        tt <- which(as.numeric(rxn_split_refine_inf$carbonNumber[rxn_reactant]) == ss)
        
        nn[1] <- rxn_reactant[1]-1+tt
      } 
  } else{
    nn[1] <- ""
  }
  ## find the base product (with largest carbon nunber)
  if (length(which( rxn_split_refine_inf$v2 %in% rxnID0 ==TRUE & rxn_split_refine_inf$type =="product")) !=0){
    rxn_product <- which( rxn_split_refine_inf$v2 %in% rxnID0 ==TRUE & rxn_split_refine_inf$type =="product")
    
    if(all(rxn_split_refine_inf$carbonNumber[rxn_product] =="")){
      nn[2] <-rxn_product[1]
    } else{
      ss0 <- max(as.numeric(rxn_split_refine_inf$carbonNumber[rxn_product]),na.rm = TRUE)
      tt0 <- which(as.numeric(rxn_split_refine_inf$carbonNumber[rxn_product]) == ss0)
      nn[2] <- rxn_product[1]-1+tt0
      
    }
  } else {
    nn[2] <-""
  }
  
  return(nn)
}


#this function is used to give the "base" sign for two of metabolites from each reaction
addBaseTypeIntoRxn <- function(rxn_split_refine_inf, metabolite_inf, currency_metabolites_inf){
  ## remove currency metabolites
  ## based on Zhengming's suggestion, we will remain these currency metabolites in each subsystem
  ## but there will different methods to format for the currecy metabolite and general metabolite
  ## for these currency metabolite, there carbon number will be set as "" to avoid the wrong defination of base reactant or product
  
  #input
  # rxn_split_refine_inf  rxn split format
  # metabolite_inf metabolite annotation information
  # currency_metabolites_inf a vector of currency metabolites
  #output
  #rxn_split_refine, rxn_split with "base" sign
  
  
  
  ## obtain other information
  rxn_split_refine_inf$metabolit_formula <- getSingleReactionFormula(metabolite_inf$`Metabolite formula`,metabolite_inf$`Metabolite description`,rxn_split_refine_inf$v3)
  
  for (i in 1:length(rxn_split_refine_inf$v2)){
    rxn_split_refine_inf$carbonCompostion[i] <- getMetaboliteComposition(rxn_split_refine_inf$metabolit_formula[i])
  }
  
  for (i in 1:length(rxn_split_refine_inf$v2)){
    rxn_split_refine_inf$carbonNumber[i] <- getCompositionNum(rxn_split_refine_inf$carbonCompostion[i])
  }
  rxn_split_refine_inf$simple_name <- str_replace_all(rxn_split_refine_inf$v3,"\\[.*?\\]", "")
  index_currency <- which(rxn_split_refine_inf$simple_name %in% currency_metabolites_inf == TRUE) ### it should be noted: currency_metabolites should be small range
  rxn_split_refine_inf$carbonNumber[index_currency] <- ""
  
  rxn_split_refine_inf$note <- ""
  rownames(rxn_split_refine_inf) <- 1:nrow(rxn_split_refine_inf)
  metabolite_type <- list()
  rxn_ID <- unique(rxn_split_refine_inf$v2)
  ## give the base definition for each reaction
  for (i in 1:length(rxn_ID)) {
    metabolite_type[[i]] <- DefineBaseMetabolite(rxn_ID[i], rxn_split_refine_inf)
  }
  
  metabolite_type <- unlist(metabolite_type)
  for (i in 1:length(rxn_split_refine_inf$v2)) {
    if (i %in% metabolite_type) {
      rxn_split_refine_inf$note[i] <- "base"
    } else {
      rxn_split_refine_inf$note[i] <- ""
    }
  }
  
  return(rxn_split_refine_inf)
}


# this function is used to get the id for the transport reaction which could connect the metabolites occured in different compartment
# for specific subsystem
getConnectedTransport <- function (id,rxn_transport_id0, rxn_transport0, met_core_carbon0){
  #input
  # id, a index of trabsport reaction id
  # rxn_transport_id0, a vector of transport id
  # rxn_transport0, a dataframe contains the annotation information for transport reaction
  # met_core_carbon0, a vector of the metabolite list from specific subsystem
  # output
  # mm, the id of transport reactions which choosed
  met_transport_index <- which(rxn_transport0$v2 %in% rxn_transport_id0[id] )
  met_transport <- rxn_transport0$v3[met_transport_index]
  ss <- vector()
  ss <- met_transport %in% met_core_carbon0
  tt <- sum(ss)
  # estimate whether the metabolite in a transport reaction all occured in a chose metabolite list
  if (tt == length(met_transport)){
    mm <- id
  } else{
    mm <- NA
    
  }
  return(mm)
}


# this function is used to choose the reaction based on subsytem definition
chooseRxnFromSubsystem <- function(rxn_split_refine_inf, subsystem0) {
  #input
  #rxn_split_refine_inf, rxn split format with the detailed annotation information
  #output
  #subsystem0, a string of defined subsystem
  
  index_combine0 <- which(str_detect(rxn_split_refine_inf$subsystem, subsystem0) == TRUE)
  
  ########### choose reaction
  rxn_core_carbon <- rxn_split_refine_inf[index_combine0, ]
  met_core_carbon <- unique(rxn_core_carbon$v3)
  
  # find the transport reactions to connect the gap in the above systerm
  index_transport <- which(str_detect(rxn_split_refine_inf$subsystem, "transport") == TRUE)
  rxn_transport <- rxn_split_refine_inf[index_transport, ]
  rxn_transport_id <- unique(rxn_transport$v2)
  connected_rxn <- vector()
  for (i in 1:length(rxn_transport_id)) {
    connected_rxn[i] <- getConnectedTransport(id = i, rxn_transport_id, rxn_transport, met_core_carbon)
  }
  connect_rxn0 <- connected_rxn[!is.na(connected_rxn)]
  ## add the connected transport reactions
  connect_rxn0
  trasport_choosed <- rxn_transport_id[connect_rxn0]
  trasport_rxn_choosed <- rxn_transport[which(rxn_transport$v2 %in% trasport_choosed == TRUE), ]
  rxn_core_carbon <- rbind.data.frame(rxn_core_carbon, trasport_rxn_choosed)
  
  
  ## remove connected reactions which is composed of currency metabolites
  rxn_id_subsytem <- unique(rxn_core_carbon$v2)
  met_inRXNsubsytem <- list()
  for (i in seq(length(rxn_id_subsytem))) {
    met_inRXNsubsytem[[i]] <- rxn_core_carbon$simple_name[which(rxn_core_carbon$v2 %in% rxn_id_subsytem[i] == TRUE)]
  }
  
  rxn_choose <- vector()
  for (i in seq(length(rxn_id_subsytem))) {
    if (all(met_inRXNsubsytem[[i]] %in% currency_metabolites == TRUE)) {
      rxn_choose[i] <- FALSE
    } else {
      rxn_choose[i] <- TRUE
    }
  }
  
  rxnID_choose0 <- rxn_id_subsytem[rxn_choose]
  rxn_core_carbon <- rxn_core_carbon[which(rxn_core_carbon$v2 %in% rxnID_choose0 == TRUE), ]
  return(rxn_core_carbon)
}



#this function is used to define the coordinate information for the metabolites
prepareMET <- function(rxn_core_carbon_inf, 
                       currency_metabolites_inf,
                       rxnID_choose_inf,  
                       x_vector = seq(100, 12000, by = 100), 
                       y_vector = seq(100, 15000, by = 50)) {
  
  #input
  # rxn_core_carbon_inf, a dataframe contained the detailed annotation of rxn_core
  # currency_metabolites_inf, a vector of currency metabolites
  # rxnID_choose_inf, a vector of choosed rxnID
  # x_vector, size of graph
  # y_vector, size of graph
  #outout
  #met_annotation, a dataframe contains the detailed annotation for met
  
  
  met_core_carbon <- select(rxn_core_carbon_inf, v2, v3, type, simple_name)
  colnames(met_core_carbon) <- c("rxnID", "name", "type", "simple_name")
  ## replace the type
  met_core_carbon$type <- "SIMPLE_MOLECULE"
  
  # divide the met into two types
  met_no_currency <- met_core_carbon[which(met_core_carbon$simple_name %in% currency_metabolites_inf == FALSE), ]
  
  # remove the duplicated metabolite in met which is not currency
  met_no_currency0 <- met_no_currency[!duplicated(met_no_currency$name), ]
  met_currency <- met_core_carbon[which(met_core_carbon$simple_name %in% currency_metabolites_inf == TRUE), ]
  
  protein_annotation <- data.frame(rxnID = character(length(rxnID_choose_inf)), stringsAsFactors = FALSE)
  gene_annotation <- data.frame(rxnID = character(length(rxnID_choose_inf)), stringsAsFactors = FALSE)
  protein_annotation$rxnID <- rxnID_choose_inf
  protein_annotation$name <- str_replace_all(protein_annotation$rxnID, "r", "p")
  protein_annotation$type <- "PROTEIN"
  
  gene_annotation$rxnID <- rxnID_choose_inf
  gene_annotation$name <- str_replace_all(gene_annotation$rxnID, "r", "g")
  gene_annotation$type <- "GENE"
  
  met_final <- rbind.data.frame(select(met_no_currency0, "rxnID", "name", "type"), select(met_currency, "rxnID", "name", "type"), protein_annotation, gene_annotation)
  
  ## sort by rxnID
  met_final <- met_final[order(met_final$rxnID), ]
  met_final0 <- met_final$name
  
  ## define the unique species
  unique_species <- data.frame(species = character(length(unique(met_final0))), stringsAsFactors = FALSE)
  unique_species$name <- unique(met_final0)
  unique_species$species <- paste("s", 1:length(unique(met_final0)), sep = "")
  met_annotation <- data.frame(species = character(length(met_final0)), stringsAsFactors = FALSE)
  met_annotation$name <- met_final0
  met_annotation$species <- getSingleReactionFormula(unique_species$species, unique_species$name, met_annotation$name)
  met_annotation$id <- paste("sa", 1:length(met_final0), sep = "")
  met_annotation$MetaID <- paste("CDMT0000", 1:length(met_final0), sep = "")
  
  # define the size of whole graph
  met_annotation$x <- rep(x_vector, each = 30)[1:length(met_annotation$id)]
  met_annotation$y <- rep(y_vector, times = 12)[1:length(met_annotation$id)]
  
  met_annotation$type <- met_final$type
  met_annotation$rxnID <- met_final$rxnID
  met_annotation$metID <- paste("m", 1:length(met_final0), sep = "")
  met_annotation <- select(met_annotation, metID, id, species, name, x, y, type, MetaID, rxnID) # will be the data source for import metabolites into graph
  
  return(met_annotation)
}



#this function is used to extract the gene and protein annotation from met_annotation
prepareGPR <- function(met_annotation_inf) {
  # input
  # met_annotation_inf a dataframe for the metabolite annotation with 9 columns which contains the metabolite, gene and proteins
  # output
  # gpr, a dataframe contains the annotation information for gene and protein
  
  # define the proteinID information
  protein_annotation0 <- filter(met_annotation_inf, type == "PROTEIN") %>%
    select(., rxnID, species, id, MetaID)
  colnames(protein_annotation0) <- c("rxnID", "protein_specie", "protein_id", "MetaID_p")
  
  # define the the gene information
  gene_annotation0 <- filter(met_annotation_inf, type == "GENE") %>%
    select(., rxnID, species, id, MetaID)
  colnames(gene_annotation0) <- c("rxnID", "gene_specie", "gene_id", "MetaID_g")
  # define the final gene-protein-reaction information
  gpr <- merge.data.frame(protein_annotation0, gene_annotation0) # will be the data source to import the reactions
  return(gpr)
}



# this function is to produce the rxn format used for celldesigner
prepareRXN <- function(rxn_core_carbon_inf, met_annotation_inf, currency_metabolites_inf) {
  #input
  # rxn_core_carbon_inf, a dataframe contained the detailed annotation of rxn_core
  # met_annotation_inf, a dataframe contained the detailed annotation of met
  # currency_metabolites_inf, a vector of currency metabolites
  # output
  # rxn_core_carbon_cellD0, a dataframe for the rxn annotation with 7 columns
  
  rxn_core_carbon_cellD <- select(rxn_core_carbon_inf, v2, v3, type, note, simple_name)
  colnames(rxn_core_carbon_cellD) <- c("rxnID", "name", "type", "note", "simple_name")
  # this metabolite in rxn should be divided into two types
  rxn_core_with_currency <- rxn_core_carbon_cellD[which(rxn_core_carbon_cellD$simple_name %in% currency_metabolites_inf == TRUE), ]
  rxn_core_with_currency$combine_rxnID_name <- paste(rxn_core_with_currency$rxnID, rxn_core_with_currency$name)
  
  # give the information of currency metabolites in rxn by mapping rhe "combine_rxnID_name"
  met_annotation2 <- met_annotation_inf
  met_annotation2$combine_rxnID_name <- paste(met_annotation2$rxnID, met_annotation2$name)
  rxn_core_with_currency$specie <- getSingleReactionFormula(met_annotation2$species, met_annotation2$combine_rxnID_name, rxn_core_with_currency$combine_rxnID_name)
  rxn_core_with_currency$id <- getSingleReactionFormula(met_annotation2$id, met_annotation2$combine_rxnID_name, rxn_core_with_currency$combine_rxnID_name)
  rxn_core_with_currency$MetaID <- getSingleReactionFormula(met_annotation2$MetaID, met_annotation2$combine_rxnID_name, rxn_core_with_currency$combine_rxnID_name)
  
  # give the information of general metabolites in rxn by mapping the name of metabolite
  rxn_core_without_currency <- rxn_core_carbon_cellD[which(rxn_core_carbon_cellD$simple_name %in% currency_metabolites_inf == FALSE), ]
  rxn_core_without_currency$specie <- getSingleReactionFormula(met_annotation_inf$species, met_annotation_inf$name, rxn_core_without_currency$name)
  rxn_core_without_currency$id <- getSingleReactionFormula(met_annotation_inf$id, met_annotation_inf$name, rxn_core_without_currency$name)
  rxn_core_without_currency$MetaID <- getSingleReactionFormula(met_annotation_inf$MetaID, met_annotation_inf$name, rxn_core_without_currency$name)
  
  names_column <- c("rxnID", "name", "specie", "id", "type", "MetaID", "note")
  rxn_core_carbon_cellD0 <- rbind.data.frame(rxn_core_without_currency[, names_column], rxn_core_with_currency[, names_column])
  rxn_core_carbon_cellD0 <- rxn_core_carbon_cellD0[order(rxn_core_carbon_cellD0$rxnID), ] # will be data source to import the reactions
  return(rxn_core_carbon_cellD0)
}




