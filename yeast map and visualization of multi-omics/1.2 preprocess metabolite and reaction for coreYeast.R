# This code is used to prepare the data for the map from model in BIGG format
source('model change.R')
source('transition for cellDesigner.R')

#part one
#prepare the reaction format
rxn <- read_excel("data/yeast_core_model.xlsx", sheet = "Reaction List")
metabolite <-  read_excel("data/yeast_core_model.xlsx", sheet = "Metabolite List")

# Unify the metabolite
metabolite <- select(metabolite, `Metabolite name`, `Metabolite formula`, Charge, KEGGID)
colnames(metabolite) <- c('Metabolite description', 'Metabolite formula', 'Charge', 'KEGGID')
#prepare the standard compartment
comparment <- unlist(str_extract_all(metabolite$`Metabolite description`, "_[:alpha:]$")) %>%
  str_replace_all(.,"_","[")
comparment <- paste(comparment, "]", sep = "")
for(i in seq_along(comparment)){
metabolite$`Metabolite description`[i] <- str_replace_all(metabolite$`Metabolite description`[i],"_[:alpha:]$", comparment[i])
}
metabolite$`Metabolite description` <- str_replace_all(metabolite$`Metabolite description`,"M_","")



rxn0 <- select(rxn, Abbreviation, Description_new)
colnames(rxn0) <- c('ID0', 'Equation')
rxn_split <- splitRxnToMetabolite(rxn0, sep0 = "<=>")
#first remove the exchange reaction
exchange_index <- which(is.na(rxn_split$MetID))
exchange_id <- rxn_split$ID[exchange_index]
others <- which(rxn_split$ID %in% exchange_id ==FALSE)
rxn_split_refine <- rxn_split[others,]
rxn_split_refine <- select(rxn_split_refine, ID, MetID, compostion)
colnames(rxn_split_refine) <- c('v2','v3','type')
rxn_split_refine$subsystem <- getSingleReactionFormula(rxn$Subsystem_new,rxn$Abbreviation,rxn_split_refine$v2)
rxn_split_refine$v2 <- str_replace_all(rxn_split_refine$v2,"R_","r_")

#prepare the standard compartment
comparment1 <- unlist(str_extract_all(rxn_split_refine$v3, "_[:alpha:]$")) %>%
  str_replace_all(.,"_","[")
comparment1 <- paste(comparment1, "]", sep = "")
for(i in seq_along(comparment1)){
  rxn_split_refine$v3[i] <- str_replace_all(rxn_split_refine$v3[i],"_[:alpha:]$", comparment1[i])
}
rxn_split_refine$v3 <- str_replace_all(rxn_split_refine$v3, "M_", "")



# choose the subsytem
subsystem1 <-  "sce00010  Glycolysis / Gluconeogenesis"

# Define the currency metabolite in each subsystem
currency_metabolites <- DefineCurrencyMet(rxn_split_refine, 
                                          subsystem0=subsystem1,
                                          numberGEM=14,
                                          numberSubsystem=1)

# remove the reactions with only one metabolite
# if we do not remove the currency metabolite in the model then this step is mainly removed exchange reaction
rxn_split_refine <- removeRxnWithSingleMet(rxn_split=rxn_split_refine)


#--------------------------------------------------------------------------------------------
## define the base reactant and product for cellDesigner
#---------------------------------------------------------------------------------------------
rxn_split_refine <- addBaseTypeIntoRxn(rxn_split_refine, metabolite, currency_metabolites)

#---------------------------------------------------
# choose reaction based on the subsystem
#-----------------------------------------------------
rxn_core_carbon <-chooseRxnFromSubsystem(rxn_split_refine, subsystem0=subsystem1) 
rxnID_choose <- unique(rxn_core_carbon$v2)

#------------------------------------------------------------------
# produce the met, rxn and gpr used for the map production
#------------------------------------------------------------------
# prepare the metabolites formula
# this funcion is used to prepare the metabolite annotation for cell designer
met_annotation <- prepareMET(rxn_core_carbon, currency_metabolites,rxnID_choose)

# prepare the rxn formula
rxn_core_carbon_cellD0 <- prepareRXN(rxn_core_carbon,met_annotation,currency_metabolites)
# prepare the protein and gene
gpr <- prepareGPR(met_annotation)

#save the exampel data format for cell designer
#write.table(met_annotation,"result/met_annotation for example.txt", row.names = FALSE, sep = "\n")
#write.table(rxn_core_carbon_cellD0,"result/rxn_core_carbon_cellD0 for example.txt", row.names = FALSE, sep = "\n")
#write.table(gpr,"result/gpr for example.txt", row.names = FALSE, sep = "\n")

#------------------------------------------------------------------
# produce the file as the input for the cellDesigner
#------------------------------------------------------------------
produceInputForCellDesigner(met_annotation, 
                            gpr,
                            rxn_core_carbon_cellD0,
                            x_size=1200, 
                            y_size=2000)
