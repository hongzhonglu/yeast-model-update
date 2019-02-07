# This code is used to prepare the data for the map
source('model change.R')
source('transition for cellDesigner.R')

# prepare the reaction format
rxn <- read_excel("data/yeastGEM_latest version.xls", sheet = "Reaction List")
metabolite <-  read_excel("data/yeastGEM_latest version.xls", sheet = "Metabolite List")
rxn_split_refine <- splitRxnToMetabolite.Yeast(rxn, metabolite)


# choose the subsytem
subsystem1 <- "glycolysis / gluconeogenesis \\( sce00010 \\)"

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

