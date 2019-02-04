### code to update the metabolite name in the model, which will be used for the yeast map
source('main_function_map.R')


# input the model and metabolite
# the model was obtained using the followed matlab code
# model = readCbModel('yeastGEM')
# outmodel = writeCbModel(model, 'format','xls', 'fileName', 'yeastGEM_october.xls')

rxn <- read_excel("data/yeastGEM_october.xls", sheet = "Reaction List")
metabolite <-  read_excel("data/yeastGEM_october.xls", sheet = "Metabolite List")


# to kegg the consistency in metabolite name, we use the old name of metablite from yeast7.6
# input the metabolite annotation from yeast7.6
metabolite_v7_6 <-  read_excel("data/yeastGEM_latest version.xls", 
                                                         sheet = "Metabolite List")
#first comparsion between the metabolites in different version
metabolite$Description_v_7_6<- getSingleReactionFormula(metabolite_v7_6$`Metabolite description`, metabolite_v7_6$`Metabolite name`, metabolite$`Metabolite name`)
metabolite$met_consistent <- metabolite$`Metabolite description` == metabolite$Description_v_7_6




# update the metabolite description
# mainly using the short name of the compartment
# also remove the space
metabolite0 <- metabolite %>% separate(`Metabolite name`, into = c('Abbreviation','compartment'), sep = '\\[')
metabolite0$compartment <- paste('[',metabolite0$compartment, sep = "")
metabolite0$Description <- str_replace_all(metabolite0$`Metabolite description`, " \\[.*?\\]", "")
metabolite0$Description <- paste(metabolite0$Description, metabolite0$compartment, sep = "")
metabolite$Description <- metabolite0$Description

#using the detailed name to take place of the short name
rxn_split <- splitAndCombine0(rxn$Reaction,rxn$Abbreviation,sep=" ")
rxn_split$v3 <- getSingleReactionFormula(metabolite$Description,metabolite$Abbreviation,rxn_split$v1)

for (i in 1:length(rxn_split$v2)){
  if(rxn_split$v3[i]=="NA"){
    rxn_split$v3[i] <- rxn_split$v1[i]
  } else{
    rxn_split$v3[i] <- rxn_split$v3[i]
  }
  
}

rxn$Description_new <- getMultipleReactionFormula(rxn_split$v3,rxn_split$v2,rxn$Abbreviation)
rxn$Description_new <- str_replace_all(rxn$Description_new,";"," ")


#merge with old and new version of subsystem
yeastGEM_old_version <- read_excel("data/yeastGEM_latest version.xls")
yeastGEM_new_version <- read_excel("/Users/luho/PycharmProjects/model/model_correction/result/subsystem for yeast8 map_V3.xlsx")
rxn$Subsystem_early_version <- getSingleReactionFormula(yeastGEM_old_version$Subsystem_new,yeastGEM_old_version$Abbreviation,rxn$Abbreviation)
rxn$Subsystem_new <-  getSingleReactionFormula(yeastGEM_new_version$subsystem_map_v3,yeastGEM_new_version$rxnID,rxn$Abbreviation)
rxn$reaction_early_version <-  getSingleReactionFormula(yeastGEM_old_version$Description_new,yeastGEM_old_version$Abbreviation,rxn$Abbreviation)

#quality check
rxn$subsystem_consistence <- rxn$Subsystem_early_version==rxn$Subsystem_new
rxn$reaction_consistence <- rxn$reaction_early_version==rxn$Description_new

#rxn formula need check
rxn_formula_need_check <- filter(rxn, reaction_consistence==FALSE)

#save the model
write.table(rxn, "result/yeastGEM_october.txt", row.names = FALSE, sep = "\t")
write.table(rxn_formula_need_check, "result/rxn_formula_need_check.txt", row.names = FALSE, sep = "\t")


#comparsion the subsytem in the new and old version
num_subsystem_new <- length(unique(rxn$Subsystem_new[str_detect(rxn$Subsystem_new,"transport")==FALSE]))

#comparsion the subsytem in the new and old version
num_subsystem_old <- length(unique(rxn$Subsystem_early_version[str_detect(rxn$Subsystem_early_version,"transport")==FALSE]))

setdiff(unique(rxn$Subsystem_new[str_detect(rxn$Subsystem_new,"transport")==FALSE]), unique(rxn$Subsystem_early_version[str_detect(rxn$Subsystem_early_version,"transport")==FALSE]))
