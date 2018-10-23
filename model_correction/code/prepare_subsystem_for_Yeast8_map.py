# -*- coding: utf-8 -*-
# -*- python 3 -*-
# -*- hongzhong Lu -*-
import os
os.chdir('/Users/luho/PycharmProjects/model/model_correction/code')
exec(open("./find_subsystem_Yeast8_using_code.py").read())
#it can be found that the reaction in different software is different in some formats
#thus the reaction list will be based on R function to keep the consistent
subsystem_map = pd.read_excel('../result/subsystem_manual check results.xlsx')
gem_dataframe['subsystem_map'] = singleMapping(subsystem_map['Subsystem_for yeast map'],subsystem_map['Abbreviation'],gem_dataframe['rxnID'])
gem_dataframe['removed_duplicate_subsystem'] = singleMapping(subsystem_map['removed_duplicate_subsystem'],subsystem_map['Abbreviation'],gem_dataframe['rxnID'])
gem_dataframe['evidence'] = singleMapping(subsystem_map['evidence'],subsystem_map['Abbreviation'],gem_dataframe['rxnID'])

#add the subsystem obtained based on geneID
for i in range(0,len(gem_dataframe['subsystem_map'])):
    if gem_dataframe['subsystem_map'][i] is None:
        gem_dataframe['subsystem_map'][i] = gem_dataframe['subsystem_xref'][i]
    else:
        gem_dataframe['subsystem_map'][i] =  gem_dataframe['subsystem_map'][i]



#add the subsystem obtained based on the keggID
for i in range(0,len(gem_dataframe['subsystem_map'])):
    if gem_dataframe['subsystem_map'][i] is '':
        gem_dataframe['subsystem_map'][i] = gem_dataframe['subsystem_rxnID'][i]
    else:
        gem_dataframe['subsystem_map'][i] =  gem_dataframe['subsystem_map'][i]


#add the information from manual check results for these reactions connected with new genes
subsystem_manual = pd.read_excel('../data/subsytem_for new genes added into Yeast8.xlsx')
subsystem_manual['inf'] = subsystem_manual.loc[:,'subpathway'] + ' @@ ' + subsystem_manual.loc[:,'note']
rxn_gene['subsystem_manual'] = multiMapping(subsystem_manual['inf'],subsystem_manual['gene'],rxn_gene['gene'],sep=" // ")
gem_dataframe['subsytem_manual_newGene'] = multiMapping(rxn_gene['subsystem_manual'] ,rxn_gene['rxnID'] ,gem_dataframe['rxnID'],sep=" // ")
gem_dataframe['subsytem_manual_newGene'] = RemoveDuplicated(gem_dataframe['subsytem_manual_newGene'].tolist())


#add the information from reaction notes for these reactions from biolog experiments
evidences_biolog = pd.read_excel('../data/classification for new reactions from biolog_result.xlsx')
evidences_biolog['inf'] = evidences_biolog.loc[:,'source'] + ' @@ ' + evidences_biolog.loc[:,'external_ID']
gem_dataframe['note'] = multiMapping(evidences_biolog['inf'], evidences_biolog['rxnID'] ,gem_dataframe['rxnID'],sep=" // ")
saveExcel(gem_dataframe,"../result/subsystem for yeast8 map.xlsx")


#refine the subsystem for the yeast map based on the reaction number and manual check results
subsystem_map_v2 = pd.read_excel('../result/subsystem for yeast8 map_V2.xlsx')












































