# -*- coding: utf-8 -*-
# -*- python 3 -*-
# -*- hongzhong Lu -*-


# Import packages
import numpy as np
import pandas as pd
import os    ##for directory
import sys
import pprint
from cobra.io import read_sbml_model, write_sbml_model,save_json_model
import re
import cobrapy
from cobra import Model, Reaction, Metabolite

os.chdir('/Users/luho/PycharmProjects/model/model_correction/code')
sys.path.append(r"/Users/luho/PycharmProjects/model/cobrapy/code")
pprint.pprint(sys.path)

# import self function
from mainFunction import *

"""step1 prepare the format of model in dataframe"""
model2 = read_sbml_model('/Users/luho/Documents/GitHub/YeastMetabolicNetwork-GEM/ModelFiles/xml/yeastGEM.xml')

# Correct metabolite and gene ids in wrong formats
model2= correctSomeWrongFormat(model2)

#produce the dataframe for the metabolites
met_dataframe = produceMetaboliteList(model2)
saveExcel(met_dataframe, '../result/met_yeastGEM.xlsx')

#produce the dataframe for the rxn
gem_dataframe = produceRxnList(model2)


'''step2 obtain the subsystem from kegg based on gene mapping'''
#get subsystem based on the gene
#input the kegg subsystem annotation
#kegg_subsystem = pd.read_table('/Users/luho/Documents/GitHub/YeastMetabolicNetwork-GEM/ComplementaryData/databases/kegg.tsv')
gene_pathwwayID = pd.read_excel('../data/sce_kegg.xlsx', sheet_name='gene_pathway')
pathway = pd.read_excel('../data/sce_kegg.xlsx', sheet_name='pathwayList')
pathway['pathwayName'] = pathway['pathwayID'] + '  ' + pathway['pathwayName']
gene_pathwwayID['pathway'] = singleMapping(pathway['pathwayName'],pathway['pathwayID'],gene_pathwwayID['pathwayID'])


rxn_gene = getRXNgeneMapping(gem_dataframe['rxnID'], gem_dataframe['GPR'])
rxn_gene['subsystem_automatic'] = multiMapping(gene_pathwwayID['pathway'], gene_pathwwayID['geneID'],rxn_gene['gene'], sep=" // ")
gem_dataframe['subsytem_automatic'] = multiMapping(rxn_gene['subsystem_automatic'],rxn_gene['rxnID'],gem_dataframe['rxnID'], sep=" // ")
gem_dataframe['subsytem_automatic'] = gem_dataframe['subsytem_automatic'].str.replace("nan //","")
gem_dataframe['subsytem_automatic'] = gem_dataframe['subsytem_automatic'].str.replace("// nan","")
gem_dataframe['subsytem_automatic'] = gem_dataframe['subsytem_automatic'].str.replace("None //","")
gem_dataframe['subsytem_automatic'] = gem_dataframe['subsytem_automatic'].str.replace("// None","")

gem_dataframe['subsytem_automatic'] = RemoveDuplicated(gem_dataframe['subsytem_automatic'].tolist())


'''step3 Add exchange and transport'''
#exchange
gem_dataframe['subsytem_automatic'] = exchange(gem_dataframe['formula'].tolist(),gem_dataframe['subsytem_automatic'].tolist())

#SLIME
gem_dataframe['subsytem_automatic'] = SLIME(gem_dataframe['name'].tolist(),gem_dataframe['subsytem_automatic'].tolist())

#transport
#note: some reactions were not found , for example,
#ATP [cytoplasm] + H2O [cytoplasm] + cadmium(2+) [cytoplasm] --> ADP [cytoplasm] + H+ [cytoplasm] + phosphate [cytoplasm] + cadmium(2+) [extracellular]	P-type cation-transporting ATPase (EC 3.6.3.3) (Cadmium resistance protein 2) (Cadmium-translocating P-type ATPase) (Cd(2+)-exporting ATPase)	r_4172
#the function 'transport' need redesigned
gem_dataframe['subsytem_automatic'] = transport(gem_dataframe['formula'].tolist(),gem_dataframe['subsytem_automatic'].tolist())




'''step4'''
#automatic find the subsystem for the reaction with genes but without subsystem
#add source
gem_dataframe1 = gem_dataframe[(gem_dataframe["subsytem_automatic"]=="None") | (gem_dataframe["subsytem_automatic"]=="nan")]
rxn_gene2 = getRXNgeneMapping(gem_dataframe1['rxnID'], gem_dataframe1['GPR'])
#find subsystem based on biocyc, uniprot and reactome
gene_subsystem = pd.read_excel('../data/geneSubsystem.xlsx')
rxn_gene2['subystem_biocyc'] = multiMapping(gene_subsystem['subsystem_biocyc'],gene_subsystem['gene'],rxn_gene2['gene'])
rxn_gene2['subystem_reactome'] = multiMapping(gene_subsystem['subsystem_reactome'],gene_subsystem['gene'],rxn_gene2['gene'])
rxn_gene2['subystem_notKEGG'] =[None]*len(rxn_gene2['rxnID'])
for i in range(0,len(rxn_gene2['gene'])):
    if rxn_gene2['subystem_biocyc'][i] =='None':
        rxn_gene2['subystem_notKEGG'][i] = nz(rxn_gene2['subystem_reactome'][i]) + '_reactome'
    elif rxn_gene2['subystem_biocyc'][i] =='nan':
        rxn_gene2['subystem_notKEGG'][i] = nz(rxn_gene2['subystem_reactome'][i]) + '_reactome'
    else:
        rxn_gene2['subystem_notKEGG'][i] = nz(rxn_gene2['subystem_biocyc'][i]) + '_biocyc'

rxn_gene2['subystem_notKEGG'] = rxn_gene2['subystem_notKEGG'].str.replace('nan_biocyc', '')
rxn_gene2['subystem_notKEGG'] = rxn_gene2['subystem_notKEGG'].str.replace('none_biocyc','')
rxn_gene2['subystem_notKEGG'] = rxn_gene2['subystem_notKEGG'].str.replace('nan_reactome','')
gem_dataframe1['subsytem_automatic'] = multiMapping(rxn_gene2['subystem_notKEGG'],rxn_gene2['rxnID'],gem_dataframe1['rxnID'])
#merge subsystem information from biocyc and reactome with kegg
s2= updateOneColumn(df1=gem_dataframe, df2=gem_dataframe1, key0='rxnID', value0='subsytem_automatic')
gem_dataframe['subsystem_xref'] =s2


'''step5'''
#correct the subsystem information based on kegg subsystem definition
rxn_mnxid = pd.read_excel('../data/rxnlistMNX.xlsx')
mnxid_xref = pd.read_table('../data/reac_xref.tsv')
mnxid_xref_kegg = mnxid_xref[mnxid_xref['XREF'].str.contains('kegg')]
mnxid_xref_rhea = mnxid_xref[mnxid_xref['XREF'].str.contains('rhea')]
rxn_mnxid['keggID2'] = singleMapping(mnxid_xref_kegg['XREF'],mnxid_xref_kegg['MNX_ID'],rxn_mnxid['MNXID'])
rxn_mnxid['rheaID'] = singleMapping(mnxid_xref_rhea['XREF'],mnxid_xref_rhea['MNX_ID'],rxn_mnxid['MNXID'])

#find the keggid based on the rheaID
rhea_kegg = pd.read_table('../data/rhea2kegg_reaction.tsv')
rhea_kegg['MASTER_ID'] = 'rhea:'+ pd.Series(rhea_kegg['MASTER_ID']).astype(str)
rxn_mnxid['keggID3'] = singleMapping(rhea_kegg['ID'],rhea_kegg['MASTER_ID'],rxn_mnxid['rheaID'])
rxn_mnxid['keggID'] = rxn_mnxid['keggID2']
rxn_mnxid['keggID1'] = 'kegg:' + rxn_mnxid['keggID1'].astype(str)
rxn_mnxid['keggID1'] = rxn_mnxid['keggID1'].replace('kegg:nan', np.nan)
for i in range(0, len(rxn_mnxid['keggID'])):
    if rxn_mnxid['keggID2'][i] is None:
        rxn_mnxid.loc[i, 'keggID']= rxn_mnxid.loc[i,'keggID1']
    else:
        rxn_mnxid.loc[i, 'keggID'] = rxn_mnxid.loc[i, 'keggID']

rxn_mnxid['keggID'] = rxn_mnxid['keggID'].str.replace('#2','')
rxn_mnxid0 = rxn_mnxid[['rxnID','MNXID','keggID','rheaID','kegg_subsystem']]

reactionID_pathway = pd.read_table('../data/reactionID_pathway.txt', header=None)
reactionID_pathway.columns=['reactionID','pathway']
reactionID_pathway0 = reactionID_pathway[reactionID_pathway['pathway'].str.contains('path:map')]
reactionID_pathway0['reactionID']=reactionID_pathway0['reactionID'].str.replace('rn:','kegg:')
reactionID_pathway0['pathway']=reactionID_pathway0['pathway'].str.replace('path:map','sce')
reactionID_pathway0['pathway0']=singleMapping(pathway['pathwayName'],pathway['pathwayID'],reactionID_pathway0['pathway'])
rxn_mnxid0['subsystem'] = multiMapping(reactionID_pathway0['pathway0'],reactionID_pathway0['reactionID'],rxn_mnxid0['keggID'], sep=" // ")

'''step6'''
#compare the subsystem obtained with gene and reactionID from kegg database
gem_dataframe['subsystem_rxnID'] = singleMapping(rxn_mnxid0['subsystem'],rxn_mnxid0['rxnID'],gem_dataframe['rxnID'])
sss = [None]*len(gem_dataframe['subsystem_rxnID'])
for i in range(len(gem_dataframe['subsystem_rxnID'])):
     sss[i] = getCommonCompartment(gem_dataframe['subsystem_rxnID'].tolist()[i],gem_dataframe['subsytem_automatic'].tolist()[i], sep0=" // ")
gem_dataframe['subsytem_automatic_common'] = sss


#save the results
saveExcel(gem_dataframe, '../result/yeastGEM_with subsystem.xlsx')










##small task for bigg
'''bigg = pd.read_excel('../data/gene_rxn_eggnog for manual check.xlsx')
s1= bigg['formula']
subsystem=list()
    for i, x0 in enumerate(s1):
        x0 = re.sub('_[^0-9]', '', x0)
        print(i, x0)
        x2 = x0
        if "<=>" in x2:
            x3 = x2.split("<=>")
        elif "<->" in x2: #bigg database format
            x3 = x2.split("<->")
        else:
            x3 = x2.split("-->")
        x3 = [x.strip() for x in x3]
        print(x3)

        if '+' in x3[0]:
            x30=x3[0].split('+')
        else:
            x30=x3[0]
        print(x30)

        x30=[x.strip() for x in x30]
        if '+' in x3[1]:
            x31 = x3[1].split('+')
        else:
            x31=x3[1]
        x31 = [x.strip() for x in x31]

        if (set(x30)-set(['atp','adp'])) == (set(x31)-set(['atp','adp'])):
            subsystem.append('Transport')
        else:
            subsystem.append('none')

bigg['subsytem'] = subsystem
#save the results
saveExcel(bigg, '../result/bigg subsystem.xlsx')'''



