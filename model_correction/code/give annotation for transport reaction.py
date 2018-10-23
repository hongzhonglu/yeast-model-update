# -*- coding: utf-8 -*-
# -*- python 3 -*-
# -*- hongzhong Lu -*-
import os
import sys
os.chdir('/Users/luho/PycharmProjects/model/model_correction/code')
sys.path.append(r"/Users/luho/PycharmProjects/model/cobrapy/code")

# import self function
from mainFunction import *

#step one
#add the genes for the dipeptide uptake
transport_rxn = pd.read_excel('../data/transport annotation for yeast8.xlsx', 'Sheet1')
#dipeptide name
dipetides = pd.read_excel('../data/transport annotation for yeast8.xlsx', 'Sheet2')
#YKR093W   Uptake of small peptides.
#if the transport reactions contains the above metabolite, then this gene will be added for the related reactions
rxn_met = splitAndCombine(transport_rxn['formula'], transport_rxn['rxnID'], sep0=" ")
rxn_met['V2'] = rxn_met['V2'].str.strip()
rxn_met0 = rxn_met[rxn_met['V2'].isin(dipetides['dipeptide'])]
transport_rxn['has_dipeptide'] = transport_rxn['rxnID'].isin(rxn_met0['V1'])
saveExcel(transport_rxn,"../result/transport annotation for yeast8.xlsx")








































