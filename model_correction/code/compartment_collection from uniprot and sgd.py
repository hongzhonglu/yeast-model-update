# Import packages
import numpy as np
import pandas as pd
import os    ##for directory
import sys
import pprint
# set the directory
os.chdir('/Users/luho/PycharmProjects/model/model_correction/code')
sys.path.append(r"/Users/luho/PycharmProjects/model/cobrapy/code")
pprint.pprint(sys.path)
# import self function
from mainFunction import *

#input the data
compartment = pd.read_excel('../data/Compartment.xlsx')
compartment_SGD = pd.read_excel('../data/Protein location SGD.xlsx')


# gene compartment annotation for SGD database
"""input the compartment annotation from SGD"""
compartment_SGD['go_type'] = compartment_SGD['go_type'].fillna('NA')
compartment_SGD0 = compartment_SGD[compartment_SGD['go_type'].str.contains("component")]
"""refine the compartment annotation from SGD"""
unique_cp_SGD = compartment_SGD0['go_term'].unique()
unique_cp_SGD = unique_cp_SGD.tolist()
unique_cp_SGD = pd.DataFrame(unique_cp_SGD)
unique_cp_SGD.columns = ['compartment']
unique_cp_SGD.columns = unique_cp_SGD.columns.str.strip()
"""standardize location information from SGD"""
ss0 = len(unique_cp_SGD)
similarTarget0 = getSimilarTarget(compartment[['description']],unique_cp_SGD[['compartment']],ss=ss0)
unique_cp_SGD['similar_target'] = similarTarget0
unique_cp_SGD.head()
"""save the standard location information from uniprot for manual check"""
writer = pd.ExcelWriter('../data/unique_compartment_SGD.xlsx')
unique_cp_SGD.to_excel(writer,'Sheet1')
writer.save()
"""replace the compartment information using standard name from model"""
xls = pd.ExcelFile('../data/Compartment.xlsx')
sgdChange = pd.read_excel(xls, 'Sheet2')
sgdChange['SGD'] = sgdChange['SGD'].str.strip()
sgdChange['model_name'] = sgdChange['model_name'].str.strip()

#for the 'membrane', we change it from "not sure" type into "membrane", which will be useful to assign the gene for these transport reactions
for i in range(len(sgdChange['SGD'])):
    print(i)
    if sgdChange.loc[i,'SGD'] == 'membrane':
        sgdChange.loc[i,'model_name'] ='membrane'
    else:
        sgdChange.loc[i, 'model_name'] =  sgdChange.loc[i, 'model_name']


compartment_SGD0['compartment'] = singleMapping(sgdChange['model_name'],sgdChange['SGD'],compartment_SGD0['go_term'])
gene_compartmentSGD = pd.DataFrame({
    'gene': compartment_SGD0['systematic_name'].str.strip().unique()
})

gene_compartmentSGD['compartment'] = multiMapping(compartment_SGD0['compartment'],compartment_SGD0['systematic_name'],gene_compartmentSGD['gene'])








# gene compartment annotation for uniprot database
"""input the gene annotation from uniprot database"""
xls = pd.ExcelFile('../data/uniprot_location.xlsx')
uniprot_location = pd.read_excel(xls, 'Sheet1')
uniprot_location['Subcellular location [CC]'] = uniprot_location['Subcellular location [CC]'].str.replace('SUBCELLULAR LOCATION: ','')
uniprot_location['Subcellular location [CC]'] = uniprot_location['Subcellular location [CC]'].str.replace(';',',')
uniprot_location['Subcellular location [CC]'] = uniprot_location['Subcellular location [CC]'].str.replace('.',',')
uniprot_location['Subcellular location [CC]'] = uniprot_location['Subcellular location [CC]'].str.replace(r'\{.+?\}','')

uniprot_location0 = uniprot_location.iloc[:,1].str.split('Note=', expand=True)
uniprot_location['location'] = uniprot_location0.iloc[:,0]
uniprot_location1 = splitAndCombine(uniprot_location['location'],uniprot_location['gene'],sep0=",")

uniprot_location2 = uniprot_location1[uniprot_location1.V2 != '']
uniprot_location2 = uniprot_location2[uniprot_location2.V2 != ' ']
compartment_uniprot = uniprot_location2[uniprot_location2.V2 != 'NA']
compartment_uniprot.columns = ['gene','compartment']
compartment_uniprot.gene = compartment_uniprot.gene.str.strip()
compartment_uniprot.compartment = compartment_uniprot.compartment.str.strip()

unique_compartment_uniprot = pd.Series(uniprot_location2.V2.unique())
unique_compartment_uniprot = unique_compartment_uniprot.str.strip()
unique_compartment_uniprot = unique_compartment_uniprot.unique()
unique_compartment_uniprot = pd.DataFrame(unique_compartment_uniprot)
unique_compartment_uniprot.columns = ['compartment']

"""standardize location information from uniprot database"""
ss0 = len(unique_compartment_uniprot)
similarTarget0 = getSimilarTarget(compartment[['description']],unique_compartment_uniprot[['compartment']],ss=ss0)
unique_compartment_uniprot['similar_target'] = similarTarget0
unique_compartment_uniprot.head()

"""save the standard location information from uniprot for manual check"""
writer = pd.ExcelWriter('../data/unique_compartment_uniprot.xlsx')
unique_compartment_uniprot.to_excel(writer,'Sheet1')
writer.save()

"""replace the compartment information using standard name from model for uniprot gene location annotation"""
xls = pd.ExcelFile('../data/Compartment.xlsx')
uniChange = pd.read_excel(xls, 'Sheet3')
uniChange['uniprot'] = uniChange['uniprot'].str.strip()
uniChange['model_name'] = uniChange['model_name'].str.strip()

#for the 'membrane', we change it from "not sure" type into "membrane", which will be useful to assign the gene for these transport reactions
for i in range(len(uniChange['uniprot'])):
    print(i)
    if uniChange.loc[i,'uniprot'] == 'Membrane':
        uniChange.loc[i,'model_name'] ='membrane'
    else:
        uniChange.loc[i, 'model_name'] =  uniChange.loc[i, 'model_name']




compartment_uniprot['model_name'] = singleMapping(uniChange['model_name'],uniChange['uniprot'],compartment_uniprot['compartment'])

gene_compartmentUNI = pd.DataFrame({
    'gene': compartment_uniprot['gene'].str.strip().unique()
})

gene_compartmentUNI['compartment'] = multiMapping(compartment_uniprot['model_name'],compartment_uniprot['gene'],gene_compartmentUNI['gene'])



saveExcel(gene_compartmentSGD,'../result/gene_compartmentSGD_updated_october.xlsx')
saveExcel(gene_compartmentUNI,'../result/gene_compartmentUNI_updated_october.xlsx')