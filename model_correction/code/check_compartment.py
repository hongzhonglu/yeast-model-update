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
yeastGEM = pd.read_excel('../data/yeastGEM_latest version.xls')
compartment = pd.read_excel('../data/Compartment.xlsx')
compartment_SGD = pd.read_excel('../data/Protein location SGD.xlsx')



#establish reaction-single gene-compartment pair
list(yeastGEM) #obtain the column names of dataframe
yeastGEM.columns #obtain the column names of dataframe
yeastGEM0 = yeastGEM.copy()
yeastGEM0 = yeastGEM0.iloc[:,[1,4,5]] #obtain the column 1, 4, 5
yeastGEM0['GPR'] = yeastGEM0['GPR'].str.replace('or',';') # string replace
yeastGEM0['GPR'] = yeastGEM0['GPR'].str.replace('and',';')
yeastGEM0['GPR'] = yeastGEM0['GPR'].str.replace('(','')
yeastGEM0['GPR'] = yeastGEM0['GPR'].str.replace(')','')
yeastGEM_s = splitAndCombine(yeastGEM0['GPR'],yeastGEM0['Abbreviation'],sep0=";")# split the genes in GPR relation, so a reaction is related to only one gene
yeastGEM_s.columns = ['Abbreviation','gene'] #name the column
yeastGEM_s.gene = yeastGEM_s.gene.str.strip() #remove the white space
yeastGEM_s['rxn'] = singleMapping(yeastGEM0.iloc[:,1].tolist(),yeastGEM0.iloc[:,0].tolist(),yeastGEM_s.iloc[:,0].tolist(), dataframe=False)#obtain the detailed rxn for the pair of gene-reactionID

"""extract the compartment information from each reaction"""
compartment['name'] = "[" + compartment['name'] + "]"
cp10 = compartment['name'].tolist()
cp20 = compartment['description'].tolist()
rxn = yeastGEM_s['rxn'].tolist()


yeastGEM_s['compartment'] = [None]*len(yeastGEM_s)
for j in range(len(rxn)):
    yeastGEM_s['compartment'][j] = getCompartment(rxn[j])




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

compartment_uniprot['model_name'] = singleMapping(uniChange['model_name'],uniChange['uniprot'],compartment_uniprot['compartment'])

gene_compartmentUNI = pd.DataFrame({
    'gene': compartment_uniprot['gene'].str.strip().unique()
})

gene_compartmentUNI['compartment'] = multiMapping(compartment_uniprot['model_name'],compartment_uniprot['gene'],gene_compartmentUNI['gene'])



#once the gene compartment information from SGD and uniprot have been given, then they will be compared with the yeastGEM
"""merge the SGD and uniprot compartment information with yeastGEM"""
yeastGEM_s['compartment_sgd'] = singleMapping(gene_compartmentSGD['compartment'],gene_compartmentSGD['gene'], yeastGEM_s['gene'])
yeastGEM_s['compartment_uni'] = singleMapping(gene_compartmentUNI['compartment'],gene_compartmentUNI['gene'],yeastGEM_s['gene'])
yeastGEM_s0 = yeastGEM_s[(yeastGEM_s['gene'] != 'NA')]
"""combine sgd and uniprot"""
yeastGEM_s0.compartment_sgd = yeastGEM_s0.compartment_sgd.astype(str)
yeastGEM_s0.compartment_uni = yeastGEM_s0.compartment_uni.astype(str)
yeastGEM_s0['combine'] = yeastGEM_s0.compartment_sgd.str.cat(yeastGEM_s0.compartment_uni, sep=';')
yeastGEM_s0['combine'] = yeastGEM_s0['combine'].astype(str)
yeastGEM_s0['combine_unique'] = [None]*len(yeastGEM_s0['combine'])
"""remove the duplicated compartment information"""
for index, row in yeastGEM_s0.iterrows():
    print(index)
    if ';' in yeastGEM_s0['combine'][index]:
        ss = yeastGEM_s0['combine'][index].split(';')
        ss = pd.Series(ss)
        ss = ss.unique()
        ss = ss.tolist()
        ss = ';'.join(ss)
        yeastGEM_s0['combine_unique'][index] = ss
    else:
        yeastGEM_s0['combine_unique'][index] = yeastGEM_s0['combine'][index]

"""find the common compartment in the model and database annotation"""
sss = [None]*len(yeastGEM_s0['gene'])
for i in range(len(yeastGEM_s0['gene'])):
     sss[i] = getCommonCompartment(yeastGEM_s0['compartment'].tolist()[i],yeastGEM_s0['combine_unique'].tolist()[i])

yeastGEM_s0['common_compartment'] = sss



# save the related results
writer = pd.ExcelWriter('../result/comparment_SGD.xlsx')
compartment_SGD0.to_excel(writer,'Sheet1')
writer.save()

writer = pd.ExcelWriter('../result/compartment_uniprot.xlsx')
compartment_uniprot.to_excel(writer,'Sheet1')
writer.save()

writer = pd.ExcelWriter('../result/gene_compartment_sgd.xlsx')
gene_compartmentSGD.to_excel(writer,'Sheet1')
writer.save()

writer = pd.ExcelWriter('../result/gene_compartment_uniprot.xlsx')
gene_compartmentUNI.to_excel(writer,'Sheet1')
writer.save()

writer = pd.ExcelWriter('../result/compartment check result.xlsx')
yeastGEM_s0.to_excel(writer,'Sheet1')
writer.save()






# other tasks
#####code to check the new added gene
"""input the data"""
newGPR = pd.read_excel('../data/compartment_newGPR.xlsx')
newGPR['compartment_uni'] = singleMapping(gene_compartmentUNI['compartment'],gene_compartmentUNI['gene'],newGPR['gene'])
newGPR['compartment_sgd'] = singleMapping(gene_compartmentSGD['compartment'],gene_compartmentSGD['gene'], newGPR['gene'])

#newGPR['compartment_uni'] = newGPR['compartment_uni'].replace(None, 'NONE')

sss0 = [None]*len(newGPR['gene'])
for i in range(len(newGPR['gene'])):
    sss0[i] = getCommonCompartment(newGPR['compartment_uni'].tolist()[i],newGPR['compartment_sgd'].tolist()[i])


newGPR['common_compartment'] = sss0
newGPR['common_compartment'] = newGPR['common_compartment'].str.replace('not sure;','')
newGPR['common_compartment'] = newGPR['common_compartment'].str.replace(';not sure','')
newGPR['common_compartment'] = newGPR['common_compartment'].str.replace('not sure','')

newGPR['compartment_uni'] = newGPR['compartment_uni'].str.replace('not sure;','')
newGPR['compartment_uni'] = newGPR['compartment_uni'].str.replace(';not sure','')
newGPR['compartment_uni'] = newGPR['compartment_uni'].str.replace('not sure','')

newGPR['compartment_sgd'] = newGPR['compartment_sgd'].str.replace('not sure;','')
newGPR['compartment_sgd'] = newGPR['compartment_sgd'].str.replace(';not sure','')
newGPR['compartment_sgd'] = newGPR['compartment_sgd'].str.replace('not sure','')
newGPR['compartment_sgd'] = newGPR['compartment_sgd'].str.replace(';None','')
newGPR['compartment_sgd'] = newGPR['compartment_sgd'].str.replace('None;','')



"""output result"""
writer = pd.ExcelWriter('../result/compartment check result for new GPR in the model update.xlsx')
newGPR.to_excel(writer,'Sheet1')
writer.save()


"""using these new GPR information to check the reaction/genes added into model"""
geneOnly = pd.read_csv('../data/DBaddonlyGenes.txt', sep="\t")
rxnOnly = pd.read_csv('../data/DBaddGenesAndRxns.txt', sep="\t")

mappingNewRXN = getRXNgeneMapping(rxnOnly['RxnID'],rxnOnly['GPR'])
mappingNewRXN['compartment_uni'] = singleMapping(newGPR['compartment_uni'],newGPR['gene'],mappingNewRXN['gene'])
mappingNewRXN['compartment_sgd'] = singleMapping(newGPR['compartment_sgd'],newGPR['gene'],mappingNewRXN['gene'])
mappingNewRXN['common_compartment'] = singleMapping(newGPR['common_compartment'],newGPR['gene'],mappingNewRXN['gene'])

mapNewGene1 = getRXNgeneMapping(geneOnly['ID'],geneOnly['GPR'])
mapNewGene2 = getRXNgeneMapping(geneOnly['ID'],geneOnly['final_GPR'])

mapNewGene1['join1'] = mapNewGene1['rxnID'] +'@' + mapNewGene1['gene']
mapNewGene2['join1'] = mapNewGene2['rxnID'] +'@' + mapNewGene2['gene']

x1 = np.array(mapNewGene1['join1'])
x2 = np.array(mapNewGene2['join1'])
x3 = np.setdiff1d(x2, x1)

newGene = pd.DataFrame({'gene':pd.Series(x3)})
newGene0 = newGene.iloc[:,0].str.split('@', expand=True)
newGene0.columns = ['ID','newGene']

#get the rxn compartment from model and get gene compartment from annotations
newGene0['rxn_compartment'] = singleMapping(yeastGEM_s0['compartment'],yeastGEM_s0['Abbreviation'],newGene0['ID'])

newGene0['compartment_uni'] = singleMapping(newGPR['compartment_uni'],newGPR['gene'],newGene0['newGene'])
newGene0['compartment_sgd'] = singleMapping(newGPR['compartment_sgd'],newGPR['gene'],newGene0['newGene'])
newGene0['common_compartment'] = singleMapping(newGPR['common_compartment'],newGPR['gene'],newGene0['newGene'])


"""output result"""
writer = pd.ExcelWriter('../result/compartment check result for new gene.xlsx')
newGene0.to_excel(writer,'Sheet1')
writer.save()

writer = pd.ExcelWriter('../result/compartment check result for new rxn.xlsx')
mappingNewRXN.to_excel(writer,'Sheet1')
writer.save()


"""update the gene compartment information for all the new genes"""
allNewGene = pd.read_csv('../data/DBnewGeneAnnotation.tsv', sep="\t")
allNewGene['compartment_uni'] = singleMapping(newGPR['compartment_uni'],newGPR['gene'],allNewGene['gene'])
allNewGene['compartment_sgd'] = singleMapping(newGPR['compartment_sgd'],newGPR['gene'],allNewGene['gene'])
allNewGene['common_compartment'] = singleMapping(newGPR['common_compartment'],newGPR['gene'],allNewGene['gene'])

writer = pd.ExcelWriter('../result/DBnewGeneAnnotation.xlsx')
allNewGene.to_excel(writer,'Sheet1')
writer.save()




####small tasks
#input the data
yeastGEM = pd.read_excel('/Users/luho/Google Drive/R application and code/yeast map and visualization of multi-omics/data/yeastGEM_october.xls')
#establish reaction-single gene-compartment pair
yeastGEM0 = yeastGEM.copy()
"""extract the compartment information from each reaction"""
rxn = yeastGEM0['Description_new'].tolist()
yeastGEM0['compartment'] = [None]*len(yeastGEM0)
for j in range(len(rxn)):
    yeastGEM0['compartment'][j] = getCompartment(rxn[j])
saveExcel(yeastGEM0,'../result/yeastGEM_october with compartment for each reaction.xlsx')