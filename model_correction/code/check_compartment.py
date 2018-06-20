# Import packages
import numpy as np
import pandas as pd
import os    ##for directory
# set the directory
os.chdir('/Users/luho/PycharmProjects/python learning/venv/project1_modelling/model_correction/code')
os.getcwd()

"""input the data"""
yeastGEM = pd.read_excel('../data/yeastGEM_latest version.xls')
compartment = pd.read_excel('../data/Compartment.xlsx')
compartment_SGD = pd.read_excel('../data/Protein location SGD.xlsx')


"""establish reaction-single gene-compartment pair"""
list(yeastGEM)
yeastGEM0 = yeastGEM.copy()
yeastGEM0 = yeastGEM0.iloc[:,[1,4,5]]

yeastGEM0['GPR'] = yeastGEM0['GPR'].str.replace('or',';')
yeastGEM0['GPR'] = yeastGEM0['GPR'].str.replace('and',';')
yeastGEM0['GPR'] = yeastGEM0['GPR'].str.replace('(','')
yeastGEM0['GPR'] = yeastGEM0['GPR'].str.replace(')','')


"""split and combine"""
def splitAndCombine(gene, rxn, sep0, moveDuplicate=False):
    
## one rxn has several genes, this function was used to splite the genes
## used for the dataframe data
    gene = gene.fillna('NA') # fill the NaN with 'NA'
    gene0 = gene.tolist()
    rxn0 = rxn.tolist()
    s1 = list()
    s2 = list()
    for i in range(len(gene0)):
       s1 = s1 + [rxn0[i]]*len(gene0[i].split(sep0))
       s2 = s2 + gene0[i].split(sep0)
    df0 = pd.DataFrame({'V1': s1,
                   'V2':s2}
                   )
    if moveDuplicate==True:
       df00 = df0.drop_duplicates()
    else:
       df00 = df0

    return df00

yeastGEM_s = splitAndCombine(yeastGEM0['GPR'],yeastGEM0['Abbreviation'],sep0=";")
yeastGEM_s.columns = ['Abbreviation','gene']


"""mapping"""
def singleMapping (description, item1, item2):
    #description = w
    #item1 = v
    #item2 = testData
    # used for the list data
    index = [None]*len(item2)
    result = [None]*len(item2)
    tt = [None]*len(item2)
    for i in range(len(item2)):
        if item2[i] in item1:
            index[i] = item1.index(item2[i])
            result[i] = description[index[i]]
        else:
            index[i] = None
            result[i] = None
    return result


def multiMapping (description, item1, item2):
     #description = w
     #item1 = v
     #item2 = testData
     #used for the list data
    index = [None]*len(item2)
    result = [None]*len(item2)
    for i in range(len(item2)):
        if item2[i] in item1:
            index0 = [description[index] for index in range(len(item1)) if item1[index] == item2[i]]
            index1 = pd.unique(index0).tolist()
            result[i] = ';'.join(str(e) for e in index1) #string cat
        else:
            result[i] = None
    return result


yeastGEM_s['rxn'] = singleMapping(yeastGEM0.iloc[:,1].tolist(),yeastGEM0.iloc[:,0].tolist(),yeastGEM_s.iloc[:,0].tolist())



"""extract the compartment information from each reaction"""
compartment['name'] = "[" + compartment['name'] + "]"
cp10 = compartment['name'].tolist()
cp20 = compartment['description'].tolist()
rxn = yeastGEM_s['rxn'].tolist()


def getCompartment(rxn,cp1 = cp10, cp2=cp20):
    cp = [None]*len(cp1)
    for i in range(len(cp1)):
       if cp1[i] in rxn:
         cp[i] = cp2[i]
       else:
          cp[i] = None
    cp1 = [x for i,x in enumerate(cp) if x is not None]
    cp0 = ';'.join(str(e) for e in cp1)
    return cp0


yeastGEM_s['compartment'] = [None]*len(yeastGEM_s)
for j in range(len(rxn)):
    yeastGEM_s['compartment'][j] = getCompartment(rxn[j])


"""input the compartment annotation"""

compartment_SGD['go_type'] = compartment_SGD['go_type'].fillna('NA')
compartment_SGD0 = compartment_SGD[compartment_SGD['go_type'].str.contains("component")]


"""refine the compartment information from SGD"""
ss = compartment_SGD0['systematic_name'].tolist()
yeastGEM_s['gene'] = yeastGEM_s['gene'].str.strip()
ss2 = yeastGEM_s['gene'].tolist()
ss3 = [None]*len(compartment_SGD0)

for i in range(len(compartment_SGD0['go_term'])):
    if ss[i] in ss2:
        ss3[i] = "Yes"
    else:
        ss3[i] = "No"

compartment_SGD0['exist'] = ss3
compartment_SGD1 = compartment_SGD0[compartment_SGD0['exist'].str.contains("Yes")]
#compartment_SGD2 = compartment_SGD1[~compartment_SGD1['go_term'].str.contains("complex")]

unique_GENE_SGD = compartment_SGD1['systematic_name'].unique()
unique_cp_SGD = compartment_SGD1['go_term'].unique()
unique_cp_SGD = unique_cp_SGD.tolist()


"""replace the compartment information using standard name from model"""
xls = pd.ExcelFile('../data/Compartment.xlsx')
sgdChange = pd.read_excel(xls, 'Sheet2')
sgdChange['SGD'] = sgdChange['SGD'].str.strip()
sgdChange['model_name'] = sgdChange['model_name'].str.strip()

compartment_SGD1['compartment'] = singleMapping(sgdChange['model_name'].tolist(),sgdChange['SGD'].tolist(),compartment_SGD1['go_term'].tolist())


gene_compartmentSGD = pd.DataFrame({
    'gene': compartment_SGD1['systematic_name'].unique()
})

gene_compartmentSGD['compartment'] = multiMapping(compartment_SGD1['compartment'].tolist(),compartment_SGD1['systematic_name'].tolist(),gene_compartmentSGD['gene'].tolist())


"""merge the SGD compartment with yeastGEM and data analysis"""
yeastGEM_s['compartment_sgd'] = singleMapping(gene_compartmentSGD['compartment'].tolist(),gene_compartmentSGD['gene'].tolist(),yeastGEM_s['gene'].tolist())
yeastGEM_s0 = yeastGEM_s[(yeastGEM_s['gene'] != 'NA')]






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

#uniprot_location2 = uniprot_location2[~uniprot_location2.V2.isnull()]
unique_compartment_uniprot = pd.Series(uniprot_location2.V2.unique())
unique_compartment_uniprot = unique_compartment_uniprot.str.strip()
unique_compartment_uniprot = unique_compartment_uniprot.unique()
unique_compartment_uniprot = pd.DataFrame(unique_compartment_uniprot)
unique_compartment_uniprot.columns = ['compartment']

def getSimilarTarget(rxn_yeast0,rxn_newGPR0,ss):
    from fuzzywuzzy import fuzz
    from fuzzywuzzy import process
    rxn_yeast1 = np.array(rxn_yeast0)  # np.ndarray()
    rxn_yeast2 = rxn_yeast1.tolist()
    rxn_yeast3 = pd.Series((v[0] for v in rxn_yeast2))
    rxn_newGPR1 = np.array(rxn_newGPR0)  # np.ndarray()
    rxn_newGPR2 = rxn_newGPR1.tolist()
    rxn_newGPR3 = pd.Series((v[0] for v in rxn_newGPR2))
    similarTarget = [None] * ss
    for i in range(ss):
        similarTarget[i] = process.extract(rxn_newGPR3[i], rxn_yeast3, limit=2)

    return similarTarget

ss0 = len(unique_compartment_uniprot)
similarTarget0 = getSimilarTarget(compartment[['description']],unique_compartment_uniprot[['compartment']],ss=ss0)
unique_compartment_uniprot['similar_target'] = similarTarget0
unique_compartment_uniprot.head()

writer = pd.ExcelWriter('../data/unique_compartment_uniprot.xlsx')
unique_compartment_uniprot.to_excel(writer,'Sheet1')
writer.save()

"""replace the compartment information using standard name from model for uniprot gene location annotation"""
xls = pd.ExcelFile('../data/Compartment.xlsx')
uniChange = pd.read_excel(xls, 'Sheet3')
uniChange['uniprot'] = uniChange['uniprot'].str.strip()
uniChange['model_name'] = uniChange['model_name'].str.strip()

compartment_uniprot['model_name'] = singleMapping(uniChange['model_name'].tolist(),uniChange['uniprot'].tolist(),compartment_uniprot['compartment'].tolist())


gene_compartmentUNI = pd.DataFrame({
    'gene': compartment_uniprot['gene'].unique()
})

gene_compartmentUNI['compartment'] = multiMapping(compartment_uniprot['model_name'].tolist(),compartment_uniprot['gene'].tolist(),gene_compartmentUNI['gene'].tolist())



"""merge the uniprot compartment with yeastGEM and data analysis"""
yeastGEM_s['compartment_uni'] = singleMapping(gene_compartmentUNI['compartment'].tolist(),gene_compartmentUNI['gene'].tolist(),yeastGEM_s['gene'].tolist())
yeastGEM_s0 = yeastGEM_s[(yeastGEM_s['gene'] != 'NA')]


"""combine sgd and uniprot"""
yeastGEM_s0.compartment_sgd = yeastGEM_s0.compartment_sgd.astype(str)
yeastGEM_s0.compartment_uni = yeastGEM_s0.compartment_uni.astype(str)
yeastGEM_s0['combine'] = yeastGEM_s0.compartment_sgd.str.cat(yeastGEM_s0.compartment_uni, sep=';')
yeastGEM_s0['combine'] = yeastGEM_s0['combine'].astype(str)
yeastGEM_s0['combine_unique'] = [None]*len(yeastGEM_s0['combine'])


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
def getCommonCompartment(c1,c2):
    c10 = c1.split(";")
    c20 = c2.split(";")
    c3 = list(set(c10).intersection(c20))
    c4 = ';'.join(str(e) for e in c3)
    return c4

sss = [None]*len(yeastGEM_s0['gene'])
for i in range(len(yeastGEM_s0['gene'])):
     sss[i] = getCommonCompartment(yeastGEM_s0['compartment'].tolist()[i],yeastGEM_s0['combine_unique'].tolist()[i])


yeastGEM_s0['common_compartment'] = sss


"""output result"""
writer = pd.ExcelWriter('../result/compartment check result.xlsx')
yeastGEM_s0.to_excel(writer,'Sheet1')
writer.save()
