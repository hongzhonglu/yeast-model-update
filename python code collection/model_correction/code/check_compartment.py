# Import packages
import numpy as np
import pandas as pd
import os    ##for directory
# set the directory
os.chdir('/Users/luho/PycharmProjects/python learning/venv/project1_modelling/check_compartment/code')
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

def getCommonCompartment(c1,c2):
    c10 = c1.split(";")
    c20 = c2.split(";")
    c3 = list(set(c10).intersection(c20))
    c4 = ';'.join(str(e) for e in c3)
    return c4

sss = [None]*len(yeastGEM_s0['gene'])
for i in range(len(yeastGEM_s0['gene'])):
     sss[i] = getCommonCompartment(yeastGEM_s0['compartment'].tolist()[i],yeastGEM_s0['compartment_sgd'].tolist()[i])


yeastGEM_s0['common_compartment'] = sss


"""output result"""
writer = pd.ExcelWriter('../result/compartment check result.xlsx')
yeastGEM_s0.to_excel(writer,'Sheet1')
writer.save()
















