# -*- coding: utf-8 -*-
# -*- python 3 -*-
# -*- hongzhong Lu -*-
import os
import sys
os.chdir('/Users/luho/PycharmProjects/model/model_correction/code')
sys.path.append(r"/Users/luho/PycharmProjects/model/cobrapy/code")

# import self function
from mainFunction import *

#refine the subsystem for the yeast map based on the reaction number and manual check results
subsystem_map_v2 = pd.read_excel('../result/subsystem for yeast8 map_V2.xlsx')
subsystem0=subsystem_map_v2.loc[:,'subsystem_map_v2'].tolist()
subsystem1=list()
for i, x in enumerate(subsystem0):
    x0=x.strip()
    if x0.startswith('sce'):
        print(x0)
        s= x0.split(' ')
        print(s)
        s1= [s0 for i, s0 in enumerate(s) if i !=0 and s0 !='']
        ss= s[0]
        print(s1)
        s2=' '.join(s1) + ' ( '  + ss + ' )'
        print(s2)
        subsystem1.append(s2)

    else:
        subsystem1.append(x)

#quality check
subsystem1 = [x.lower() for x in subsystem1]
subsystem2=list()
for i, x in enumerate(subsystem1):
    x0=x.strip()
    if 'sce' in x0:
        print(x0)
        s= x0.split(' ( ')
        print(s)
        subsystem2.append(s[0])
    else:
        subsystem2.append(x)

#remove the transport
subsystem3 = [x for x in subsystem2 if 'transport' not in x]

#count the reactions contained in each subsystem
#subsystem_summary = dict((i, subsystem2.count(i)) for i in subsystem3)
#here we only choose subsytem with reaction number larger than 3, if smaller than three, it will be classified into 'other'
#or classified into the related subsystem
#lysine biosynthesis  ==> lysine metabolism
#lysine degradation   ==> lysine metabolism
#fatty acid elongation ( sce00062 ) ==> biosynthesis of unsaturated fatty acids ( sce01040 )
#peroxisome ==> other
subsystem_summary= pd.Series(subsystem3).value_counts()
subsystem_summary[subsystem_summary >=3]



#save the result
subsystem_map_v2['subsystem_map_v3'] = subsystem2

# adjust the format of the removed subsystem
subsystem_map_v2['removed_duplicate_subsystem'] = subsystem_map_v2['removed_duplicate_subsystem'].str.replace('sce','// sce').str.strip()

# find the standard ID for all the subsystem based on kegg annotation
pathway = pd.read_excel('../data/sce_kegg.xlsx', sheet_name='pathwayList')
pathway['name']=pathway['pathwayName'].str.lower()
subsystem_map_v2['pathwayID'] = singleMapping(pathway['pathwayID'],pathway['name'],subsystem_map_v2['subsystem_map_v3'])
for i in range(0,len(subsystem2)):
    if subsystem_map_v2.loc[i,'pathwayID'] is None:
        subsystem_map_v2.loc[i, 'subsystem_map_v3'] = subsystem_map_v2.loc[i, 'subsystem_map_v3']
    else:
        subsystem_map_v2.loc[i, 'subsystem_map_v3'] = subsystem_map_v2.loc[i, 'subsystem_map_v3'] +' ( ' + subsystem_map_v2.loc[i,'pathwayID'] + ' )'

for i in range(0,len(subsystem2)):
    if  'glycerolipid metabolism /  glycerophospholipid metabolism' not in subsystem_map_v2.loc[i, 'subsystem_map_v3']:
        subsystem_map_v2.loc[i, 'subsystem_map_v3'] = subsystem_map_v2.loc[i, 'subsystem_map_v3']
    else:
        subsystem_map_v2.loc[i, 'subsystem_map_v3'] = subsystem_map_v2.loc[i, 'subsystem_map_v3'] + ' ( sce00561,sce00564 )'
        print(subsystem_map_v2.loc[i, 'subsystem_map_v3'])

saveExcel(subsystem_map_v2,"../result/subsystem for yeast8 map_V3.xlsx")








































