import os
import re
import pandas as pd
import numpy as np
# set the directory
os.chdir('/Users/luho/PycharmProjects/python learning/venv/project1_modelling/check_compartment/code')
os.getcwd()
subsystem = pd.read_table('../data/subsystem_correction.tsv')

subsystem['Subsystem_new'] = subsystem['Subsystem_new'].str.replace('(', '@')  # string replace
subsystem['num'] = [None]*len(subsystem['Subsystem_new']) # construct the null column
for i in range(len(subsystem['Subsystem_new'])):
    subsystem['num'][i] = subsystem['Subsystem_new'][i].count('@') # count the number of specific character in each row

subsystem1 = subsystem[subsystem['num'] >=2]

subsystem['Subsystem_new'] = subsystem['Subsystem_new'].str.replace('@gpi', 'gpi')
subsystem['Subsystem_new'] = subsystem['Subsystem_new'].str.replace('@tca', 'gpi')

subsystem2 = subsystem['Subsystem_new'].str.split('@', expand=True) # split one column into multiple
subsystem2.iloc[:,1] = subsystem2.iloc[:,1].str.replace(')','')
subsystem2['new'] = [None]*len(subsystem['Subsystem_new'])
subsystem2['new'] = subsystem2.iloc[:,1] + ' ' + subsystem2.iloc[:,0].map(str) # combine multiple column into one


'''replace nan with subysystem'''
for i in range(len(subsystem['Subsystem_new'])):
    if pd.isnull(subsystem2['new'][i]):
        subsystem2['new'][i] = subsystem2.iloc[:,0][i]
    else:
        subsystem2['new'][i] = subsystem2['new'][i]

subsystem['Subsystem_new0'] = subsystem2['new']


'''correct the remark part'''
subsystem3 = subsystem.iloc[:,2].str.split(';', expand=True)

subsystem['removed_duplicate_subsystem'] = subsystem3.iloc[:,0]
subsystem['evidence'] = subsystem3.iloc[:,1] + ';' + subsystem3.iloc[:,2].map(str)
subsystem['evidence'] = subsystem['evidence'].str.replace(';None','')

for i in range(len(subsystem['Subsystem_new'])):
    if pd.isnull(subsystem['evidence'][i]):
        subsystem['evidence'][i] = ''
    else:
        subsystem['evidence'][i] = subsystem['evidence'][i]

ss0 = subsystem['evidence'].tolist()
z1 = [i for i, s in enumerate(ss0) if 'sce' in s]
for i in z1:
    subsystem['Subsystem_new0'][i] = subsystem['Subsystem_new0'][i] + subsystem['evidence'][i]
    subsystem['evidence'][i] = ''

subsystem['removed_duplicate_subsystem'] = subsystem['removed_duplicate_subsystem'].str.replace(r'/', ' ')
subsystem['removed_duplicate_subsystem'] = subsystem['removed_duplicate_subsystem'].str.strip() # remove the space



ss1 = subsystem['removed_duplicate_subsystem'].tolist()
z2 = [i for i, s in enumerate(ss1) if 'sce' not in str(s)]
for i in z2:
    subsystem['evidence'][i] = subsystem['removed_duplicate_subsystem'][i]
    subsystem['removed_duplicate_subsystem'][i] = ''

subsystem['evidence'] = subsystem['evidence'].str.strip()
writer = pd.ExcelWriter('../result/subsystem.xlsx')
subsystem.to_excel(writer,'Sheet1')
writer.save()