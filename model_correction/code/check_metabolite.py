import os
import re
import pandas as pd
import numpy as np
# set the directory
os.getcwd()
os.chdir('/Users/luho/PycharmProjects/model/model_correction/code')

met_fulllist = pd.read_table('../data/metabolite_manual_curation_full_list.tsv')
met_fulllist = met_fulllist.iloc[0:1059]
met_part = pd.read_table('../data/metabolite_manual_curation.tsv')

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

met_fulllist['remark2'] = singleMapping(met_part['remark'].tolist(),met_part['metID'].tolist(),met_fulllist['metID'].tolist())

writer = pd.ExcelWriter('../result/metabolite correction.xlsx')
met_fulllist.to_excel(writer,'Sheet1')
writer.save()