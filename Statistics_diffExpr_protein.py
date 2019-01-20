#!/usr/bin/env python
#-*-coding:utf8-*-
# Author: Todzhu
# Date: 2019/1/20 11:58

import xlsxwriter
import pandas as pd

df_protExp = pd.read_table('MS_identified_information.txt',sep='\t',low_memory=False)

headers = []
compare = []
sheet = {}
for i in df_protExp.columns:
    if '/' not in i:
        headers.append(i)
    else:
        compare.append(i)
        sheet[i.split()[0].replace('/','vs')] = 1

for fold in [1.2,1.3,1.5,2]:
    workbook = xlsxwriter.Workbook(str(fold)+'_fold_differentially_expressed_protein.xlsx')
    workbook.add_worksheet('Differentially')
    for key in sheet:
        workbook.add_worksheet(key)


















