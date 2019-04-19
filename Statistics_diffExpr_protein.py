#!/usr/bin/env python
#-*- coding:utf-8 -*-
# Author: Todzhu
# Date: 2019/12/20 11:58

import pandas as pd

df_protExp = pd.read_table('MS_identified_information.txt',sep='\t',low_memory=False)

headers = []
compare = {}
for i in df_protExp.columns:
    if '/' not in i:
        headers.append(i)
    else:
        compare[i.split()[0]] = 1

def regulatedType(ratio,pvalue,fold):
    if ratio > fold and pvalue < 0.05:
        return 'Up'
    if ratio < 1/fold and pvalue < 0.05:
        return 'Down'

for fold in [1.2,1.3,1.5,2]:
    writer = pd.ExcelWriter('fold_'+str(fold)+' '+'differentially_expressed_protein.xlsx', engine='xlsxwriter')
    group = []
    up = []
    down = []
    for key in compare:
        headers.extend([key+' Ratio',key+' P value','Regulated Type'])
        df_protExp['Regulated Type'] = df_protExp.apply(lambda x: regulatedType(x[key+' Ratio'],x[key+' P value'],fold),axis=1)
        group.append(key)
        up.append(list(df_protExp['Regulated Type']).count('Up'))
        down.append(list(df_protExp['Regulated Type']).count('Down'))
        df = pd.DataFrame({'Compared sample name': group, 'Up-regulated': up, 'Down-regulated': down})
        df.to_excel(writer, sheet_name='Diff distribution',index=None)
        df_protExp[((df_protExp[key+' Ratio']>fold) | (df_protExp[key+' Ratio']< 1/fold)) & (df_protExp[key+' P value']< 0.05)][headers].to_excel(writer, sheet_name=key.replace('/','vs'),index=None)
        headers.pop()
        headers.pop()
        headers.pop()


    writer.save()

