#!/usr/bin/env python
#-*-coding:utf8-*-
# author: Todzhu
# Date: 2019/1/14 10:26

import numpy as np
import pandas as pd
from scipy import stats
from plotnine import *
import xlsxwriter


class proteinQuant:
    def __init__(self, sample, evidence, peptides, proteinGroups, summary):
        self.sample = sample
        self.evidence = evidence
        self.peptides = peptides
        self.proteinGroups = proteinGroups
        self.summary = summary
        self.label = {}
        self.compare = {}
        self.repeat = {}
        self.coefficient = {}
        self.df_peptides = pd.read_table(self.peptides, sep='\t', low_memory=False)
        self.df_peptides = self.df_peptides[(self.df_peptides['Reverse'] != '+') & (self.df_peptides['Potential contaminant'] != '+')][['Sequence','Unique (Groups)']]
        self.df_evidence = pd.read_table(self.evidence, sep='\t', low_memory=False)
        self.df_evidence = self.df_evidence[(self.df_evidence['Reverse'] != '+') & (self.df_evidence['Potential contaminant'] != '+')]
        self.df_proteinGroups = pd.read_table(self.proteinGroups,sep='\t',low_memory=False)
        self.df_proteinGroups = self.df_proteinGroups[(self.df_proteinGroups['Reverse'] != '+') & (self.df_proteinGroups['Potential contaminant'] != '+')]


    def readSample(self):

        with open(self.sample, 'r') as sample:
            for line in sample.readlines():
                line = line.strip().split('\t')
                if '/' not in line[0]:
                    self.label['Reporter intensity ' + line[0]] = line[1]
                elif len(line[0].split('/')[0])>2 and len(line[0].split('/')[0])>2:
                    treat = [self.label['Reporter intensity '+x] for l in line[0].split('/')[0] for x in l]
                    contr = [self.label['Reporter intensity '+x] for l in line[0].split('/')[1] for x in l]
                    self.repeat.setdefault(line[1],[]).append(treat)
                    self.repeat.setdefault(line[1], []).append(contr)
                else:
                    self.compare[line[1]] = line[1]


    def peptideQuant(self):

        pep_out_header = ['id','Sequence','Proteins','Leading Razor Protein','Gene Names','Protein Names','Charge','m/z','Uncalibrated - Calibrated m/z [ppm]','PEP','Score','Unique (Groups)']
        pep_out_header.extend(list(self.label.keys()))
        df_pep_out = pd.merge(self.df_evidence,self.df_peptides,on='Sequence',how='outer')
        df_pep_out = df_pep_out[pep_out_header]
        df_pep_out[list(self.label.keys())] = df_pep_out[list(self.label.keys())].replace(0,np.nan)
        df_pep_out[list(self.label.keys())] = df_pep_out[list(self.label.keys())].apply(lambda x: x/np.mean(x),axis=1)
        #df_pep_out[list(self.label.keys())] = df_pep_out[list(self.label.keys())].apply(lambda x: round(x,3),axis=1).replace(np.nan,'NAN')
        pep_out_header_new = {'id':'ID','Leading Razor Protein':'Protein accession','Gene Names':'Gene name','Protein Names':'Protein description','Uncalibrated - Calibrated m/z [ppm]':'mass error [ppm]','Unique (Groups)':'Unique [yes/no]'}
        pep_out_header_new.update(self.label)
        df_pep_out.rename(columns=pep_out_header_new, inplace=True)
        self.df_pep_out = df_pep_out
        self.df_prot = df_pep_out


    def proteinQuant(self):

        self.df_proteinGroups = self.df_proteinGroups.replace(r';.*','',regex=True)
        self.df_proteinGroups = self.df_proteinGroups[['Majority protein IDs','Protein names','Gene names','Mol. weight [kDa]','Score','Sequence coverage [%]','Peptides','Unique peptides']]
        self.df_proteinGroups.rename(columns={'Majority protein IDs':'Protein accession','Protein names':'Protein description','Gene names':'Gene name','Mol. weight [kDa]':'MW [kDa]','Sequence coverage [%]':'Coverage [%]'},inplace=True)
        # Normalization coefficient
        for i in self.label.values():
            self.coefficient[i] = self.df_prot[i].median()
            self.df_prot[i] = round(self.df_prot[i] / self.coefficient[i], 3)
            #print(i,self.coefficient[i])

        group = self.df_prot[self.df_prot['Unique [yes/no]']=='yes'].groupby('Protein accession')
        self.df_prot_out = round(group[list(self.label.values())].median(),3)
        self.df_prot_out.insert(0,'PSMs',self.df_prot.groupby('Protein accession').size())

        def ttest(arr1,arr2):
            arr1 = arr1.dropna()
            arr2 = arr2.dropna()
            if arr1.size > 2  and arr2.size > 2 and stats.levene(arr1,arr2)[1] < 0.05:
                pvalue = stats.ttest_ind(np.log2(arr1),np.log2(arr2),equal_var=False)
                return pvalue[1]
            elif arr1.size > 1  and arr2.size > 1:
                pvalue = stats.ttest_ind(np.log2(arr1), np.log2(arr2), equal_var=True)
                return pvalue[1]

        for key in self.compare.values():
            fenzi = key.split('/')[0]
            fenmu = key.split('/')[1]
            ratio = self.df_prot.groupby(['Protein accession']).apply(lambda x: x[fenzi].dropna().median()/x[fenmu].dropna().median())
            pvalue = self.df_prot.groupby(['Protein accession']).apply(lambda x: ttest(x[fenzi], x[fenmu]))
            self.df_prot_out[key+' Ratio'] = round(ratio,3)
            self.df_prot_out[key+' P value'] = pvalue.apply(lambda x:'%e' % x).replace('nan','')

        self.df_prot_out = pd.merge(self.df_proteinGroups,self.df_prot_out,on='Protein accession',how='right')

        for key, value in self.repeat.items():
            self.df_prot_out[key+' Ratio'] = round(self.df_prot_out.apply(lambda x: x[value[0]].mean()/x[value[1]].mean(), axis=1),3)
            self.df_prot_out[key+' P value'] = self.df_prot_out.apply(lambda x:ttest(x[value[0]].astype(np.float),x[value[1]].astype(np.float)),axis=1).apply(lambda x:'%e' % x).replace('nan','')

        self.df_prot_out.to_csv('MS_identified_information.txt',sep='\t',index=None)

    def Statistics(self):

        print(ggplot(self.df_prot_out)
              + aes('MW [kDa]','Coverage [%]')
              + geom_point(colour='#1C86EE', size=4,alpha=0.5)
              + xlim(0, 500)
              + ylim(0, 100)
              + theme_linedraw()
              + labs(x="Protein mass [kDa]", y="Protein sequence coverage [%]", title="Protein mass and coverage distribution")
              )


    def writeOut(self):
        workbook = xlsxwriter.Workbook('MS_identified_information.xlsx')
        summary = workbook.add_worksheet('Summary')
        protein_quant = workbook.add_worksheet('Protein_quant')
        peptide_quant = workbook.add_worksheet('Peptide_quant')
        statistics = workbook.add_worksheet('Statistics')

        protein_quant.set_row(0,25.5)
        protein_quant.set_column('B:B',30)
        header_format = workbook.add_format({
            'text_wrap': 1,   #自动换行
            'font_size':10,
            'font_name':'Times New Roman',
            'align':'center',
            'valign':'vcenter',
            'bold': True,
            'fg_color': '#00CD00',
            'border':1,
            'bottom':1,
            'border_color':'#000000'
            })

        for col_num,value in enumerate(self.df_prot_out.columns):
            protein_quant.write(0,col_num,value,header_format)

        workbook.close()

if __name__ == '__main__':
    P = proteinQuant('sample', 'evidence.txt', 'peptides.txt', 'proteinGroups.txt', 'summary.txt')
    P.readSample()
    P.peptideQuant()
    P.proteinQuant()
    P.Statistics()
    P.writeOut()