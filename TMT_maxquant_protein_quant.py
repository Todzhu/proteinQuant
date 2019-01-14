#!/usr/bin/env python
#-*-coding:utf8-*-
# author: Todzhu
# Date: 2019/1/14 10:26

import numpy as np
import pandas as pd
from scipy import stats
from plotnine import *


class proteinQuant:
    def __init__(self, sample, evidence, peptides, proteinGroups, summary):
        self.sample = sample
        self.evidence = evidence
        self.peptides = peptides
        self.proteinGroups = proteinGroups
        self.summary = summary
        self.label = {}
        self.compare = {}
        self.coefficient = {}

        self.df_peptides = pd.read_table(self.peptides, sep='\t', low_memory=False)
        self.df_peptides = self.df_peptides[(self.df_peptides['Reverse'] != '+') & (self.df_peptides['Potential contaminant'] != '+')][['Sequence','Unique (Groups)']]

        self.df_evidence = pd.read_table(self.evidence, sep='\t', low_memory=False)
        self.df_evidence = self.df_evidence[(self.df_evidence['Reverse'] != '+') & (self.df_evidence['Potential contaminant'] != '+')]


    def readSample(self):
        with open(self.sample, 'r') as sample:
            for line in sample.readlines():
                line = line.strip().split('\t')
                if '/' not in line[0]:
                    self.label['Reporter intensity ' + line[0]] = line[1]
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

        self.df_prot = df_pep_out


    def proteinQuant(self):
        # Normalization coefficient
        for i in self.label.values():
            self.coefficient[i] = self.df_prot[i].median()
            self.df_prot[i] = round(self.df_prot[i] / self.coefficient[i],3)

        def divid(arr1,arr2):
            arr1 = arr1.dropna()
            arr2 = arr2.dropna()
            return arr1.median()/arr2.median()

        def ttest(arr1,arr2):
            arr1 = arr1.dropna()
            arr2 = arr2.dropna()
            if arr1.size > 2  and arr2.size > 2 and stats.levene(arr1,arr2)[1] < 0.05:
                pvalue = stats.ttest_ind(arr1,arr2,equal_var=False)
                return pvalue[1]
            elif arr1.size > 1  and arr2.size > 1:
                pvalue = stats.ttest_ind(arr1, arr2, equal_var=True)
                return pvalue[1]

        group = self.df_prot[self.df_prot['Unique [yes/no]']=='yes'].groupby('Protein accession')
        self.df_prot_out = round(group[list(self.label.values())].median(),3)

        for key in self.compare.values():
            fenzi = key.split('/')[0]
            fenmu = key.split('/')[1]
            ratio = self.df_prot.groupby(['Protein accession']).apply(lambda x: divid(x[fenzi], x[fenmu]))
            pvalue = self.df_prot.groupby(['Protein accession']).apply(lambda x: ttest(np.log2(x[fenzi]), np.log2(x[fenmu])))
            self.df_prot_out[key+' Ratio'] = ratio
            self.df_prot_out[key+' P value'] = pvalue

        self.df_prot_out.to_csv('test.xls',sep='\t')



    def Statistics(self):
        pass


if __name__ == '__main__':
    P = proteinQuant('sample', 'evidence.txt', 'peptides.txt', 'proteinGroups.txt', 'summary.txt')
    P.readSample()
    P.peptideQuant()
    P.proteinQuant()
