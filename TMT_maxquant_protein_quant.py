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

        self.df_summary = pd.read_table(self.summary,sep='\t',low_memory=False)
        self.df_summary = self.df_summary.iloc[-1:]

        self.df_peptides = pd.read_table(self.peptides, sep='\t', low_memory=False)
        self.num_uniq_peptide = list(self.df_peptides['Unique (Groups)']).count('yes')
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
            self.df_prot_out[key+' P value'] = pvalue.apply(lambda x:'%.2E' % x).replace('NAN','')

        self.df_prot_out = pd.merge(self.df_proteinGroups,self.df_prot_out,on='Protein accession',how='right')

        for key, value in self.repeat.items():
            self.df_prot_out[key+' Ratio'] = round(self.df_prot_out.apply(lambda x: x[value[0]].mean()/x[value[1]].mean(), axis=1),3)
            self.df_prot_out[key+' P value'] = self.df_prot_out.apply(lambda x:ttest(x[value[0]].astype(np.float),x[value[1]].astype(np.float)),axis=1).apply(lambda x:'%.2E' % x).replace('NAN','')

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

        merge_format = workbook.add_format({
            'font_name': 'Times New Roman',
            'font_size': 11,
            'bold': 1,
            'border': 1,
            'align': 'center',
            'valign': 'vcenter',
            'fg_color': '00CD00'})

        summary.set_default_row(18)
        summary.set_column('A:A',15.75)
        summary.set_column('B:B', 50.15)
        summary.set_column('D:I', 20.15)
        summary.merge_range('A2:B2','Header line description of table 1',merge_format)
        summary.merge_range('A10:B10', 'Header line description of worksheet "Protein_quant"', merge_format)
        summary.merge_range('A20:B20', 'Header line description of worksheet "Peptide_quant"', merge_format)
        summary.merge_range('D2:I2', 'Table 1. MS/MS spectrum database search analysis summary', merge_format)
        summary.merge_range('D8:I8', 'Table 2. Differentially expressed protein summary (Filtered whith threshold value of expression fold change and P vlaue < 0.05)', merge_format)

        header_description = {'Total spectrum':'Number of spectrum produced by mass spectrometer',
                                           'Matched spectrum':'Number of spectrum matched with alignment protein',
                                           'Peptides':'Number of peptides which spectrum hit',
                                           'Unique peptides':'Number of unique peptides which spectrums hit',
                                           'Identified proteins':'Number of proteins detected by spectrum search analysis',
                                           'Quantifiable proteins':'Number of proteins quantifiable'}

        quant_description = {'Protein accession':'Protein accession number of database used for search',
                                          'Protein description': 'Protein functional description',
                                          'Gene name': 'Indicates the name of the gene that code for the protein sequence',
                                          'MW [kDa]': 'Protein molecular weight ,unit [kDa]',
                                          'Score': 'A simple rule to be used to judge whether a result is significant or not',
                                          'Coverage [%]': 'Percent of identified peptdie sequence covering in protein sequence',
                                          'Peptides': 'Number of peptides whitch spectrums hit',
                                          'Unique peptides': 'Number of unique peptides which spectrums hit'}

        peptide_description = {'ID':'Identified peptide index number',
                                            'Sequence': 'Identified peptide amino acid sequence',
                                            'Protein accession': 'Protein accession number of database used for search',
                                            'Protein description': 'Protein functional description',
                                            'Charge': 'Carried charge of peptide',
                                            'm/z': 'Mass-to-charge ratio of peptide',
                                            'mass error [ppm]': 'Parts per million ratio of peptide mass error between the theory and practice',
                                            'PEP': 'The maximal posterior error probability for peptides',
                                            'Score': 'A simple rule to be used to judge whether a result is significant or not',
                                            'Unique [yes/no]': 'When marked with "yes", this particular peptide is unique to a single identified protein group'}

        format1 = workbook.add_format({
            'bold': True,
            'fg_color': '#EAEAEA',
            'font_size': 11,
            'font_name': 'Times New Roman'})

        format2 = workbook.add_format({
            'fg_color': '#EAEAEA',
            'font_size': 10,
            'font_name': 'Times New Roman'})

        format3 = workbook.add_format({
            'border': 1,
            'valign': 'vcenter',
            'bold': True,
            'align': 'center',
            'fg_color': '#EAEAEA',
            'font_size': 11,
            'font_name': 'Times New Roman'})

        row = 2
        for key, value in header_description.items():
            summary.write(row,0,key,format1)
            summary.write(row,1,value,format2)
            row += 1

        for key, value in quant_description.items():
            summary.write(row+2,0,key,format1)
            summary.write(row+2,1,value,format2)
            row += 1

        for key, value in peptide_description.items():
            summary.write(row+4,0,key,format1)
            summary.write(row+4,1,value,format2)
            row += 1

        summary.write(3, 3, 'Total spectrums', format3)
        summary.write(4, 3, list(self.df_summary['MS/MS'])[0], format3)
        summary.write(3, 4, 'Matched spetrums', format3)
        summary.write(4, 4, list(self.df_summary['MS/MS Identified'])[0], format3)
        summary.write(3, 5, 'Peptides', format3)
        summary.write(4, 5, list(self.df_summary['Peptide Sequences Identified'])[0], format3)
        summary.write(3, 6, 'Unique peptides', format3)
        summary.write(4, 6, self.num_uniq_peptide, format3)
        summary.write(3, 7, 'Identified proteins', format3)
        summary.write(4, 7, self.df_prot_out['Protein accession'].count(), format3)
        summary.write(3, 8, 'Quantifiable proteins', format3)


        def regulatedType(ratio, pvalue, fold):
            if ratio > fold and pvalue < 0.05:
                return 'Up'
            if ratio < 1 / fold and pvalue < 0.05:
                return 'Down'

        l = []
        row = 11
        for key in self.df_prot_out:
            loc = 'D' + str(row) + ':' + 'D' + str(row + 1)
            if ' Ratio' in key:
                l.append(key)
                summary.merge_range(loc,key.replace(' Ratio',''),format3)
                summary.write(row-1, 4 ,'Up-regulated',format3)
                summary.write(row, 4, 'Down-regulated', format3)
                row += 2


        summary.write(4, 8, self.df_prot_out[l].dropna(axis=0, how='all').shape[0], format3)

        summary.write(9,3,'Compare group',format3)
        summary.write(9, 4, 'Regulated type', format3)
        summary.write(9, 5, 'fold change >1.2', format3)
        summary.write(9, 6, 'fold change >1.3', format3)
        summary.write(9, 7, 'fold change >1.5', format3)
        summary.write(9, 8, 'fold change >2', format3)

#
        protein_quant = workbook.add_worksheet('Protein_quant')
        peptide_quant = workbook.add_worksheet('Peptide_quant')
        statistics = workbook.add_worksheet('Statistics')

        protein_quant.set_row(0,25.5)
        protein_quant.set_column('B:B',30)
        peptide_quant.set_column('F:F', 30)
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

        content_format = workbook.add_format({
            'font_size': 10,
            'font_name': 'Times New Roman',
        })


        self.df_prot_out = self.df_prot_out.fillna('')
        n = 0
        for col_name in self.df_prot_out:
            protein_quant.write(0,n,col_name,header_format)
            for index,row in self.df_prot_out.iterrows():
                protein_quant.write(index+1,n,row[col_name],content_format)
            n += 1

        self.df_pep_out = self.df_pep_out.fillna('')
        n = 0
        for col_name in self.df_pep_out:
            peptide_quant.write(0,n,col_name,header_format)
            for index,row in self.df_pep_out.iterrows():
                peptide_quant.write(index+1,n,row[col_name],content_format)
            n += 1

        workbook.close()

if __name__ == '__main__':
    P = proteinQuant('sample', 'evidence.txt', 'peptides.txt', 'proteinGroups.txt', 'summary.txt')
    P.readSample()
    P.peptideQuant()
    P.proteinQuant()
    #P.Statistics()
    P.writeOut()