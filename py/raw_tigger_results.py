import os
import sys
import logging

import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

import allele_basic

class RawTiggerResult:
    def __init__(self, genotype_df, ind_id, logger):
        self.ind_id = ind_id
        self.genotype_df = genotype_df
        self.logger = logger
        self._CreateGeneDict()

    def _CreateGeneDict(self):
        self.gene_dict = dict()
        for i in self.genotype_df.index.values:
            gene_id = self.genotype_df['GENE'][i]
            if gene_id.find('IGHV') == -1:
                continue
#            print('Processing ' + gene_id)
            self.gene_dict[self.genotype_df['GENE'][i]] = allele_basic.CountedGeneAlleles()
            allele_str = str(self.genotype_df['GENOTYPED_ALLELES'][i]) # str(self.genotype_df['ALLELES'][i])
            counts_str = str(self.genotype_df['Freq_by_Clone'][i]) # str(self.genotype_df['COUNTS'][i])
            if allele_str in ['nan', 'Deletion']: # or counts_str == 'nan':
                continue
#            print allele_str, counts_str
            allele_splits = allele_str.split(',')
#            counts_str = ','.join(['100'] * len(allele_splits))
            count_splits = counts_str.split(',') if ',' in counts_str else counts_str.split(';')
            if len(allele_splits) != len(count_splits):
                self.logger.warning("ERROR: fields " + str(allele_splits) + ' & ' + str(count_splits) +
                                    ' are not consistent')
                sys.exit(1)
            for a, c in zip(allele_splits, count_splits):
                count = int(c)
                if count != 0:
                    inferred_allele = allele_basic.InferredAllele()
                    inferred_allele.ParseFromString(a)
                    self.gene_dict[self.genotype_df['GENE'][i]].Add(inferred_allele, count)
#            for a, c in self.gene_dict[self.genotype_df['GENE'][i]]:
#                print(a, a.Type(), c)
#        print('Alleles for ' + str(len(self.gene_dict)) + ' genes were extracted')

    def Individual(self):
        return self.ind_id

    def GetCounterAllelesByGene(self, gene):
        if gene in self.gene_dict:
            return self.gene_dict[gene]
        return allele_basic.CountedGeneAlleles()

    def GeneIter(self):
        for gene in self.gene_dict:
            yield gene

    def Visualize(self, output_dir):
        if len(self.gene_dict) == 0:
            self.logger.warning("Alleles for individual " + str(self.ind_id) + ' are empty')
            return
        max_num_alleles = max([len(self.gene_dict[v]) for v in self.gene_dict])
        matrix = []
        annot_matrix = []
        v_genes = []
        vis_config = allele_basic.AlleleVisConfig()
        for v in sorted(self.gene_dict):
            counts = sum([p[1] for p in self.gene_dict[v]])
            if counts < 20:
                continue
            row = []
            annot_row = []
            alleles = self.gene_dict[v]
            for a, c in alleles:
                color, label = vis_config.GetColorLabelByAllele(a)
                row.append(color)
                annot_row.append(label)
            row.extend([vis_config.GetEmptyAlleleColor()] * (max_num_alleles - len(alleles)))
            annot_row.extend([''] * (max_num_alleles - len(alleles)))
            matrix.append(row)
            annot_matrix.append(annot_row)
            v_genes.append(v)
        sns.heatmap(matrix, annot = np.array(annot_matrix), fmt = '', yticklabels = v_genes,
                    cmap = vis_config.ColorMap(), cbar = False, annot_kws = {'size' : 6},
                    xticklabels = [], vmin = vis_config.MinValue(), vmax = vis_config.MaxValue())
        plt.yticks(fontsize = 6)
        plt.savefig(os.path.join(output_dir, self.ind_id + '_alleles.pdf'))
        plt.clf()

    def VisualizeCounts(self, output_dir):
        x = [] # v position
        y = [] # counts
        counts = []
        v_genes = []
        for v in sorted(self.gene_dict):
            alleles = self.gene_dict[v]
            if alleles.Empty():
                continue
            max_count = max([p[1] for p in alleles])
            sum_count = sum([p[1] for p in alleles])
            for a, c in alleles:
                x.append(len(v_genes))
                y.append(float(c)) # / max_count)
            v_genes.append(v[4 : ])
            counts.append(sum_count)
        fig, axes = plt.subplots(nrows = 2)
        plt.sca(axes[0])
        plt.plot(x, y, marker = 'o', linestyle = 'None', alpha = 0.25)
        plt.xticks(range(len(v_genes)), [''] * len(v_genes))
        plt.ylabel('allele counts')
#        plt.ylabel('allele count / max allele count')
        plt.sca(axes[1])
        plt.bar(range(len(v_genes)), counts)
        plt.xticks(range(len(v_genes)), v_genes, fontsize = 6, rotation = 90)
        plt.ylabel('total gene count')
        plt.savefig(os.path.join(output_dir, self.ind_id + '_counts.pdf'))
        plt.clf()
