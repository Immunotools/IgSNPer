import logging

import allele_basic
import annotated_gene

import os
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

class AnnotatedGeneVisualizer:
    def __init__(self, logger):
        self.logger = logger
        self.nucls = {'-' : 0, '*' : 0, 'N' : 0, 'A' : 1, 'C' : 2, 'G' : 3, 'T' : 4}
        self.nucl_heatmap = 'coolwarm'
        self.nucl_min = min(self.nucls.values())
        self.nucl_max = max(self.nucls.values())
        self.num_ind_heatmap = 'summer'
        self.allele_empty_value = 15
        self.allele_max_value = 19
        self.allele_min_value = 0
        self.allele_cmap = 'tab20'

    def _ProcessNucleotideMatrix(self, annot_matrix):
        matrix = []
        for annot_row in annot_matrix:
            row = [self.nucls[n] for n in annot_row]
            matrix.append(row)
        return matrix

    def _Savefig(self, output_base):
        plt.savefig(output_base + '.pdf')
        plt.savefig(output_base + '.svg')


    def VisualizeAlleles(self, gene, output_base):
        """
        @type gene : annotated_gene.AnnotatedGene
        """
        if len(gene.SNPs()) == 0 or gene.NumAlleles() == 1:
            return
#        allowed_alleles = [str(a) for a in [01, 02, 04, 05, 06, 9, 10, 12, 15, 17]]
#        pos_to_remove = []
#        if gene.V() == 'IGHV1-69':
#            pos_to_remove = [3, 23, 44, 215, 218]
        annot_matrix = []
        snps = gene.SNPs()
        snps = [snp for snp in snps if snp[0]] # not in pos_to_remove]
        snp_str = [str(snp[0]) for snp in snps]
        df = {'Descr' : [], 'NumInd' : [], 'Index' : []}
        for a in sorted(gene.AlleleIter(), key = lambda x : gene.NumIndividualsByAllele(x), reverse = True):
#            if gene.V() == 'IGHV1-69' and a.ToShortString() not in allowed_alleles:
#                continue
            annot_row = []
            for snp in snps:
                annot_row.append(a.GetNucleotideByAlignmentPos(snp[0]))
            annot_matrix.append(annot_row)
            df['Descr'].append(a.ToShortString())
            df['NumInd'].append(gene.NumIndividualsByAllele(a))
            df['Index'].append(len(df['Descr']) - 1)
        # output
        fig, axes = plt.subplots(ncols = 2, figsize = (14, 6))
        # heatmap
        matrix = self._ProcessNucleotideMatrix(annot_matrix)
        sns.heatmap(matrix, annot = np.array(annot_matrix), fmt = '', yticklabels = df['Descr'], xticklabels = snp_str,
                    cbar = False, ax = axes[0], cmap = self.nucl_heatmap, vmin = self.nucl_min,
                    vmax = self.nucl_max, linewidth = 0.2)
        plt.sca(axes[0])
        plt.xticks(rotation = 90)
        plt.yticks(rotation = 0, fontsize = 8)
        sns.barplot(y = df['Index'], x = df['NumInd'], orient = 'h', ax = axes[1], palette = self.num_ind_heatmap)
        plt.sca(axes[1])
        plt.yticks([])
        plt.xlabel('# individuals')
        plt.suptitle(gene.V())
        self._Savefig(output_base)
        plt.clf()
        plt.close()

    def _ProcessGenotypeMatrix(self, annot_matrix):
        matrix = []
        for row_ind, row in enumerate(annot_matrix):
            m_row = [0] * len(row)
            matrix.append(m_row)
            for col_ind, e in enumerate(row):
                annot_matrix[row_ind][col_ind] = allele_basic.NuclsToStr(e)
        for j in range(len(matrix[0])):
            row_nucls = [annot_matrix[i][j] for i in range(len(matrix))]
            row_nucls = sorted(row_nucls, key = lambda x : len(x))
            for i in range(len(matrix)):
                if len(annot_matrix[i][j]) == 1:
                    matrix[i][j] = row_nucls.index(annot_matrix[i][j])
                else:
                    matrix[i][j] = 0.5
        return matrix, annot_matrix

    def VisualizeMultiSNPGenotypes(self, gene, output_base):
        """
        @type gene : annotated_gene.AnnotatedGene
        """
        if len(gene.SNPs()) == 0:
            return
        annot_matrix = []
        # SNPs
        pos_to_remove = []
        if gene.V() == 'IGHV1-2':
            pos_to_remove = [254, 266, 293]
        snps = [snp for snp in gene.SNPs() if snp[0] not in pos_to_remove]
        snp_str = [str(snp[0] + 1) for snp in snps]
        df = {'Descr' : [], 'NumInd' : [], 'Index' : []}
        for h_ind, h in enumerate(sorted(gene.HaplotypeIter(), key = lambda x : gene.NumIndividualsByHaplotype(x),
                                         reverse = True)):
            if gene.NumIndividualsByHaplotype(h) == 1:
                break
            annot_row = []
            for snp in snps:
                pos_nucls = h.GetNucleotidesByPosition(snp[0])
                annot_row.append(pos_nucls)
            annot_matrix.append(annot_row)
            df['Descr'].append('G' + str(h_ind + 1))
            df['NumInd'].append(gene.NumIndividualsByHaplotype(h))
            df['Index'].append(len(df['Index']))
        if len(annot_matrix) == 0:
            return
        # output
        fig, axes = plt.subplots(ncols = 2, figsize = (14, 6))
        # heatmap
        matrix, annot_matrix = self._ProcessGenotypeMatrix(annot_matrix)
        sns.heatmap(matrix, annot = np.array(annot_matrix), fmt = '', yticklabels = df['Descr'], xticklabels = snp_str,
                    cbar = False, ax = axes[0], cmap = self.nucl_heatmap, vmin = 0, vmax = 1,
                    annot_kws = {'size' : 6}, linewidth = 0.2)
        plt.sca(axes[0])
        plt.xticks(rotation = 90)
        plt.yticks(rotation = 0)
        sns.barplot(y = df['Index'], x = df['NumInd'], orient = 'h', ax = axes[1], palette = self.num_ind_heatmap)
        plt.sca(axes[1])
        plt.yticks([])
        plt.xlabel('# individuals')
        plt.suptitle(gene.V())
        self._Savefig(output_base)
        plt.clf()
        plt.close()
        # heatmap separately
        plt.figure()
        sns.heatmap(matrix, annot=np.array(annot_matrix), fmt='', yticklabels=df['Descr'], xticklabels=snp_str,
                    cbar=False, cmap=self.nucl_heatmap, vmin=0, vmax=1, annot_kws={'size': 8}, linewidth=0.2)
        plt.yticks(rotation = 0)
        plt.savefig(output_base + '_hmap.pdf')
        plt.clf()
        plt.close()

    def VisualizeHaplotypes(self, gene, output_base):
        if len(gene.SNPs()) == 0:
            return
        annot_matrix = []
        matrix = []
        df = {'NumInd' : [], 'Index' : [], 'Descr' : []}
        haplotypes = [len(h) for h in gene.HaplotypeIter() if gene.NumIndividualsByHaplotype(h) > 1]
        if len(haplotypes) == 0:
            return
        max_num_alleles = max(haplotypes)
        for h_ind, h in enumerate(sorted(gene.HaplotypeIter(), key = lambda x : gene.NumIndividualsByHaplotype(x),
                                         reverse = True)):
            if gene.NumIndividualsByHaplotype(h) == 1:
                break
            row = [self.allele_empty_value] * max_num_alleles
            annot_row = [''] * max_num_alleles
            df['NumInd'].append(gene.NumIndividualsByHaplotype(h))
            df['Index'].append(len(df['Index']))
            df['Descr'].append('G' + str(h_ind + 1) + ' - ' + str(gene.NumIndividualsByHaplotype(h)))
            for a_ind, allele in enumerate(h):
                if allele.Ambiguous():
                    row[a_ind] = self.allele_min_value
                else:
                    row[a_ind] = allele.MainAllele()
                annot_row[a_ind] = str(allele)
            annot_matrix.append(annot_row)
            matrix.append(row)
        plt.figure(figsize = (8, 6))
        sns.heatmap(matrix, annot = np.array(annot_matrix), fmt = '', cmap = self.allele_cmap,
                    vmin = self.allele_min_value, vmax = self.allele_max_value, cbar = False, annot_kws = {'size' : 12},
                    xticklabels = [], yticklabels = df['Descr'], linewidth = 0.2)
        plt.yticks(rotation = 0)
        plt.title(gene.V())
        self._Savefig(output_base)
        plt.clf()
        plt.close()

    def OutputHaplotypesByGene(self, gene, output_fname):
        """
        @type gene : annotated_gene.AnnotatedGene
        """
        fh = open(output_fname, 'w')
        fh.write('Project\tIndividual\tHaplotypeRank\tHaplotypeDescr\n')
        for h_ind, h in enumerate(sorted(gene.HaplotypeIter(), key = lambda x : gene.NumIndividualsByHaplotype(x),
                                         reverse = True)):
            individuals = gene.GetIndividualsByHaplotype(h)
            for pr_id, ind_id in individuals:
                fh.write(pr_id + '\t' + ind_id + '\t' + str(h_ind + 1) + '\t' + str(h) + '\n')
        fh.close()

    def VisualizeSNPStates(self, gene, output_dir):
        """
        @type gene : annotated_gene.AnnotatedGene
        """
        genebase = allele_basic.ModifyGeneName(gene.V())
        for snp_pos in gene.PolymorphismPositionIter():
            nucl_counter = gene.ComputePositionState(snp_pos)
            x = []
            y = []
            for n in sorted(nucl_counter, key = lambda x : nucl_counter[x], reverse = True):
                x.append(n)
                y.append(nucl_counter[n])
            plt.figure()
            sns.barplot(x = x, y = y, palette = self.num_ind_heatmap)
            plt.ylabel('# individuals')
            plt.title(gene.V() + ', position ' + str(snp_pos))
            self._Savefig(os.path.join(output_dir, genebase + '_' + str(snp_pos)))
            plt.clf()
            plt.close()