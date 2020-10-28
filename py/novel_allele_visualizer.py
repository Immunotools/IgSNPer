import allele_basic
import aligned_alleles
import annotated_gene
import annotated_tigger_results

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

class NovelAlleleVisualizer:
    def __init__(self, annotated_results, allele_storage):
        """
        @type annotated_results : annotated_tigger_results.AnnotatedTiggerResults
        @type allele_storage : aligned_alleles.AlignedAlleleStorage
        """
        self.annotated_results = annotated_results
        self.allele_storage = allele_storage
        self._ProcessNovelAlleles()
        self.bar_color = '#4987DD'

    def _ProcessNovelAlleles(self):
        """
        @type gene : annotated_gene.AnnotatedGene
        """
        self.df = {'VGene' : [], 'Position' : [], 'DistToEnd' : [], 'AlleleID' : []}
        for gene in self.annotated_results.GeneIterator():
            vgene = gene.V()
            v_aligned_alleles = self.allele_storage.GetAllelesByVGene(allele_basic.ModifyGeneName(vgene))
            alignment_len = v_aligned_alleles.AlignmentLength()
            for a_ind, a in enumerate(gene.AlleleIter()):
                modifications = a.MainModification()
                for m in modifications:
                    self.df['VGene'].append(vgene[4 : ])
                    self.df['Position'].append(m[0])
                    self.df['DistToEnd'].append(alignment_len - m[0])
                    self.df['AlleleID'].append(vgene + '_' + str(a_ind))
        self.df = pd.DataFrame(self.df)

    def OutputAllelePositions(self, output_fname):
        sorted_pos = sorted(set(self.df['DistToEnd']))
        sns.countplot(x = 'DistToEnd', data = self.df, order = sorted_pos, color = self.bar_color)
        plt.xticks(rotation = 90, fontsize = 7)
        plt.xlabel('distance to the gene end (nt)')
        plt.ylabel('# novel mismatches')
        plt.savefig(output_fname)
        plt.clf()
        #
#        print 3, len(self.df.loc[self.df['DistToEnd'] <= 3]), '/', len(self.df)
#        print 10, len(self.df.loc[self.df['DistToEnd'] <= 10]), '/', len(self.df)

    def OutputPutativePositions(self, output_fname):
        put_df = self.df.loc[self.df['DistToEnd'] <= 10]
        allele_df = put_df.groupby(['VGene', 'AlleleID']).agg({'Position' : 'count'})
        allele_df.reset_index(inplace=True)
        v_genes = sorted(set(allele_df['VGene']))
        sns.countplot(x = 'VGene', data = allele_df, order = v_genes, color = self.bar_color)
        plt.xticks(rotation = 90, fontsize = 8)
        plt.ylabel('# novel alleles')
        plt.xlabel('')
        plt.savefig(output_fname)
        plt.clf()
#        print len(allele_df), len(allele_df.loc[allele_df['VGene'] == '3-23']), len(set(allele_df['VGene']))