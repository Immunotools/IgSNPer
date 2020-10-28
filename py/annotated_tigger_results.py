import logging

import aligned_alleles
import annotated_tigger_result
import annotated_allele_merger
import annotated_gene
import annotated_allele_unifier

class AnnotatedTiggerResults:
    def __init__(self, allele_storage, imgt_numbered_alleles, allele_merger, allele_unifier, logger):
        """
        @type allele_storage: aligned_alleles.AlignedAlleleStorage
        @type imgt_numbered_alleles : aligned_alleles.IMGTNumberedAlleles
        @type allele_merger : annotated_allele_merger.AgressiveAlleleMerger
        @type allele_unifier : annotated_allele_unifier.AlleleUnifier
        """
        self.allele_storage = allele_storage
        self.imgt_numbered_alleles = imgt_numbered_alleles
        self.allele_merger = allele_merger
        self.allele_unifier = allele_unifier
        self.logger = logger
        self.tigger_dict = dict() # project_id, individual_id -> annotated_tigger_result.AnnotatedTiggerResult
        self.gene_dict = dict() # gene -> list of annotated_gene.AnnotatedGene

    def AddAllele(self, project_id, ind_id, tigger_result):
        """
        @type tigger_result : annotated_tigger_result.AnnotatedTiggerResult
        """
        individual_id = (project_id, ind_id)
        self.logger.info('Processing individual ' + str(individual_id))
        self.tigger_dict[individual_id] = tigger_result
        for gene in tigger_result.GeneIterator():
            alleles = tigger_result.GetAnnotatedAllelesByGene(gene)
            self.logger.debug('Processing ' + gene + ', individual ' + str(individual_id) +
                         ', # alleles: ' + str(len(alleles)))
            alleles = self.allele_unifier.UnifySNPs(gene, alleles)
            merged_alleles = self.allele_merger.MergeAlleles(alleles)
            if gene not in self.gene_dict:
                self.gene_dict[gene] = annotated_gene.AnnotatedGene(gene, self.allele_storage, self.allele_unifier)
            self.gene_dict[gene].AddIndividualAlleles(individual_id, merged_alleles)

    def ComputeAllelesByGene(self):
        return

    def GeneIterator(self):
        for gene in self.gene_dict:
            yield self.gene_dict[gene]

    def PrintGluingStats(self):
        self.allele_merger.PrintGluingStats()