import os

import allele_basic
import annotated_gene
import annotated_tigger_results

class SNPResultWriter:
    def __init__(self):
        return

    def OutputSNPStates(self, gene, output_fname):
        """
        @type gene : annotated_gene.AnnotatedGene
        """
        fh = open(output_fname, 'w')
        fh.write('Gene\tPos\tState\tDataset\tIndividual\n')
        for snp_pos in gene.PolymorphismPositionIter():
            pos_states = gene.ComputeIndividualsByPositionState(snp_pos)
            for st in pos_states:
                for ind in pos_states[st]:
                    fh.write(gene.V() + '\t' + str(snp_pos) + '\t' + st + '\t' + ind[0] + '\t' + ind[1] + '\n')
        fh.close()

    def OutputSNPsPerSubject(self, annotated_results, ind_id, output_fname):
        """
        @type annotated_results : annotated_tigger_results.AnnotatedTiggerResults
        """
        fh = open(output_fname, 'w')
        fh.write('Dataset\tIndividual\tGene\tPos\tState\n')
        for annot_gene in sorted(annotated_results.GeneIterator(), key = lambda x : x.V()):
            if not annot_gene.ContainsIndividual(ind_id):
                continue
            for snp_pos in annot_gene.PolymorphismPositionIter():
                snp_state = annot_gene.GetSNPsByIndividual(ind_id, snp_pos)
                fh.write(ind_id[0] + '\t' + ind_id[1] + '\t' + annot_gene.V() + '\t' + str(snp_pos) + '\t' +
                         snp_state + '\n')
        fh.close()