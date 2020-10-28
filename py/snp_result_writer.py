import allele_basic
import aligned_alleles
import annotated_gene
import annotated_tigger_results

class SNPResultWriter:
    def __init__(self, allele_storage):
        """
        @type allele_storage: aligned_alleles.AlignedAlleleStorage
        """
        self.allele_storage = allele_storage

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
        fh.write('Dataset\tIndividual\tGene\tHaplotype\tPos\tState\tAllelePos\n')
        for annot_gene in sorted(annotated_results.GeneIterator(), key = lambda x : x.V()):
            if not annot_gene.ContainsIndividual(ind_id):
                continue
            v_name = allele_basic.ModifyGeneName(annot_gene.V())
            aligned_alleles = self.allele_storage.GetAllelesByVGene(v_name)
            main_allele_id = aligned_alleles.MinimalAlleleId()
            ind_haplotype = annot_gene.GetHaplotypeByIndividual(ind_id)
            for snp_pos in annot_gene.PolymorphismPositionIter():
                snp_state = annot_gene.GetSNPsByIndividual(ind_id, snp_pos)
                fh.write(ind_id[0] + '\t' + ind_id[1] + '\t' + annot_gene.V() + '\t' + str(ind_haplotype) + '\t' +
                         str(snp_pos) + '\t' + snp_state + '\t' + str(main_allele_id) + '\n')
        fh.close()