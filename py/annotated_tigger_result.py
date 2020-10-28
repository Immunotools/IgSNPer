import sys
import logging

import aligned_alleles
import raw_tigger_results
import annotated_allele_constructor

class AnnotatedTiggerResult:
    def __init__(self, raw_tigger, allele_storage, imgt_numbered_alleles, logger):
        """
        @type raw_tigger: raw_tigger_results.RawTiggerResult
        @type allele_storage: aligned_alleles.AlignedAlleleStorage
        @type imgt_numbered_alleles : aligned_alleles.IMGTNumberedAlleles
        """
        self.raw_tigger = raw_tigger
        self.allele_storage = allele_storage
        self.imgt_numbered_alleles = imgt_numbered_alleles
        self.logger = logger
        self.allele_constructor = annotated_allele_constructor.AnnotatedAlleleConstructor(allele_storage,
                                                                                          imgt_numbered_alleles)
        self._ComputeAnnotatedAlleles()

    def _ComputeAnnotatedAlleles(self):
        self.gene_allele_dict = dict()
        for gene in self.raw_tigger.GeneIter():
            counted_alleles = self.raw_tigger.GetCounterAllelesByGene(gene)
            if len(counted_alleles) == 0:
                continue
            self.gene_allele_dict[gene] = []
            for allele, count in counted_alleles:
                self.gene_allele_dict[gene].append(self.allele_constructor.CreateAnnotateAllele(gene, allele))

    def _CheckGeneFatal(self, gene):
        if gene not in self.gene_allele_dict:
            self.logger.warning("ERROR: gene " + gene + 'is not found in output')
            sys.exit(1)

    def GeneIterator(self):
        for gene in self.gene_allele_dict:
            yield gene

    def GetAnnotatedAllelesByGene(self, gene):
        self._CheckGeneFatal(gene)
        return self.gene_allele_dict[gene]

    def WriteToTxt(self, output_fname):
        fh = open(output_fname, 'w')
        self.logger.info('Annotated TIgGER results were written to ' + output_fname)
        fh.write('Gene\tDescription\tAlleleType\tSequence\tAlignedSequence\tSNPs\tAmbiguousPositions\tIMGTAlleles\t'
                 'NovelModifications\n')
        for gene in self.GeneIterator():
            alleles = self.GetAnnotatedAllelesByGene(gene)
            for a in alleles:
                # SNPs
                snps = [str(snp[0]) + ':' + str(snp[1]) for snp in a.SNPs()]
                snp_str = '-'
                if len(snps) != 0:
                    snp_str = ','.join(snps)
                # allele list
                allele_list = a.AlleleList()
                allele_str = '-'
                if len(allele_list) != 0:
                    allele_str = ','.join([str(a_id) for a_id in allele_list])
                # ambiguous positions
                amb_pos_list = a.AmbiguousPositions()
                amb_pos_str = '-'
                if len(amb_pos_list) != 0:
                    amb_pos_str = ','.join([str(p) for p in amb_pos_list])
                # novel modifications
                modifications = a.MainModification()
                modification_str = '-'
                if len(modifications) != 0:
                    modification_str = ','.join([str(m[0]) + ':' + str(m[1]) for m in modifications])
                fh.write(gene + '\t' + a.Type().name + '\t' + str(a) + '\t' + a.Sequence() + '\t' +
                         a.AlignedSequence() + '\t' + snp_str + '\t' + amb_pos_str + '\t' + allele_str +
                         '\t' + modification_str + '\n')
        fh.close()