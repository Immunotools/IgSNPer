import logging

import aligned_alleles
import allele_basic

class AlleleUnifier:
    def __init__(self, allele_storage, logger):
        """
        @type allele_storage: aligned_alleles.AlignedAlleleStorage
        """
        self.allele_storage = allele_storage
        self.logger = logger

    def UnifySNPs(self, v_gene, alleles):
        v_alleles = self.allele_storage.GetAllelesByVGene(allele_basic.ModifyGeneName(v_gene))
        main_allele = v_alleles.MinimalAlleleId()
        snp_position_union = set()
        for allele in alleles:
            snps = allele.SNPs()
            for pos, nucl in snps:
                snp_position_union.add(pos)
        self.logger.debug("Combined position set: " + str(snp_position_union))
        snp_position_union = sorted(list(snp_position_union))
        for a_ind, allele in enumerate(alleles):
            snps = allele.SNPs()
            snp_dict = dict()
            for pos, nucl in snps:
                snp_dict[pos] = nucl
            updated_snps = []
            for pos in snp_position_union:
                if pos in snp_dict:
                    updated_snps.append((pos, snp_dict[pos]))
                else:
                    main_allele_nucl = v_alleles.GetNucleotideByAlleleAndPos(main_allele, pos)
                    updated_snps.append((pos, main_allele_nucl))
            alleles[a_ind].UpdateSNPSet(updated_snps)
        return alleles