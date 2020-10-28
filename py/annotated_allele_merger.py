import logging

import allele_basic
import annotated_alleles
import aligned_alleles

from collections import Counter

class AgressiveAlleleMerger:
    def __init__(self, allele_storage, imgt_numbered_alleles, logger):
        """
        @type allele_storage: aligned_alleles.AlignedAlleleStorage
        @type imgt_numbered_alleles : aligned_alleles.IMGTNumberedAlleles
        """
        self.allele_storage = allele_storage
        self.imgt_numbered_alleles = imgt_numbered_alleles
        self.logger = logger
        # aux fields
        self.glued_pairs = dict() # V gene -> pairs (non amb allele, amb allele)

    def _UpdateGluingStats(self, non_amb_allele, amb_allele):
        v_gene = non_amb_allele.V()
        if v_gene not in self.glued_pairs:
            self.glued_pairs[v_gene] = []
        self.glued_pairs[v_gene].append((non_amb_allele, amb_allele))

    def MergeAlleles(self, annotated_alleles):
        for allele_ind, allele in enumerate(annotated_alleles):
            self.logger.debug(str(allele_ind) + ', ' + allele.ToLongString())
        putative_allele_indices = []
        # select non-ambiguous alleles
        for ind, allele in enumerate(annotated_alleles):
            if allele.Type() not in [allele_basic.AlleleType.AMBIGUOUS_NOVEL, allele_basic.AlleleType.AMBIGUOUS_GERMLINE]:
                putative_allele_indices.append(ind)
        self.logger.debug(putative_allele_indices)
        # remove excessive good alleles
        to_remove = set()
        for ind, allele in enumerate(annotated_alleles):
            if ind in putative_allele_indices:
                continue
            for p_ind in putative_allele_indices:
                self.logger.debug('Gluing check: ' + str(ind) + '>' + str(p_ind))
                can_be_glued = self._AllelesCanBeGlued(annotated_alleles[p_ind], allele)
                if can_be_glued:
                    self.logger.debug('Glued: ' + str(ind) + '>' + str(p_ind))
                    self._UpdateGluingStats(annotated_alleles[p_ind], allele)
                    to_remove.add(ind)
                    break
        # collapse identical sequences
        seq_dict = dict()
        for ind, allele in enumerate(annotated_alleles):
            if ind in to_remove:
                continue
            seq = allele.Sequence()
            if seq not in seq_dict:
                seq_dict[seq] = []
            seq_dict[seq].append(ind)
        # output good alleles
        filtered_annotated_alleles = []
        for s in seq_dict:
            filtered_annotated_alleles.append(annotated_alleles[seq_dict[s][0]])
        self.logger.debug("# good alleles: " + str(len(filtered_annotated_alleles)))
        return filtered_annotated_alleles

    def _AllelesCanBeGlued(self, non_amb_allele, amb_allele):
        snps_1 = non_amb_allele.SNPs()
        snps_2 = amb_allele.SNPs()
        amb_pos = amb_allele.ambiguous_positions
        confirmed_snps = set()
        for snp in snps_2:
            if snp[1] == '*':
                continue
            if snp not in snps_1:
                return False
            confirmed_snps.add(snp)
        for snp in snps_1:
            if snp not in confirmed_snps and snp[0] not in amb_pos:
                return False
        return True

    def PrintGluingStats(self):
        self.logger.debug("Ambiguous alleles...")
        for v in self.glued_pairs:
            self.logger.debug('== ' + v + ' (' + str(len(self.glued_pairs[v])) + ' individuals)')
            counter = Counter(self.glued_pairs[v])
            for a in sorted(counter, key = lambda x : counter[x], reverse = True):
                self.logger.debug('Main allele - ' + str(a[0]) + ', ambiguous allele - ' + str(a[1]) + ': ' + \
                      str(counter[a]) + ' individuals')