import allele_basic
import aligned_alleles
import annotated_alleles

class AnnotatedAlleleConstructor:
    def __init__(self, allele_storage, imgt_numbered_alleles):
        """
        @type allele_storage: aligned_alleles.AlignedAlleleStorage
        @type imgt_numbered_alleles : aligned_alleles.IMGTNumberedAlleles
        """
        self.allele_storage = allele_storage
        self.imgt_numbered_alleles = imgt_numbered_alleles

    def CreateAnnotateAllele(self, v_gene, inferred_allele):
        """
        @type inferred_allele : InferredAllele
        """
        modified_v_gene = allele_basic.ModifyGeneName(v_gene)
        if inferred_allele.Ambiguous():
            return self._CreateAnnotatedAlleleFromAmbiguousAllele(v_gene, inferred_allele)
        aligned_imgt_alleles = self.allele_storage.GetAllelesByVGene(modified_v_gene)
        main_allele_id = aligned_imgt_alleles.MinimalAlleleId()
        inferred_allele_id = inferred_allele.MainAllele()
        inferred_imgt_allele = self.imgt_numbered_alleles.GetIMGTAllele(v_gene, inferred_allele_id)
        known_snps = []
        if len(aligned_imgt_alleles) != 0:
            known_snps = [(s[0], s[2]) for s in
                        aligned_imgt_alleles.GetDifferencesBetweenAlleles(main_allele_id, inferred_allele_id)]
        # add modifications
        novel_descr, novel_snps = self._GetNovelModifications(inferred_imgt_allele, inferred_allele_id,
                                                                    aligned_imgt_alleles,
                                                                    inferred_allele.MainModification())
        # modifying allele sequence
        aligned_seq_list = list(aligned_imgt_alleles.GetSequenceByAllele(main_allele_id))
        snp_list = self._MergeSNPLists(known_snps, novel_snps)
        modified_aligned_seq = self._CreateSequenceBySNPList(aligned_seq_list, snp_list, [])
        # type and description
        allele_descr = str(inferred_allele_id)
        if len(novel_descr) != 0:
            allele_descr += '_' + ','.join(novel_descr)
        allele_type = inferred_allele.Type()
#        print '===='
#        print "V gene: " + v_gene + ", allele: " + str(inferred_allele)
#        print "IMGT SNPs: " + str(known_snps)
#        print "Combined SNP: " + str(snp_list)
#        print 'Original aligned sequence: ' + aligned_imgt_alleles.GetSequenceByAllele(main_allele_id)
#        print 'Modified aligned sequence: ' + modified_aligned_seq
        return annotated_alleles.AnnotatedVGeneAllele(allele_type, v_gene, allele_descr, snp_list, [], modified_aligned_seq,
                                    [inferred_allele_id], novel_snps)

    def _ComputeSNPList(self, allele_ids, aligned_imgt_alleles, main_allele_id):
        snp_list = []
        ambiguous_positions = []
        for p in aligned_imgt_alleles.PolymoprhismIter():
            min_allele_nucl = aligned_imgt_alleles.GetNucleotideByAlleleAndPos(main_allele_id, p)
            amb_nucls = dict()
            for allele_id in allele_ids:
                allele_nucl = aligned_imgt_alleles.GetNucleotideByAlleleAndPos(allele_id, p)
                if allele_nucl not in amb_nucls:
                    amb_nucls[allele_nucl] = 0
                amb_nucls[allele_nucl] += 1
            final_nucl = list(amb_nucls.keys())[0]
            if len(amb_nucls) == 1 and final_nucl != min_allele_nucl:
                snp_list.append((p, final_nucl))
            elif len(amb_nucls) > 1:
                snp_list.append((p, '*'))
                ambiguous_positions.append(p)
        return snp_list, ambiguous_positions

    def _CreateSequenceBySNPList(self, aligned_seq_list, snp_set, ambiguous_positions):
        for snp in snp_set:
            aligned_seq_list[snp[0]] = snp[1]
        for p in ambiguous_positions:
            aligned_seq_list[p] = '*'
        return ''.join(aligned_seq_list)

    def _GetNovelModifications(self, imgt_allele, allele_id, aligned_imgt_alleles,
                                     inferred_modifications):
        novel_descr = []
        novel_snps = []
        for m in inferred_modifications:
            allele_pos = imgt_allele.GetPositionByIMGTPosition(m.Pos())
            aligned_pos = aligned_imgt_alleles.GetAlignmentPosByAllelePos(allele_id, allele_pos)
            novel_snps.append((aligned_pos, m.DstNucl()))
            novel_descr.append(str(aligned_pos) + ':' + m.DstNucl())
        return novel_descr, novel_snps

    def _MergeSNPLists(self, known_snps, novel_snps):
        known_snps.extend(novel_snps)
        merged_snps = []
        for snp in sorted(known_snps, key = lambda x : x[0]):
            merged_snps.append(snp)
        return merged_snps

    def _GetAlleleAndDescription(self, inferred_allele, v_gene, ambiguous_positions, allele_ids, snp_list, novel_descr):
        allele_type = allele_basic.AlleleType.AMBIGUOUS_GERMLINE
        allele_descr = ','.join([str(a) for a in sorted(allele_ids)])
        allele_list = []
        if len(ambiguous_positions) == 0: # allele is not ambiguous
            closest_allele_id = self._FindClosestAllele(v_gene, allele_ids, snp_list)
            allele_descr = str(closest_allele_id)
            allele_list.append(closest_allele_id)
            if len(novel_descr) == 0: # allele is known
                allele_type = allele_basic.AlleleType.GERMLINE
            else: # allele is novel
                allele_descr += '_' + ','.join(novel_descr)
                allele_type = allele_basic.AlleleType.NOVEL
        else: # allele is ambiguous
            allele_list = inferred_allele.AlleleList()
            if len(novel_descr) != 0:
                allele_descr += '_' + ','.join(novel_descr)
                allele_type = allele_basic.AlleleType.AMBIGUOUS_NOVEL
#        print "Type: " + allele_type.name
        return allele_type, allele_descr, allele_list

    def _CreateAnnotatedAlleleFromAmbiguousAllele(self, v_gene, inferred_allele):
        """
        @type inferred_allele : InferredAllele
        """
        allele_ids = [allele_id for allele_id in inferred_allele.AlleleList()]
        aligned_imgt_alleles = self.allele_storage.GetAllelesByVGene(allele_basic.ModifyGeneName(v_gene))
        main_allele_id = aligned_imgt_alleles.MinimalAlleleId()
#        print '==== Ambiguous allele'
#        print "V gene: " + v_gene + ", allele: " + str(inferred_allele)
#        print allele_ids
        # 1. generate set of SNPs wrt to the first allele of the gene
        known_snps, ambiguous_positions = self._ComputeSNPList(allele_ids, aligned_imgt_alleles, main_allele_id)
#        print 'Known SNPs: ' + str(known_snps)
#        print 'Ambiguous positions: ' + str(ambiguous_positions)
        # 2. create sequence based on this set
        aligned_seq_list = list(aligned_imgt_alleles.GetSequenceByAllele(main_allele_id))
        # 3. add modifications
        main_imgt_allele = self.imgt_numbered_alleles.GetIMGTAllele(v_gene, main_allele_id)
        novel_descr, novel_snps = self._GetNovelModifications(main_imgt_allele, main_allele_id,
                                                                    aligned_imgt_alleles,
                                                                    inferred_allele.MainModification())
        snp_list = self._MergeSNPLists(known_snps, novel_snps)
        modified_sequence = self._CreateSequenceBySNPList(aligned_seq_list, snp_list, ambiguous_positions)
        # 4. computing allele type
        allele_type, allele_descr, allele_list = self._GetAlleleAndDescription(inferred_allele, v_gene,
                                                                               ambiguous_positions, allele_ids,
                                                                               snp_list, novel_descr)
        return annotated_alleles.AnnotatedVGeneAllele(allele_type, v_gene, allele_descr, snp_list, ambiguous_positions,
                                    modified_sequence, allele_list, novel_snps)

    def _FindClosestAllele(self, v_gene, imgt_allele_ids, snp_set):
        position_dict = dict()
        for p in snp_set:
            position_dict[p[0]] = p[1]
        aligned_imgt_alleles = self.allele_storage.GetAllelesByVGene(allele_basic.ModifyGeneName(v_gene))
        best_num_diff = len(snp_set)
        best_allele = -1
        for allele_id in imgt_allele_ids:
            allele_num_diff = 0
            for p in aligned_imgt_alleles.PolymoprhismIter():
                if p not in position_dict:
                    continue
                allele_nucl = aligned_imgt_alleles.GetNucleotideByAlleleAndPos(allele_id, p)
                if allele_nucl == position_dict[p]:
                    allele_num_diff += 1
            if allele_num_diff <= best_num_diff:
                best_num_diff = allele_num_diff
                best_allele = allele_id
        return best_allele