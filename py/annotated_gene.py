import allele_basic
import aligned_alleles
import annotated_alleles
import annotated_allele_unifier

from collections import Counter

class AnnotatedGene:
    def __init__(self, v_gene, allele_storage, allele_unifier):
        """
        @type allele_storage : aligned_alleles.AlignedAlleleStorage
        @type allele_unifier : annotated_allele_unifier.AlleleUnifier
        """
        self.v_gene = v_gene
        self.allele_storage = allele_storage
        self.allele_unifier = allele_unifier
        # inner fields
        self.allele_list = [] # all alleles
        self.individuals = [] # individual list
        self.individual_pointers = [] # index of individuals in allele_list
        self.alleles_unified = False
        #
        self.haplotype_dict = dict() # haplotype -> list of individuals
        self.allele_dict = dict() # allele -> list of individuals
        self.ind_haplotype_dict = dict() # ind -> haplotype

    def AddIndividualAlleles(self, ind_id, ind_alleles):
        """
        @type haplotype : annotated_alleles.Haplotype
        """
        self.alleles_unified = False
        self.individuals.append(ind_id)
        self.individual_pointers.append(len(self.allele_list))
        for a in ind_alleles:
            self.allele_list.append(a)

    def _UnifyAlleles(self):
        if self.alleles_unified:
            return
        self.alleles_unified = True
        unified_alleles = self.allele_unifier.UnifySNPs(self.v_gene, self.allele_list)
        self.allele_list = unified_alleles
        if self.individual_pointers[-1] != len(self.allele_list):
            self.individual_pointers.append(len(self.allele_list))
        for i in range(len(self.individual_pointers) - 1):
            start_pos = self.individual_pointers[i]
            end_pos = self.individual_pointers[i + 1]
            ind_id = self.individuals[i]
            haplotype = annotated_alleles.Haplotype(self.v_gene, ind_id)
            for j in range(start_pos, end_pos):
                haplotype.AddAllele(unified_alleles[j])
            if haplotype not in self.haplotype_dict:
                self.haplotype_dict[haplotype] = []
            self.haplotype_dict[haplotype].append(ind_id)
            for a in haplotype:
                if a not in self.allele_dict:
                    self.allele_dict[a] = []
                self.allele_dict[a].append(a)
            self.ind_haplotype_dict[ind_id] = haplotype

    def V(self):
        return self.v_gene

    def AlleleIter(self):
        self._UnifyAlleles()
        for a in self.allele_dict:
            yield a

    def NumAlleles(self):
        return len(self.allele_dict)

    def GetIndividualsByAllele(self, allele):
        self._UnifyAlleles()
        return self.allele_dict[allele]

    def NumIndividualsByAllele(self, allele):
        self._UnifyAlleles()
        return len(self.allele_dict[allele])

    def HaplotypeIter(self):
        self._UnifyAlleles()
        for h in self.haplotype_dict:
            yield h

    def GetIndividualsByHaplotype(self, haplotype):
        self._UnifyAlleles()
        return self.haplotype_dict[haplotype]

    def NumIndividualsByHaplotype(self, haplotype):
        self._UnifyAlleles()
        return len(self.haplotype_dict[haplotype])

    def SNPs(self):
        self._UnifyAlleles()
        return self.allele_list[0].SNPs()

    def PolymorphismPositionIter(self):
        self._UnifyAlleles()
        for snp in self.allele_list[0].SNPs():
            yield snp[0]

    def ComputePositionState(self, snp_pos):
        snp_counter = []
        for h in self.HaplotypeIter():
            pos_nucls = allele_basic.NuclsToStr(h.GetNucleotidesByPosition(snp_pos))
            snp_counter.extend([pos_nucls] * self.NumIndividualsByHaplotype(h))
        return Counter(snp_counter)

    def ComputeIndividualsByPositionState(self, snp_pos):
        pos_states = dict() # nucl -> individuals
        for h in self.HaplotypeIter():
            individuals = self.GetIndividualsByHaplotype(h)
            nucls = allele_basic.NuclsToStr(h.GetNucleotidesByPosition(snp_pos))
            if nucls not in pos_states:
                pos_states[nucls] = []
            pos_states[nucls].extend(individuals)
        return pos_states

    def ContainsIndividual(self, ind_id):
        return ind_id in self.ind_haplotype_dict

    def GetHaplotypeByIndividual(self, ind_id):
        return self.ind_haplotype_dict[ind_id]

    def GetSNPsByIndividual(self, ind_id, snp_pos):
        haplotype = self.ind_haplotype_dict[ind_id]
        nucls = allele_basic.NuclsToStr(haplotype.GetNucleotidesByPosition(snp_pos))
        return nucls