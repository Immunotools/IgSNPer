import allele_basic

class AnnotatedVGeneAllele:
    def __init__(self, allele_type, v_gene, description, snp_set, ambiguous_positions, aligned_sequence,
                 allele_list, modification_list):
        self.allele_type = allele_type
        self.v_gene = v_gene
        self.description = description
        self.snp_set = snp_set
        self.ambiguous_positions = ambiguous_positions
        self.sequence = ''
        self.aligned_sequence = aligned_sequence
        self.allele_list = allele_list
        self.modification_list = modification_list

    def __eq__(self, other):
        return self.Sequence() == other.Sequence()

    def __hash__(self):
        return hash(self.Sequence())

    def Sequence(self):
        if self.sequence == '':
            for c in self.aligned_sequence:
                if c != '-':
                    self.sequence += c
        return self.sequence

    def AlignedSequence(self):
        return self.aligned_sequence

    def Type(self):
        return self.allele_type

    def SNPIter(self):
        for snp in self.snp_set:
            yield snp

    def SNPs(self):
        return self.snp_set

    def UpdateSNPSet(self, new_snp_set):
        self.snp_set = new_snp_set

    def AmbiguousPositions(self):
        return self.ambiguous_positions

    def __str__(self):
        return self.description

    def __repr__(self):
        return self.description

    def Germline(self):
        return self.allele_type == allele_basic.AlleleType.GERMLINE

    def Novel(self):
        return self.allele_type == allele_basic.AlleleType.NOVEL

    def AmbiguousGermline(self):
        return self.allele_type == allele_basic.AlleleType.AMBIGUOUS_GERMLINE

    def AmbiguousNovel(self):
        return self.allele_type == allele_basic.AlleleType.AMBIGUOUS_NOVEL

    def Ambiguous(self):
        return self.AmbiguousGermline() or self.AmbiguousNovel()

    def MainAllele(self):
        return self.allele_list[0]

    def AlleleList(self):
        return self.allele_list

    def MainModification(self):
        return self.modification_list

    def ToLongString(self):
        return 'Type: ' + str(self.allele_type) + ', description: ' + self.description + '\nSequence: ' + self.Sequence()

    def ToShortString(self):
        # todo refactor
        splits = self.description.split('_')
#        print splits
        to_str = splits[0]
        if len(splits) > 1:
            to_str += '*'
        return to_str

    def V(self):
        return self.v_gene

    def GetNucleotideByAlignmentPos(self, align_pos):
        return self.aligned_sequence[align_pos]


class Haplotype:
    def __init__(self, v_gene, individual):
        self.v_gene = v_gene
        self.individual = individual
        self.alleles = []

    def AddAllele(self, allele):
        """
        @type allele : AnnotatedVGeneAllele
        """
        self.alleles.append(allele)

    def __getitem__(self, ind):
        return self.alleles[ind]

    def __iter__(self):
        for a in self.alleles:
            yield a

    def __len__(self):
        return len(self.alleles)

    def __hash__(self):
        hash_value = 1
        for a in self.alleles:
            hash_value *= hash(a)
        return hash_value

    def __eq__(self, other):
        if len(self) != len(other):
            return False
        seq_set = set([a.Sequence() for a in self.alleles])
        for a in other:
            if a.Sequence() not in seq_set:
                return False
        return True

    def __str__(self):
        h_str = '-'.join([str(a) for a in self.alleles])
        return h_str

    def GetNucleotidesByPosition(self, aligned_pos):
        nucls = []
        for a in self.alleles:
            nucls.append(a.AlignedSequence()[aligned_pos])
        return nucls