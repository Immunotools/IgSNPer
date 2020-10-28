from enum import Enum

class AlleleType(Enum):
    GERMLINE = 0
    NOVEL = 1
    AMBIGUOUS_GERMLINE = 2
    AMBIGUOUS_NOVEL = 3


class AlleleModification:
    def __init__(self, mod_str):
        self.src_nucl = mod_str[0]
        self.dst_nucl = mod_str[-1]
        self.pos = int(mod_str[1 : -1]) - 1

    def Pos(self):
        return self.pos

    def DstNucl(self):
        return self.dst_nucl

    def SrcNucl(self):
        return self.src_nucl

    def __str__(self):
        return str(self.Pos()) + ':' + self.DstNucl()

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        return self.Pos() == other.Pos() and self.DstNucl() == self.DstNucl()

    def __hash__(self):
        return hash(self.Pos()) * hash(self.DstNucl())


class InferredAllele:
    def __init__(self):
        self.allele_ids = []
        self.modifications = []

    def ParseFromString(self, allele_str):
        splits = allele_str.split('_')
        for s in splits:
            if s.isdigit():
                self.allele_ids.append(int(s))
            else:
                self.modifications.append(AlleleModification(s))

    def AddAlleleIdWithModifications(self, allele_id, mod_list):
        self.allele_ids.append(allele_id)
        self.modifications.append(mod_list)

    def MainAllele(self):
        return self.allele_ids[0]

    def MainModification(self):
        return self.modifications

    def NumModifications(self):
        return len(self.modifications)

    def Unambiguous(self):
        return len(self.allele_ids) == 1

    def Ambiguous(self):
        return not self.Unambiguous()

    def AmbiguousGermline(self):
        if self.Unambiguous():
            return False
        return len(self.modifications) == 0

    def AmbiguousNovel(self):
        if self.Unambiguous():
            return False
        return len(self.modifications) != 0

    def Germline(self):
        return self.Unambiguous() and len(self.modifications) == 0

    def Novel(self):
        return self.Unambiguous() and len(self.modifications) != 0

    def Type(self):
        if self.Germline():
            return AlleleType.GERMLINE
        elif self.Novel():
            return AlleleType.NOVEL
        elif self.AmbiguousGermline():
            return AlleleType.AMBIGUOUS_GERMLINE
        return AlleleType.AMBIGUOUS_NOVEL

    def __repr__(self):
        return str(self.allele_ids) + ' - ' + str(self.modifications)

    def AlleleList(self):
        return self.allele_ids

    def __hash__(self):
        return hash(self.MainAllele()) * hash(','.join([str(m) for m in self.MainModification()]))

    def __eq__(self, other):
        if self.Ambiguous() or other.Ambiguous():
            return False
        if self.MainAllele() != other.MainAllele():
            return False
        if self.NumModifications() != other.NumModifications():
            return False
        for m in other.MainModification():
            if m not in self.MainModification():
                return False
        return True


class CountedGeneAlleles:
    def __init__(self):
        self.alleles = []
        self.counts = []

    def Add(self, allele, count):
        self.alleles.append(allele)
        self.counts.append(count)

    def __iter__(self):
        for a, c in zip(self.alleles, self.counts):
            yield a, c

    def __len__(self):
        return len(self.alleles)

    def Empty(self):
        return len(self) == 0


class AlleleVisConfig:
    def __init__(self):
        self.min_value = 0
        self.max_value = 20
        self.cmap = 'tab20c'

    def MinValue(self):
        return self.min_value

    def MaxValue(self):
        return self.max_value

    def ColorMap(self):
        return self.cmap

    def GetEmptyAlleleColor(self):
        return self.MaxValue()

    def GetColorLabelByAllele(self, allele):
        if allele.Germline():
            return allele.MainAllele(), str(allele.MainAllele())
        elif allele.Novel():
            return allele.MainAllele(), str(allele.MainAllele()) + ' + ' + ','.join([str(m) for m in allele.MainModification()])
        allele_str = ', '.join([str(i) for i in allele.AlleleList()])
        if allele.AmbiguousNovel():
            allele_str += ' + ' + ','.join([str(m) for m in allele.MainModification()])
        return 0, allele_str


def ModifyGeneName(v_gene):
    return v_gene.replace('/', '_')


def NuclsToStr(pos_nucls):
    nucl_set = set(pos_nucls)
    if len(nucl_set) == 1:
        return pos_nucls[0]
    return '/'.join(sorted(pos_nucls))
