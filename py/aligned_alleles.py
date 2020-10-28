import os
import sys
import logging

from Bio import SeqIO

from collections import Counter

class AlignedAlleles:
    def __init__(self, v_gene, align_fasta):
        self.v_gene = v_gene
        self.align_fasta = align_fasta
        self.aligned_length = 0
        self._ReadAlignedAlleles()
        self._ComputeAlignmentBounds()
        self._FindPolymorphicPositions()
        self._MapPolymorphicPositions()

    def _ReadAlignedAlleles(self):
#        print 'Reading ' + self.align_fasta + '...'
        self.aligned_dict = dict()
        for r in SeqIO.parse(self.align_fasta, 'fasta'):
#            print r.id
            allele_id = int(r.id.split('*')[1])
            self.aligned_dict[allele_id] = str(r.seq).upper()
            self.aligned_length = len(self.aligned_dict[allele_id])

    def _ComputeAlignmentBounds(self):
        self.alignment_bounds = dict() # allele id -> alignment start, alignment end (both inclusive)
        for allele_id in self.aligned_dict:
            allele_seq = self.aligned_dict[allele_id]
            bounds = [-1, -1]
            for i in range(len(allele_seq)):
                if allele_seq[i] != '-':
                    bounds[0] = i
                    break
            for i in range(len(allele_seq)):
                pos = len(allele_seq) - i - 1
                if allele_seq[pos] != '-':
                    bounds[1] = pos
                    break
            self.alignment_bounds[allele_id] = bounds

    def _FindPolymorphicPositions(self):
        self.polymorphic_positions = [] # positions in multiple alignment
        for i in range(self.aligned_length):
            nucleotides = []
            for allele_id in self.aligned_dict:
                if i >= self.alignment_bounds[allele_id][0] and i <= self.alignment_bounds[allele_id][1]:
                    nucleotides.append(self.aligned_dict[allele_id][i])
            nucl_dict = Counter(nucleotides)
            if len(nucl_dict) == 1:
                continue
#            if len(nucl_dict) == 2 and '-' in nucl_dict.keys():
#                continue
            self.polymorphic_positions.append(i)

    def _MapPolymorphicPositions(self):
        self.position_map = dict() # alignment position -> dict: allele id -> allele position
        for pos in self.polymorphic_positions:
            self.position_map[pos] = dict()
        for allele_id in self.aligned_dict:
            aligned_allele_seq = self.aligned_dict[allele_id]
            seq_pos = -1
            for i in range(len(aligned_allele_seq)):
                if aligned_allele_seq[i] != '-1':
                    seq_pos += 1
                if i in self.polymorphic_positions:
                    self.position_map[i][allele_id] = seq_pos

    def MinimalAlleleId(self):
        return min(self.aligned_dict.keys())

    def AlignmentLength(self):
        return self.aligned_length

    def __len__(self):
        return len(self.aligned_dict)

    def PolymoprhismIter(self):
        # iterator over polymorphic positions in the multiple alignment
        for p in self.polymorphic_positions:
            yield p

    def AlleleIter(self):
        for allele_id in self.aligned_dict:
            yield allele_id

    def _CheckAlleleFatal(self, allele_id):
        if not self.ContainsAllele(allele_id):
            logging.warning("ERROR: gene " + self.v_gene + ' does not have allele ' + str(allele_id))
            sys.exit(1)

    def GetAllelePosByAlignmentPos(self, allele_id, pos):
        self._CheckAlleleFatal(allele_id)
        aligned_allele_seq = self.aligned_dict[allele_id]
        seq_pos = -1
        for i in range(len(aligned_allele_seq)):
            if aligned_allele_seq[i] != '-':
                seq_pos += 1
            if i == pos:
                return seq_pos
        return -1

    def GetNucleotideByAlleleAndPos(self, allele_id, alignment_pos):
        return self.aligned_dict[allele_id][alignment_pos]
#        return self.aligned_dict[allele_id][self.GetAllelePosByAlignmentPos(allele_id, alignment_pos)]

    def GetNucleotideByAlleleAndAllelePos(self, allele_id, allele_pos):
        return self.aligned_dict[allele_id][self.GetAlignmentPosByAllelePos(allele_id, allele_pos)]

    def ContainsAllele(self, allele_id):
        return allele_id in self.aligned_dict

    def GetAlignmentPosByAllelePos(self, allele_id, pos):
        self._CheckAlleleFatal(allele_id)
        aligned_allele_seq = self.aligned_dict[allele_id]
        seq_pos = -1
        for i in range(len(aligned_allele_seq)):
            if aligned_allele_seq[i] != '-':
                seq_pos += 1
            if seq_pos == pos:
                return i
        return -1

    def GetPolymorhismsByAlleles(self, allele_list):
        polymorhisms = []
        for pos in self.PolymoprhismIter():
            nucls = set()
            for allele_id in allele_list:
                nucl = self.GetNucleotideByAlleleAndPos(allele_id, pos)
                nucls.add(nucl)
#            print pos, nucls
            if len(nucls) > 1:
                polymorhisms.append(pos)
        return polymorhisms

    def PositionIsKnownPolymorphism(self, alignment_pos):
        return alignment_pos in self.position_map

    def GetSequenceByAllele(self, allele_id):
        self._CheckAlleleFatal(allele_id)
        return self.aligned_dict[allele_id]

    def GetAlleleNucleotidesByAlignmentPositions(self, alignment_pos):
        nucls = set()
        for allele_id in self.AlleleIter():
            nucls.add(self.GetNucleotideByAlleleAndPos(allele_id, alignment_pos))
        return nucls

    def GetDifferencesBetweenAlleles(self, allele1, allele2):
        self._CheckAlleleFatal(allele1)
        self._CheckAlleleFatal(allele2)
        snps = []
        for p in self.PolymoprhismIter():
            nucl1 = self.GetNucleotideByAlleleAndPos(allele1, p)
            nucl2 = self.GetNucleotideByAlleleAndPos(allele2, p)
            if nucl1 != nucl2:
                snps.append((p, nucl1, nucl2))
        return snps


class AlignedAlleleStorage:
    def __init__(self):
        self.v_dict = dict() # V gene -> AlignedAlleles

    def Add(self, v_gene, v_fasta):
        self.v_dict[v_gene] = AlignedAlleles(v_gene, v_fasta)

    def GetAllelesByVGene(self, v):
        return self.v_dict[v]

    def __len__(self):
        return len(self.v_dict)


def ReadAlignedAlleles(allele_dir):
    allele_storage = AlignedAlleleStorage()
    files = os.listdir(allele_dir)
    gene_dict = dict()
    suffix = '_align.fasta'
    for f in files:
        if f.find(suffix) == -1:
            continue
        v_gene = f[ : len(f) - len(suffix)]
        full_path = os.path.join(allele_dir, f)
        gene_dict[v_gene] = (full_path, os.path.join(allele_dir, v_gene + '.fasta'))
    for v_gene in gene_dict:
        file_to_read = gene_dict[v_gene][0]
        if len(open(file_to_read).readlines()) == 0:
            file_to_read = gene_dict[v_gene][1]
        allele_storage.Add(v_gene, file_to_read)
    logging.info("Alleles for " + str(len(allele_storage)) + ' V genes were extracted from ' + allele_dir)
    return allele_storage


class IMGTNumberedAllele:
    def __init__(self, allele_name, allele_seq):
        self.allele_name = allele_name
        self.allele_seq = allele_seq.upper() # seq with dots

    def GetNucleotideByIMGTPosition(self, imgt_pos):
        return self.allele_seq[imgt_pos]

    def GetPositionByIMGTPosition(self, imgt_pos):
        pos = -1
        for i in range(len(self.allele_seq)):
            if self.allele_seq[i] != '.':
                pos += 1
            if i == imgt_pos:
                return pos
        return -1


class IMGTNumberedAlleles:
    def __init__(self, imgt_numbered_fasta):
        self.v_gene_dict = dict() # V gene -> allele -> IMGTNumberedAllele
        for r in SeqIO.parse(imgt_numbered_fasta, 'fasta'):
            v_gene = r.id.split('|')[1]
            gene_base = v_gene.split('*')[0]
            v_allele = int(v_gene.split('*')[1])
            if gene_base not in self.v_gene_dict:
                self.v_gene_dict[gene_base] = dict()
            self.v_gene_dict[gene_base][v_allele] = IMGTNumberedAllele(v_gene, str(r.seq))

    def GetIMGTAllele(self, gene, allele_id):
        return self.v_gene_dict[gene][allele_id]