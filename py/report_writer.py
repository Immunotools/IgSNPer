import logging

import allele_basic
import annotated_gene

import os
import sys

class HTMLWriter:
    def __init__(self, output_fname):
        self.fhandler = open(output_fname, "w")

    def _WriteHeader(self, level, text, align):
        self.fhandler.write("<h" + str(level) + " align = " + align + ">" + text + "</h" + str(level) + ">\n")

    def _GetAlignment(self, align_to_center = True):
        if not align_to_center:
            return "left"
        return "center"

    def _WriteTableCell(self, elem):
        self.fhandler.write("<td>" + str(elem) + "</td>")

    def _WriteTableRow(self, row):
        self.fhandler.write("<tr>\n")
        for r in row:
            self._WriteTableCell(r)
        self.fhandler.write("</tr>\n")

    def WriteH1(self, text, align_to_center = True):
        self._WriteHeader(1, text, self._GetAlignment(align_to_center))

    def WriteH2(self, text, align_to_center = True):
        self._WriteHeader(2, text, self._GetAlignment(align_to_center))

    def WriteH3(self, text, align_to_center = True):
        self._WriteHeader(3, text, self._GetAlignment(align_to_center))

    # width in percent
    def WriteImage(self, path_to_image, width = 60):
        self.fhandler.write("<p align = center>\n")
        self.fhandler.write("<image src = " + path_to_image + " width = " + str(width) + "%></image>\n")
        self.fhandler.write("</p>\n")

    def WriteTable(self, col_names, row_names, values, width = 60):
        if len(values) == 0:
            return
        if len(row_names) != len(values):
            logging.warning("# rows in table and # row names are not consistent")
            sys.exit(1)
        if len(col_names) != len(values[0]):
            logging.warning("# columns in table and # column names are not consistent")
            sys.exit(1)
        self.fhandler.write("<table width = " + str(width) + "% align = center>\n")
        header_row = [""]
        header_row.extend(col_names)
        self._WriteTableRow(header_row)
        for i in range(0, len(row_names)):
            cur_row = [row_names[i]]
            cur_row.extend(values[i])
            self._WriteTableRow(cur_row)
        self.fhandler.write("</table>\n")

    def WriteEmptyLine(self):
        self.fhandler.write("<br>\n")

    def WriteHorizontalLine(self, width = 100):
        self.fhandler.write("<hr width=" + str(width) + "%>\n")

    def WriteParagraph(self, text, alignment = "justify"):
        self.fhandler.write("<p align = \"" + alignment + "\">\n")
        self.fhandler.write(text + '\n')
        self.fhandler.write('</p>\n')

    def CloseFile(self):
        self.fhandler.close()

    def WriteImageWithTitle(self, fname, title, width = 60):
        if title != '':
            self.WriteH2(title)
        self.WriteImage(fname, width=width)


class ReportCompiler:
    def __init__(self, dir_config, output_fname, gene):
        """
        @type gene : annotated_gene.AnnotatedGene
        """
        self.image_format = '.svg'
        self.dir_config = dir_config
        self.output_fname = output_fname
        self.gene = gene
        self.v = gene.V()
        self.mod_v = allele_basic.ModifyGeneName(self.v)

    def _GetAlleleImage(self):
        return os.path.join(self._CreateSubPath(self.dir_config['alleles']), self.mod_v + self.image_format)

    def _GetGenotypeImage(self):
        return os.path.join(self._CreateSubPath(self.dir_config['genotypes']), self.mod_v + self.image_format)

    def _GetHaplotypeImage(self):
        return os.path.join(self._CreateSubPath(self.dir_config['haplotypes']), self.mod_v + self.image_format)

    def _GetSNPImages(self):
        snp_images = os.listdir(self.dir_config['snps'])
        return [os.path.join(self._CreateSubPath(self.dir_config['snps']), f) for f in snp_images if
                f.find(self.mod_v) != -1 and f.find(self.image_format) != -1]

    def _CreateSubPath(self, path):
        return os.path.join('..', os.path.basename(path))

    def Output(self):
        html_writer = HTMLWriter(self.output_fname)
        html_writer.WriteH1(self.v)
        # haplotype
        h_img = self._GetHaplotypeImage()
        html_writer.WriteImageWithTitle(h_img, 'Haplotypes', 50)
        # alleles
        allele_img = self._GetAlleleImage()
        html_writer.WriteImageWithTitle(allele_img, 'Alleles')
        # genotypes
        g_img = self._GetGenotypeImage()
        html_writer.WriteImageWithTitle(g_img, 'Multi-SNP genotypes')
        # snps
        snp_images = self._GetSNPImages()
        html_writer.WriteH2('SNPs')
        for snp_img in sorted(snp_images, key = lambda x : int(x.split('_')[-1].split('.')[0])):
            html_writer.WriteImage(snp_img, 40)
        # close
        html_writer.CloseFile()