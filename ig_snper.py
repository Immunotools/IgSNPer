import os
import sys
import logging
import getopt
import pandas as pd

sys.path.append('py/')
import utils
import aligned_alleles
import allele_basic
import raw_tigger_results
import annotated_tigger_result
import annotated_tigger_results
import annotated_allele_merger
import annotated_allele_unifier
import annotated_gene_vis
import report_writer
import novel_allele_visualizer
import snp_result_writer

def ProcessRawTiggerResults(config_df, allele_storage, imgt_numbered_alleles, output_dir, logger):
    annotated_results = dict()
    for i in range(len(config_df)):
        project_id = config_df['ProjectID'][i]
        project_dir = os.path.join(output_dir, project_id + '_processed')
        os.mkdir(project_dir)
        tigger_dir = config_df['TiggerOutputDir'][i]
        tigger_files = os.listdir(tigger_dir)
        for f in tigger_files:
            if f == '.DS_Store':
                continue
            full_path = os.path.join(tigger_dir, f)
            basename = os.path.basename(f)
            ind_id = basename[ : len(basename) - len('_geno_H_binom.tab')]
            tigger_df = pd.read_csv(full_path, sep = '\t')
            raw_result = raw_tigger_results.RawTiggerResult(tigger_df, logger)
            annotated_result = annotated_tigger_result.AnnotatedTiggerResult(raw_result, allele_storage,
                                                                             imgt_numbered_alleles, logger)
            output_fname = os.path.join(project_dir, basename + '.txt')
            annotated_result.WriteToTxt(output_fname)
            annotated_results[(project_id, ind_id)] = annotated_result
    return annotated_results


def OutputIndividualSNPs(subject_list, annotated_results, snp_writer, output_dir):
    for project_id, ind_id in subject_list:
        ind_id_str = project_id + '_' + ind_id
        snp_writer.OutputSNPsPerSubject(annotated_results, (project_id, ind_id),
                                        os.path.join(output_dir, ind_id_str + '.txt'))


def OutputResultsPerGene(annotated_results, snp_writer, gene_vis, config, logger):
    """
    @type config : IgSNPerConfig
    """
    for annot_gene in sorted(annotated_results.GeneIterator(), key = lambda x : x.V()):
        logger.info(annot_gene.V() + '...')
#        for h in sorted(annot_gene.HaplotypeIter(), key = lambda x : annot_gene.NumIndividualsByHaplotype(x),
#                        reverse = True):
#            print h, annot_gene.NumIndividualsByHaplotype(h)
        modified_name = allele_basic.ModifyGeneName(annot_gene.V())
        gene_vis.VisualizeAlleles(annot_gene, os.path.join(config.AlleleDir(), modified_name))
        gene_vis.VisualizeMultiSNPGenotypes(annot_gene, os.path.join(config.GenotypeDir(), modified_name))
        gene_vis.VisualizeHaplotypes(annot_gene, os.path.join(config.HaplotypeDir(), modified_name))
        gene_vis.VisualizeSNPStates(annot_gene, config.SNPPlotDir())
        gene_vis.OutputHaplotypesByGene(annot_gene, os.path.join(config.HaplotypeTxtDir(), modified_name + '.txt'))
        report_wr = report_writer.ReportCompiler(config, os.path.join(config.HTMLDir(), modified_name + '.html'),
                                                 annot_gene)
        snp_writer.OutputSNPStates(annot_gene, os.path.join(config.SNPTxtDir(), modified_name + '_snp.txt'))
        report_wr.Output()

class IgSNPerConfig:
    def __init__(self, args):
        self.config_fname = ''
        self.output_dir = ''
        self.aligned_allele_dir = os.path.abspath('data/imgt_alleles/IGHV_allele_alignments')
        self.imgt_numbering_fasta = os.path.abspath('data/imgt_alleles/IGHV_imgt_numbering.fa')
        self._ParseArgs(args)
        self._PrepareOutputDirs()

    def _ParseArgs(self, args):
        try:
            options, remainder = getopt.getopt(args[1:], 'c:o:', ['help'])
        except getopt.GetoptError as err:
            print str(err)  # will print something like "option -a not recognized"
            sys.exit(2)
        for opt, arg in options:
            if opt == "-c":
                self.config_fname = os.path.abspath(arg)
            elif opt == '-o':
                self.output_dir = os.path.abspath(arg)
            elif opt == '--help':
                self._PrintHelp()
                sys.exit(0)
            else:
                assert False, "unhandled option"

    def _PrintHelp(self):
        print('python ig_snper.py -c config.txt -o output_dir')
        sys.exit(0)

    def _CreateDirs(self):
        dirs = ['alleles', 'genotypes', 'haplotypes', 'snp_plots', 'html', 'haplotypes_txt', 'snp_txt', 'subject_snp_txt']
        dir_config = dict()
        for d in dirs:
            dir_config[d] = os.path.join(self.output_dir, d)
            os.mkdir(dir_config[d])
        return dir_config

    def _PrepareOutputDirs(self):
        utils.PrepareOutputDir(self.output_dir)
        self.dir_config = self._CreateDirs()

    def ConfigureLogger(self):
        log_formatter = logging.Formatter('%(asctime)s\t%(levelname)s\t%(message)s')
        self.logger = logging.getLogger()
        self.logger.setLevel(logging.INFO)
        fileHandler = logging.FileHandler(os.path.join(self.output_dir, 'log.info'))
        fileHandler.setFormatter(log_formatter)
        self.logger.addHandler(fileHandler)
        consoleHandler = logging.StreamHandler()
        consoleHandler.setFormatter(log_formatter)
        self.logger.addHandler(consoleHandler)
        return self.logger

    def SubjectSNPDir(self):
        return self.dir_config['subject_snp_txt']

    def AlleleDir(self):
        return self.dir_config['alleles']

    def GenotypeDir(self):
        return self.dir_config['genotypes']

    def HaplotypeDir(self):
        return self.dir_config['haplotypes']

    def HaplotypeTxtDir(self):
        return self.dir_config['haplotypes_txt']

    def SNPPlotDir(self):
        return self.dir_config['snp_plots']

    def SNPTxtDir(self):
        return self.dir_config['snp_txt']

    def HTMLDir(self):
        return self.dir_config['html']

def main(args):
    config = IgSNPerConfig(args)
    config_df = pd.read_csv(config.config_fname, delim_whitespace = True)
    logger = config.ConfigureLogger()
    allele_storage = aligned_alleles.ReadAlignedAlleles(config.aligned_allele_dir)
    imgt_numbered_alleles = aligned_alleles.IMGTNumberedAlleles(config.imgt_numbering_fasta)
    annotated_dict = ProcessRawTiggerResults(config_df, allele_storage, imgt_numbered_alleles, config.output_dir, logger)
    allele_merger = annotated_allele_merger.AgressiveAlleleMerger(allele_storage, imgt_numbered_alleles, logger)
    allele_unifier = annotated_allele_unifier.AlleleUnifier(allele_storage, logger)
    annotated_results = annotated_tigger_results.AnnotatedTiggerResults(allele_storage, imgt_numbered_alleles,
                                                                        allele_merger, allele_unifier, logger)
    for project_id, ind_id in annotated_dict:
        annotated_results.AddAllele(project_id, ind_id, annotated_dict[(project_id, ind_id)])
    annotated_results.PrintGluingStats()
    logger.info("Writing haplotypes info...")
    novel_allele_vis = novel_allele_visualizer.NovelAlleleVisualizer(annotated_results, allele_storage)
#    novel_allele_vis.OutputAllelePositions(os.path.join(config.output_dir, '_novel_pos.pdf'))
#    novel_allele_vis.OutputPutativePositions(os.path.join(config.output_dir, '_putative_novel_pos.pdf'))
    gene_vis = annotated_gene_vis.AnnotatedGeneVisualizer(logger)
    snp_writer = snp_result_writer.SNPResultWriter()
    OutputIndividualSNPs(annotated_dict.keys(), annotated_results, snp_writer, config.SubjectSNPDir())
    OutputResultsPerGene(annotated_results, snp_writer, gene_vis, config, logger)


if __name__ == "__main__":
    main(sys.argv)

