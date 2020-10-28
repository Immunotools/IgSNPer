import os
import sys
import logging
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


def CreateDirs(output_dir):
    dirs = ['alleles', 'genotypes', 'haplotypes', 'snps', 'html', 'haplotypes_txt', 'snp_txt', 'subject_snp_txt']
    dir_config = dict()
    for d in dirs:
        dir_config[d] = os.path.join(output_dir, d)
        os.mkdir(dir_config[d])
    return dir_config


def OutputIndividualSNPs(subject_list, annotated_results, snp_writer, output_dir):
    for project_id, ind_id in subject_list:
        ind_id_str = project_id + '_' + ind_id
        snp_writer.OutputSNPsPerSubject(annotated_results, (project_id, ind_id),
                                        os.path.join(output_dir, ind_id_str + '.txt'))


def OutputResultsPerGene(annotated_results, snp_writer, gene_vis, dir_config, logger):
    for annot_gene in sorted(annotated_results.GeneIterator(), key = lambda x : x.V()):
        logger.info('== ' + annot_gene.V())
#        for h in sorted(annot_gene.HaplotypeIter(), key = lambda x : annot_gene.NumIndividualsByHaplotype(x),
#                        reverse = True):
#            print h, annot_gene.NumIndividualsByHaplotype(h)
        modified_name = allele_basic.ModifyGeneName(annot_gene.V())
        gene_vis.VisualizeAlleles(annot_gene, os.path.join(dir_config['alleles'], modified_name))
        gene_vis.VisualizeMultiSNPGenotypes(annot_gene, os.path.join(dir_config['genotypes'], modified_name))
        gene_vis.VisualizeHaplotypes(annot_gene, os.path.join(dir_config['haplotypes'], modified_name))
        gene_vis.VisualizeSNPStates(annot_gene, dir_config['snps'])
        gene_vis.OutputHaplotypesByGene(annot_gene, os.path.join(dir_config['haplotypes_txt'], modified_name + '.txt'))
        report_wr = report_writer.ReportCompiler(dir_config, os.path.join(dir_config['html'], modified_name + '.html'),
                                                 annot_gene)
        snp_writer.OutputSNPStates(annot_gene, os.path.join(dir_config['snp_txt'], modified_name + '_snp.txt'))
        report_wr.Output()


def ConfigureLogger(output_dir):
    log_formatter = logging.Formatter('%(asctime)s\t%(levelname)s\t%(message)s')
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    fileHandler = logging.FileHandler(os.path.join(output_dir, 'log.info'))
#    fileHandler.setFormatter(log_formatter)
    logger.addHandler(fileHandler)

    consoleHandler = logging.StreamHandler()
#    consoleHandler.setFormatter(log_formatter)
    logger.addHandler(consoleHandler)
    return logger

def main(config_fname, allele_dir, output_dir):
    config_df = pd.read_csv(config_fname, delim_whitespace = True)
    utils.PrepareOutputDir(output_dir)
    dir_config = CreateDirs(output_dir)
    logger = ConfigureLogger(output_dir)
    allele_storage = aligned_alleles.ReadAlignedAlleles(allele_dir)
    imgt_numbered_alleles = aligned_alleles.IMGTNumberedAlleles('IGHV_imgt_numbering.fa')
    annotated_dict = ProcessRawTiggerResults(config_df, allele_storage, imgt_numbered_alleles, output_dir, logger)
    allele_merger = annotated_allele_merger.AgressiveAlleleMerger(allele_storage, imgt_numbered_alleles, logger)
    allele_unifier = annotated_allele_unifier.AlleleUnifier(allele_storage, logger)
    annotated_results = annotated_tigger_results.AnnotatedTiggerResults(allele_storage, imgt_numbered_alleles,
                                                                        allele_merger, allele_unifier, logger)
    for project_id, ind_id in annotated_dict:
        annotated_results.AddAllele(project_id, ind_id, annotated_dict[(project_id, ind_id)])
    annotated_results.PrintGluingStats()
    logger.info("Writing haplotypes info...")
    novel_allele_vis = novel_allele_visualizer.NovelAlleleVisualizer(annotated_results, allele_storage)
    novel_allele_vis.OutputAllelePositions(os.path.join(output_dir, '_novel_pos.pdf'))
    novel_allele_vis.OutputPutativePositions(os.path.join(output_dir, '_putative_novel_pos.pdf'))
    gene_vis = annotated_gene_vis.AnnotatedGeneVisualizer(logger)
    snp_writer = snp_result_writer.SNPResultWriter()
    OutputIndividualSNPs(annotated_dict.keys(), annotated_results, snp_writer, dir_config['subject_snp_txt'])
    OutputResultsPerGene(annotated_results, snp_writer, gene_vis, dir_config, logger)


if __name__ == "__main__":
    # python process_tigger_results.py config.txt allele_dir output_dir
    main(sys.argv[1], sys.argv[2], sys.argv[3])

