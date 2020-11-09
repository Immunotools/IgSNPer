# IgSNPer

a tool for transforming immunoglobulin V gene (IGHV) alleles into Single Nucleotide Polymorphisms (SNPs).

## Requirements
- python 2 or 3
- Biopython
- Matplotlib
- NumPy
- pandas
- seaborn

## Usage
```
python ig_snper.py -c config.txt -o output_dir
```

Example:
```
python ig_snper.py -c test_datasets/config_p1.txt -o ig_snper_p1_processed
```

## Input files
### Config format
```
TiggerOutputDir       ProjectID
Genotype_dir1         Project_id1
Genotype_dir2         Project_id2
Genotype_dir3         Project_id3
```
`Genotype_dirN` is a directory containing individual IGHV alleles inferred by TigGER tool (Gadala-Maria et al., 2019). Project identifier `Project_idN` is used in IgSNPer output files. An example of TiggerOutputDir can be found in `data/vdjbase/tigger_results_p1/Genotypes`, it corresponds to the project P1 (PRJEB26509; Gidoni et al., 2019) at the VDJbase (https://www.vdjbase.org/). The config file for the P1 project can be found in `test_datasets/config_p1.txt` and can take either of the following formats:
```
TiggerOutputDir	                          ProjectID
data/vdjbase/tigger_results_p1/Genotypes  p1
```
or
```
TiggerFilePath                           ProjectID     SubjectID
samples/P1/P1_I100_S1/P1_I100_S1.tsv     P1            P1_I100
```

### TigGER genotype files
IgSNPer scans the directories described in the config and seeks for genotype files reported by TIgGER with names in the following format:
```
individualId_geno_H_binom.tab
```
E.g., for file `P1_I1_S1_geno_H_binom.tab`, individual ID is `P1_I1_S1`. 
Individual identifiers are extracted by IgSNPer and used in output files.


## Output files
IgSNPer creates an output directory `output_dir` with the following structure:
```
output_dir
|_ html_reports
|_ plot_alleles
|_ plot_genotypes
|_ plot_haplotypes
|_ plot_snps
|_ txt_haplotypes_by_gene
|_ txt_snps_by_gene
|_ txt_snps_by_subject
```

* Directories starting with `plot_` contain plots showing individual SNPs, SNPs combining alleles, haplotypes, and genotypes of IGHV genes. 
* The directory `html_reports` contains HTML files, each of which corresponds to a single IGHV gene and combines plots corresponding to it.
* The directory `txt_haplotypes_by_gene` contains TXT files corresponding to IGHV genes and describing haplotypes across all individuals. 
* The directory `txt_snps_by_gene` contains TXT files corresponding to IGHV genes and describing SNPs across all individuals.
* The directory `txt_snps_by_subject` contains TXT files corresponding to individuals and desribing their SNPs across all IGHV genes. 

### Individual SNP file
A fragment of the SNP file for individual `P1_I1_S1` is shown below:
```
Dataset	Individual	Gene	    Haplotype   Pos	    State   AllelePos
p1	P1_I1_S1	IGHV1-18    1	        3	    G	    1
p1	P1_I1_S1	IGHV1-18    1	        5	    T	    1
p1	P1_I1_S1	IGHV1-18    1	        11	    G	    1
p1	P1_I1_S1	IGHV1-18    1	        23	    A	    1
p1	P1_I1_S1	IGHV1-18    1	        24	    G	    1
p1	P1_I1_S1	IGHV1-18    1	        95	    T	    1
p1	P1_I1_S1	IGHV1-18    1	        293	    A	    1
...
p1	P1_I1_S1	IGHV3-53    1-4	        15	    G	    1
p1	P1_I1_S1	IGHV3-53    1-4	        18	    T	    1
p1	P1_I1_S1	IGHV3-53    1-4	        33	    A/G	    1
p1	P1_I1_S1	IGHV3-53    1-4	        213	    C/G	    1
p1	P1_I1_S1	IGHV3-53    1-4	        260	    C/T	    1
...
```
Each line corresponds to a single V gene SNP described by a position (the column `Pos`) and a state (the column `State`). Positions are 0-based and given with respect to the IMGT allele with the minimum identifier. The identifier of an allele is the number coming after `*` in the allele name, e.g., identifier of allele `IGHV1-69*10` is `10`. The identifier of the allele determining the SNP position is provided in the column `AllelePos`. E.g., positions of SNPs of `IGHV1-18` are provided with respect to allele `IGHV1-18*01`.
