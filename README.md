# IgSNPer

a tool for transforming immunoglobulin V gene (IGHV) alleles into Single Nucleotide Polymorphisms (SNPs).

## Usage: 
```
python ig_snper.py -c config.txt -o output_dir
```

## Config format:
```
TiggerOutputDir       ProjectID
Genotype_dir1         Project_id1
Genotype_dir2         Project_id2
Genotype_dir3         Project_id3
```
where Genotype_dir1-3 are directories containing individual IGHV alleles inferred by TigGER tool (Gadala-Maria et al., 2019). Project identificators Project_id1-3 will be used in IgSNPer output files. An example of TiggerOutputDir can be found in data/vdjbase/tigger_results_p1/Genotypes, it corresponds to the project P1 (PRJEB26509; Gidoni et al., 2019) at the VDJbase (https://www.vdjbase.org/). The config file for the P1 project can be found in test_datasets/config_p1.txt and looks like this:
```
TiggerOutputDir	                          ProjectID
data/vdjbase/tigger_results_p1/Genotypes	p1
```

## Requirements:
- python 2 or 3
- Biopython
- Seaborn
- NumPy
- pandas
