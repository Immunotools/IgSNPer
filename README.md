# IgSNPer

a tool for transformation of immunoglobulin V gene (IGHV) alleles into Single Nucleotide Polymorphisms (SNPs).

## Usage: 
```
python ig_snper.py -c config.txt -o output_dir
```

## Config format:
```
TiggerOutputDir       ProjectID
genotype_dir1         project_id1
genotype_dir2         project_id2
genotype_dir3         project_id3
```
where genotype_dir1-3 are directories containing individual IGHV alleles inferred by TigGER tool (Gadala-Maria et al., 2019). 

## Requirements:
- python 2 or 3
- Biopython
- Seaborn
- NumPy
- pandas
