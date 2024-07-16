# annopipe
A bash script for functional annotation (Trinotate) and metabolic pathway prediction (Pathwaytools)

The script is a wrapper for:
Trinotate - https://github.com/Trinotate/Trinotate
emapper2gbk - https://github.com/AuReMe/emapper2gbk
mpwt - https://github.com/AuReMe/mpwt

```
annopipe [OPTIONS]
```

## Options
| Short     | Long      | Description     |
| ------------- | ------------- | -------- |
| -n          | --name         | Name of the run  |
| -g         | --genome         | Genome FASTA file  |
| -f          | --gtff         | Structural annotation file in gtf/gff3 format  |
| -o          | --outdir         | Output directory (Default = Name_annopipe)  |
| -e          | --egg         | Annotation file produced by EggnogMapper (use '--egg path_to_file' for egg2gbk only; taken from Trinotate if omitted)  |
| -p          | --prot         | Proteome file (for egg2gbk)  |
| -s          | --species         | Full species name (e.g. \"Arabidopsis thaliana\") (for egg2gbk)  |
| -r          | --run         | Trinotate analyses to run (swissprot_blastp,swissprot_blastx,pfam,signalp6,tmhmmv2,infernal,EggnogMapper; run all with '--run 0')  |
|	| --bingo         | Create BiNGO-style GO mapping from GOmap file (default when running Trinotate) |
| -gt          | --gff-type         | Type of GFF CDS retrieval passed to emapper2gbk (see README)  |
| -c          | --patho         | Create Pathway tools database using Pathologic  |
| -k          | --keep         | Keep all temp files in outdir  |
| -t           | --threads #        | Number of CPU threads to use (Default = Detected processors or 1)  |
| -h           | --help       | Display help  |


## Examples

### Run Trinotate and Pathologic
```
	annopipe --threads 32 --rm consensi.fa.classified --run 0 --name A.thaliana --genome TAIR10_genome.fa --outdir At_anno
```

### Skip functional annotation and use a protein file to create a gene-bank file with eggnog data
```
	annopipe --threads 32 --name A.thaliana --genome TAIR10_genome.fa --gtff TAIR10.gff3 --prot At.proteins.fa --egg At.eggnog.annotations --outdir At_anno
```
### Create functional annotation with BLASTp and BLASTx alone, and create a gene-bank file with that data
```
	annopipe --threads 32 --name A.thaliana --genome TAIR10_genome.fa --gtff TAIR10.gff3 --egg 0 --run swissprot_blastp,swissprot_blastx,EggnogMapper --outdir At_anno
```

## Comments
- Please read https://github.com/AuReMe/emapper2gbk for compatibility of GFF, Protein and eggnog files. This script runs with `-gt mRNA` by default, change this with `-gt | --gff-type`.
- The `-n | --name`, `-g | --genome` and `-f | --gtff` are mandatory.
- Excluding `-r | --run` skips functional annotation entirely.
- The `-p | --prot` option uses the designated file instead of the one produced by Trinotate. To use the proteome file from Trinotate, just exclude `-p | --prot`.
