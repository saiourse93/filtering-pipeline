# IMSA on Nextflow
The **IMSA on Nextflow** is a modified version of the [Intergrated Metagenomic Sequence Analysis (IMSA)](https://sourceforge.net/projects/arron-imsa/) metagenomics pipeline developed by [Aaron Lab](http://dermatology.ucsf.edu/arronlab/Arron_Lab.html). The IMSA pipeline takes as input reads from high throughput sequencing and filters out exogeneous sequences in a host-genomic background. Assembly pipeline for exogenous reads filtered from human RNA-seq data.

# Dependencies
The **IMSA on Nextflow** pipeline depends on the following programs and databases:
## Programs:
1. Modified `IMSA pipeline` scripts (available on this repository)
2. [`Python 2.6 - 2.7`](https://www.python.org/)
2. [`bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
3. [`blat`](http://hgdownload.cse.ucsc.edu/downloads.html)
4. [`blastn` (NCBI-BLAST+ )](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
5. [`nextflow`](https://www.nextflow.io/)
5. [`Trinity`](https://github.com/trinityrnaseq/trinityrnaseq/wiki)

## Databases and Reference Sequences:
1. Reference genome, including:
   * `bowtie2` index
   * `blat` index (and associated ooc file)
   * `blastn` index
2. Comprehensive database for sequence characterization (nr database from NCBI)

# Docker Image
*Coming soon...*

# Configuration
In order to use the **IMSA on Nextflow** pipeline, the following `config.py` python script of IMSA has to be configured with program and database paths required by IMSA to run:
## `config.py`

After changing the line, edit the following variables in the `config.py` python script.

| Variable | Definition |
| :----- | :----- |
| `SRC_DIRECTORY` | Full path of the directory that contains the `imsa` directory. |
| `BOWTIE2_DATABASES` | Full path to the `bowtie2` ebwt index files. |
| `BLAT_DATABASES` | Full path to the `blat` 2bit index files. |
| `BLAT_OOC_FILES` | Full path to the `blat`  ooc file for the database. |
| `BLAST_DATABASES` | Full path to the BLAST databases. |
| `BLAST_TAX_DB` | Full path the BLAST taxonomy database. |
| `PATH_TO_BOWTIE2` | R |
| `PATH_TO_BLASTN` | R |
| `PATH_TO_BLAT` | R |
| `PIPELINE_DIRECTORY` | R |

After configuring the `config.py` script of IMSA, change the first line of every `python` script in the IMSA pipeline folder to point to the default `python` (2.6 or 2.7) interpreter. Using terminal, "`cd`" into the `imsa` folder packaged with this **IMSA on Nextflow** pipeline and run the following command (replace </path/to/python> with the full path to your `python`):

`for script in $(ls *.py); do sed -i ’s|"#!/opt/exp_soft/python27/bin/python"|"#!<path/to/your/python>"|g’ $script; done`


## ```main.nf```

# Pipeline Execution

# Pipeline Output
