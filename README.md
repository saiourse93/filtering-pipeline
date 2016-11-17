# IMSA on Nextflow
The **IMSA on Nextflow** is a modified version of the [Intergrated Metagenomic Sequence Analysis (IMSA)](https://sourceforge.net/projects/arron-imsa/) metagenomics pipeline developed by [Aaron Lab](http://dermatology.ucsf.edu/arronlab/Arron_Lab.html). Through user-defined databases and applications, the IMSA pipeline takes as input reads from high throughput sequencing and filters out the host reads. The remaining reads are then characterised using a comprehensive nucleotide and taxonomy databases, allowing for the identification of microbial/pathogen genomes within host organisms. The **IMSA on Nextflow** pipeline presented has been modified to take advantage of todays powerful computing clusters, which allows for each application in the IMSA pipeline to run with user specified resources. This is particulalry advantageous for when analysing a large number of samples as the jobs are submitted in parallel to each other, thus reducing the overall runtime.

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
2. Comprehensive nucleotide database (preferably `nr`) and BLAST Taxonomy database for sequence characterization (download from [NCBI](ftp://ftp.ncbi.nlm.nih.gov/) `ftp` site)
3. 

# Docker Image
*Coming soon...*

# Configuration
In order to use the **IMSA on Nextflow** pipeline, the following `config.py` python script of IMSA has to be configured with program and database paths required by IMSA to run:
## `config.py`

| Variable | Definition |
| :----- | :----- |
| `SRC_DIRECTORY` | Full path of the directory that contains the `imsa` directory. |
| `BOWTIE2_DATABASES` | Full path to the `bowtie2` ebwt index files. |
| `BLAT_DATABASES` | Full path to the `blat` 2bit index files. |
| `BLAT_OOC_FILES` | Full path to the `blat`  ooc file for the database. |
| `BLAST_DATABASES` | Full path to the BLAST databases. |
| `BLAST_TAX_DB` | Full path the BLAST taxonomy database. |
| `PATH_TO_BOWTIE2` | Full path to your `bowtie2` executable (use program name if installed in you path). |
| `PATH_TO_BLASTN` | Full path to your `blastn` executable (use program name if installed in you path)|
| `PATH_TO_BLAT` | Full path to yout `blat` executable (use program name if installed in you path)|

After configuring the `config.py` script of IMSA, change the first line of every `python` script in the IMSA pipeline folder to point to the default `python` (2.6 or 2.7) interpreter using the following command on terminal (replace </path/to/python> with the full path to your `python`):

````
$ find ./imsa -iname *.py -exec sed -i 's|^#!.*python$|#!<path/to/your/python>|' {} \;
```

# Pipeline Execution
Assuming you have a folder `/home/myhome/fasqdata` in you machine with a set of paired-end reads in `fastq` format, the IMSA on Nextflow pipwline can be executed as follows:

```bash
$ nextflow run main.nf --data /home/myhome/fastqdata --actions actions.txt --out /home/myhome/output
```

Parameters:
`--imsa`: Optional parameter - specifies the loction of the `imsa` folder.
`--data`: Full path to the directory with reads to be analysed.
`--actions`: Action file containing all the actions to be taken in the analysis.
`--out`: Output directory.

# Pipeline Output
