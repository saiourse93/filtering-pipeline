#!/opt/exp_soft/python27/bin/python
#

########################  The following variables will have to be changed for all installations #######

# this is the full path of the directory that contains the 'imsa' directory with your code.  
SRC_DIRECTORY = "/home/phelelani/applications/"

# BOWTIE_DATABASES specifies where the ebwt index lives for each bowtie database you want to use in
# IMSA.  The path should be the full path that you would use on a bowtie command line when referencing
# the ebwt file.  BOWTIE_DB_ABBREV links the name of the database to the short code used internally.
BOWTIE_DATABASES = {"human":"/global/blast/reference_genomes/Homo_sapiens/UCSC/hg19/Sequence/BowtieIndex/genome"}
BOWTIE_DB_ABBREV = {"human":"hg"}

# BOWTIE2 DATABASES
BOWTIE2_DATABASES = {"human":"/global/blast/reference_genomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome"}
BOWTIE2_DB_ABBREV = {"human":"hg"}

# BLAT_DATABASES specifies where the 2bit blat index lives for each blat database you want to use in
# IMSA.  The BLAT_OOC_FILES list the ooc file for the database, if it exists.
BLAT_DATABASES = {"human":"/global/blast/reference_genomes/Homo_sapiens/UCSC/hg19/Sequence/BlatIndex/genome.2bit"} 

BLAT_OOC_FILES = {"human":"/global/blast/reference_genomes/Homo_sapiens/UCSC/hg19/Sequence/BlatIndex/genome.11.ooc"}
BLAT_DB_ABBREV = {"human":"hg"}

# BLAST_DATABASES list where each blast database lives for the databases you want to use within imsa.
BLAST_DATABASES = {"nt":"/global/blast/blast_dbs/nt", 
                   "human":"/global/blast/reference_genomes/Homo_sapiens/UCSC/hg19/Sequence/BlastDB/genome.fa"}
BLAST_DB_ABBREV = {"nt":"nt", "human":"hg", 
                   "humanRNA":"hr", "humanRNA-1":"hr1", "humanRNA-2":"hr2", "humanRNA-3":"hr3",
                   "humanRepeat":"rep"}
BLAST_TAX_DB = "/global/blast/gi_taxid_nucl.dmp"


###################  The following variables may need to be changed, depending on your computer configuration  ###############

# The PATH_TO variables are what you would type on the command line to run the following programs.
# If the programs have been installed in your path, you'll just need the program name, as shown below.
PATH_TO_BOWTIE = "bowtie"
PATH_TO_BOWTIE2 = "bowtie2"
PATH_TO_BLASTN = "blastn"
PATH_TO_BLAT = "blat"


# PIPELINE_DIRECTORY is the name of the python module.  It will be imsa unless you change the
# name of the directory after unpacking the imsa tarball
PIPELINE_DIRECTORY = "imsa"

# USE_GRID should be true if you have SUN GRID ENGINE installed.
USE_GRID = False

# The default number of processors to use if no value is specified
DEFAULT_NUM_PROCESSORS = 12

# The delimiter is used when analyzing paired end files to find the pairs.
DEFAULT_DELIMITER = " "

# The quality offset depends on how the fastq file was generated
DEFAULT_QUALITY_OFFSET = 33

###################  The following variables will not need to be changed in most cases ###############


RANK_TO_COLOR = {"species":'"#bae4b3"', "genus":'"#74c476"', "family":'"#31a354"', "division":'"#006d2c"'}
MAX_NODE_SIZE = 3.0
MIN_NODE_SIZE = 1.5
MAX_NODE_SIZE_CUTOFF = 0.2
WIDTH_TO_HEIGHT = 3.0

# to make it easier to override settings, put each setting into a list
# ex: DEFAULT_BLAT_SETTINGS = ["-minIdentity=80", "-fastMap"]
DEFAULT_BLAT_PARAMS = ["-minIdentity=80"]
# If nothing is specified for BLAT, right now a match with 80% identity
#    across the ENTIRE query (total percent) is required.  Change in the code
#    below if you want a different default action, like minCoverage.
DEFAULT_BLAT_THRESHOLD = 0.8  #the total percent id across the whole query length is 80 


# to make it easier to override settings, put each setting into a list
# ex: DEFAULT_BLAST_SETTINGS = ["-task blastn", "-ungapped", "-perc_identity 0.8"]
# only return a maximum of 5 hits per sequence!  Awesome -- this should speed things up considerably!
DEFAULT_BLAST_PARAMS = ["-max_target_seqs 5"]
#DEFAULT_BLAST_PARAMS = [""]
# currently the BLAST_THRESHOLD is applied to the EVAL.  To have a different action,
# change the code in the Blast function below where this constant is used.
DEFAULT_BLAST_THRESHOLD = 0.0001  #TBD

DEFAULT_LZW_THRESHOLD = 0.5
DEFAULT_REMOVE_N = 5
DEFAULT_HOMOPOLYMER_THRESHOLD = 15

DEFAULT_QUALITY_NUM_BASES = 3
DEFAULT_QUALITY_THRESHOLD = 15

DEFAULT_ACTIONS = """removeN 
blat human -fastMap
blat human
blast human
lzw 0.5"""


