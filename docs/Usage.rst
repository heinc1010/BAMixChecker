======================================
Usage
======================================

Required arguments
--------------------
::

    -d –-DIR    Directory path of the .BAM files 
    or
    -l --List   A file with the list of files ( The form is referred above in the 'Input' section )

    -r --Ref    Reference file


    additionally for the Targeted sequencing dataset

    -b --BEDfile    Targeted bed file for Targeted sequencing data mode.

Optional arguments
------------------------------

::

    -v –-RefVer ['hg38','hg19'] Default is “hg38”. If the reference is hg19, give this option ‘-r hg19’.

    -o --OutputDIR  Output directory path. BAMixChecker creates the new directory '/BAMixChecker' under the current directory as a default.

    -p --MaxProcess The max number of process. Default = 1

    --FullPATH  Use to report with the full path of the file. BAMixChecker resports with the only file name as a default.

    --RemoveVCF Use this option to remove called VCF files after running

    --OFFFileNameMatching   Use this option to compare files only by genotype.

    -nhSNP --NonHumanSNPlist    A SNP list for non-human organism sample matching check-up. 

    -pld' --Ploidy  Ploidy of sample. Default = 2 for human.


Usage for each data type
---------------------------------------

BAMixChecker runs a mode for WES and RNA-seq as a default without bed file.  

If a bed file is given with –b option, it runs as a targeted sequencing mode.


1)	Whole genome or Whole Exome data or RNA sequencing data 

::
    
    $ python BAMixChecker.py \
    -d /path/aligned/files/ \
    -r /path/reference/HG38/genome.fa \
    -o /path/new/directory 

Or

::

    $ python BAMixChecker.py \
    -l /path/aligned/file_list.txt \
    -p 4 \
    -r /path/reference/HG38/genome.fa \
    -o /path/new/directory


2)	Targeted sequencing data

::

    $ python BAMixChecker.py \
    -d /path/aligned/files/ \
    -p 2 \
    -r /path/reference/HG19/genome.fa \
    -o /path/new/directory \
    -v hg19 \
    -b /path/targeted.bed


Or

::

    $ python BAMixChecker.py \
    -l /path/aligned/file_list.txt \
    -r /path/reference/HG19/genome.fa \
    -o /path/new/directory \
    -v hg19 \
    -b /path/targeted.bed

.. note:: If the dataset consists of both of WES/RNA-seq and Targeted sequencing data mapping with the same reference, run as targeted sequencing data mode with the targeted bed file for the targeted sequencing data.


Tutorial
---------------------------------------

1. Set the required variants::
    
    $ cd BAMixChecker/tutorialData
    $ vim script_example.sh
    
        BAMixChecker_PATH=/path/to/BAMixChecker
        OUT_DIR=/output/dir
        REF=/path/of/reference/hg38/genome.fa


2. Run the script::

    $ sh script_example.sh
