======================================
Input
======================================

BAM files
-----------------

To call variants by running GATK HaplotypeCaller, each bam file should be indexed.

::
    
    samtools index /path/Tumor_01.bam

BAMixChecker calls variants in GVCF file format which can be called for a a bam file with single sample ID.

If input is multi-sample BAM file, it needs to replace a read group with Picard AddOrReplaceReadGroups.

::

    java -jar picard.jar AddOrReplaceReadGroups \
    I=RNA_T_01.bam \
    O=RNA_T_01.rg_added_sorted.bam \
    SO=coordinate \
    RGID=project \
    RGLB=library \
    RGPL=platform \
    RGPU=machine \
    RGSM=sample

.. seealso:: Additional recommanded processing for accurate variant discovery with GATK is instrcted in https://software.broadinstitute.org/gatk/best-practices/workflow?id=11165.

To run BAMixChecker, indicate the bam files directory path with –d option or a list of bam files with –l option is required.

The acceptable list form can be two types.


* One bam file on each line. BAMixChecker checks the file names and evaluates whether the files are paired based on the file name.

::
     
    /path/Tumor_01.bam
    /path/Normal_01.bam
    /path/Tumor_02.bam
    /path/Normal_02.bam
    /path/Tumor_03.bam
    /path/Normal_03.bam
    /path/Tumor_04.bam
    /path/Normal_04.bam

.. note:: If you want to compare files only by genotype, you can use '--OFFFileNameMatching' option.

* Tab-separated matched files on each line. BAMixChecker matches files based on the pair information instead of file-name-based matching.

::

    /path/Tumor_01.bam  /path/Normal_01.bam
    /path/Tumor_02.bam  /path/Normal_02.bam
    /path/Tumor_03.bam  /path/Normal_03.bam
    /path/Tumor_04.bam  /path/Normal_04.bam


::
 
        /path/Tumor_01.bam	/path/Normal_01.bam	/path/Meta_01.bam
        /path/Tumor_02.bam	/path/Normal_02.bam	/path/Meta_02.bam
        /path/Tumor_03.bam	/path/Normal_03.bam	/path/Meta_03.bam
        /path/Tumor_04.bam	/path/Normal_04.bam	/path/Meta_04.bam

  


.. note:: If the number of files is under 6 or the file names don’t contain common regulation when it is divided by the delimiters, BAMixChecker only matches by genotype, not by file name and skips to make ‘Mismatched_sample.txt’ which is the same using '--OFFFileNameMatching' option.


RNA-seq BAM file
~~~~~~~~~~~~~~~~~

BAMixChecker calls variants in GVCF file format which can be called for a a bam file with single sample ID.

So to run GATK HapplotypeCaller, RNA-seq bam file needs to replace a read group with Picard AddOrReplaceReadGroups.

::

    java -jar picard.jar AddOrReplaceReadGroups \
    I=RNA_T_01.bam \
    O=RNA_T_01.rg_added_sorted.bam \
    SO=coordinate \
    RGID=project \
    RGLB=library \
    RGPL=platform \
    RGPU=machine \
    RGSM=sample

.. seealso:: Additional proper processing for RNA-seq data is instructed in https://gatkforums.broadinstitute.org/gatk/discussion/3891/calling-variants-in-rnaseq.

Reference sequence
------------------------------

To call with GATK HaplotypeCaller, it requires proper reference sequence file with '.fai' file and '.dict' file for the reference which is the **same** reference used to align your bam files.

iGenomes provides 'Ready-To-Use' reference sequence file of various species including human with the annotation files.

https://support.illumina.com/sequencing/sequencing_software/igenome.html

Or, you can create the annotation files by your self with samtools and Picard.::

    samtools faidx Homo_sapiens.GRCh38.fa

::

    java -jar picard.jar CreateSequenceDictionary \
    R=Homo_sapiens.GRCh38.fa \
    O=Homo_sapiens.GRCh38.dict

.. seealso:: see more details in https://gatkforums.broadinstitute.org/gatk/discussion/1601/how-can-i-prepare-a-fasta-file-to-use-as-reference .


SNP list for non-human organism
---------------------------------------

BAMixChecker is developed for a human NGS dataset. However, it, also, can be applied to other species with ‘--NonHumanSNPlist’ option for customized SNP list and proper reference. 

To extract only informative region, SNPs observed in a large population is required. 

It’s hard to generalize the method to select informative SNPs because annotation in each database is various.

However, a mandatory annotation is minor allele frequency (MAF) in the population of the species. 

In addition to MAF, annotation about region affecting mappability is useful to select informative loci. 

Recommended SNP loci selection steps are below: 

    1. Filter out uncertain variants from a list of SNP observed in the interested species. It is recommended to filter variants using mapping quality, quality by depth condition, etc. which the database offers, among calling filter passed variants.

    2. Remove SNPs in a low mapping region like a low complex region, a segment duplicated regions, and a simple repeat region, etc. 
    It can be annotated on the database or you can get the region information in UCSC genome browser for example of simple repeat region. 

    3. For SNPs located in high mapping rate region, select only higher MAF SNP loci. For human, global MAF over 0.45 and under 0.55 and MAF over 0.35 and under 0.65 within each population are applied. 
    If the database doesn’t have MAF information with interested study population, GMAF information only can be applied. 
    However, we recommend to consider higher GMAF condition if the SNP set is too large to produce expected result of the program due to lack of proper filtering annotation from earlier steps. 

    For targeted sequencing dataset, the SNP set is considered not only higher MAF but also the number of SNPs. 

    To compare genotype of samples, enough number of SNP loci to compare is required to compare between samples. 

    For human data, BAMixChecker adjusts MAF condition to contain SNPs over 200 for a dataset with the target region information from BED file. Even 
    
    though SNPs under 50 could be discriminated in RNA-seq with the condition which is global MAFs over 0.45 and under 0.55 and MAFs over 0.35 and under 0.65 within each population, we recommend SNPs set to have over 200 loci because the possibility of mutation is decreased by decreasing MAF. 
    
    These steps can’t be automated for non-human organism because of a uncertainness of each database annotation. 
    
    Instead of it, users can check the number of SNPs in the targeted region with bedtools. 

    The command is following, ::

	bedtools intersect –a SNP_LIST.BED –b TARGETED.BED | wc –l

   
 
    If the number is too small, we recommend adjusting MAF condition. 
    
    
    To reduce calling time, we suggest to give the intersected SNPs creating with a command
    
    ::

        bedtools intersect –a SNP_LIST.BED –b TARGETED.BED > snp_list.targeted_only.bed

  
Another precaution is that the contigs in the generated SNP list should be the same contigs in the reference.

    ex) SNP list contigs : [chr1, chr2,...] , Reference contigs : [chr1, chr2, ... ]                        .... Working
        SNP list contigs : [1, 2,...] , Reference contigs : [chr1, chr2, ... ]                              .... ERROR

        SNP list contigs : [chr1, chr2,...,chrY, hs37d5] , Reference contigs : [chr1, chr2,...,chrY]        .... ERROR
    
Also, it should be the same with contigs in BAM files. (If the reference is the same with the one to align the BAM files, it would be the same with the reference contigs.)

.. seealso:: Additionally, the user can refer http://evodify.com/gatk-in-non-model-organism/ for bam file processing for non-human organism. 


