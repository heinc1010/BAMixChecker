======================================
Input
======================================

BAM files
-----------------

To call variants by running GATK HaplotypeCaller, each BAM file should be indexed.

::
    
    samtools index /path/Tumor_01.bam

BAMixChecker calls variants in GVCF file formats.
Typically, there will be one BAM file with single sample ID, but if the input is a multiple-sample BAM file, it needs to contain the read group informaiton, which can be added with the Picard AddOrReplaceReadGroups tool. (The tool instruction is described at https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.0.0/picard_sam_AddOrReplaceReadGroups.php)

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

.. seealso:: recommended processing for accurate variant discovery using GATK is described at https://software.broadinstitute.org/gatk/best-practices/workflow?id=11165.

To run BAMixChecker, indicating the BAM files directory path with the -d option or a list of BAM files with the -l option is required.

The list can be formatted either as:


* One BAM file on each line. BAMixChecker checks the file names and evaluates whether the files are paired based on the file name.

::
     
    /path/Tumor_01.bam
    /path/Normal_01.bam
    /path/Tumor_02.bam
    /path/Normal_02.bam
    /path/Tumor_03.bam
    /path/Normal_03.bam
    /path/Tumor_04.bam
    /path/Normal_04.bam

.. note:: If you want to compare files only by genotype, you can use the '--OFFFileNameMatching' option.

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

  


.. note:: If the number of files is less than 6 or the file names do not contain common pattern when it is divided by the delimiters, BAMixChecker only matches by genotype, not by file name and skips the creation of the ‘Mismatched_sample.txt’ file, similarly to what happens using the --OFFFileNameMatching option.

RNA-seq BAM file
~~~~~~~~~~~~~~~~~
.. seealso:: Additional proper processing for RNA-Seq data is described on GATK forum: https://gatkforums.broadinstitute.org/gatk/discussion/3891/calling-variants-in-rnaseq.

Reference sequence
------------------------------

To perform the variant calling using GATK HaplotypeCaller requires an indexed reference sequence file ('.fa' and '.fai' format) and '.dict' format. The reference needs to be the same reference, which was used to align the reads.

iGenomes provides 'Ready-To-Use' reference sequence file of various species including human, along with their annotation files.

https://support.illumina.com/sequencing/sequencing_software/igenome.html

Or, you can create the annotation files by your self using samtools and Picard tools.::

    samtools faidx Homo_sapiens.GRCh38.fa

::

    java -jar picard.jar CreateSequenceDictionary \
    R=Homo_sapiens.GRCh38.fa \
    O=Homo_sapiens.GRCh38.dict

.. seealso:: see more details in https://gatkforums.broadinstitute.org/gatk/discussion/1601/how-can-i-prepare-a-fasta-file-to-use-as-reference .


SNP list for non-human organism
---------------------------------------

BAMixChecker has been implemented for human NGS datasets. However, it, also, can be applied to other species using the option --NonHumanSNPlist to provide a customized SNP list and reference.

To extract only informative region, SNPs observed in a large population is required.

It is hard to generalize the method to select informative SNPs because annotation varies across databases.

However, a mandatory annotation is minor allele frequency (MAF) in the population of the species.

In addition to MAF, mappability information is useful to select informative loci.

The recommended SNP loci selection steps are described below:

Filter out uncertain variants from a list of SNP observed in the species of interest. It is recommended to filter variants using mapping quality, quality by depth condition, etc.

Remove SNPs in low mapping region such as regions of low complexity, segmentally duplicated regions, simple repeat regions, etc. Such information may be available from the annotation database or you can get these regions information from e.g the UCSC genome browser simple repeat regions track.

For SNPs located in high mapping rate regions, select SNP loci with only a higher MAF. For human, we apply a global MAF (GMAF) over 0.45 and below 0.55 and a MAF over 0.35 and below 0.65 within each population. If the database does not have MAF information for the study population of interest, GMAF information can be applied. However, we recommend considering higher GMAF condition if the SNP set is too large to produce expected results. This may happen due to lack of proper annotation filtering at an earlier step.

For targeted sequencing dataset, the SNP set should take into consideration not only a higher MAF but also the number of SNPs.

To compare sample genotype, a sufficient number of SNP loci is required for samples comparison.

For human data, we adjusted MAF condition to contain over 200 SNPs for a dataset within the target region information from the BED file. Even though as few as 50 SNPs could be discriminative in RNA-Seq data with 0.45 < MAF < 0.55 and 0.35 < within population < 0.65, we recommend SNPs set to have over 200 loci as the possibility of mutation decreases by decreasing MAF.

These steps cannot be automated for non-human organism due to the variability inherent to every annotation database.

Instead, users can check the number of SNPs in the targeted region using bedtools.

The command is as follows, ::

	bedtools intersect -a SNP_LIST.BED -b TARGETED.BED | wc -l

   
 
    If the number is too small, we recommend adjusting the MAF conditional selection. 
    To reduce variants calling time, we suggest to save the intersected SNPs using the command

 ::

        bedtools intersect -a SNP_LIST.BED -b TARGETED.BED > snp_list.targeted_only.bed

  
Another precaution to ensure that the contigs in the generated SNP list are the same as in the reference.

    ex) SNP list contigs : [chr1, chr2,...] , Reference contigs : [chr1, chr2, ... ]                        .... Working
        SNP list contigs : [1, 2,...] , Reference contigs : [chr1, chr2, ... ]                              .... ERROR

        SNP list contigs : [chr1, chr2,...,chrY, hs37d5] , Reference contigs : [chr1, chr2,...,chrY]        .... ERROR
    
Also, they should be the same as the contigs in BAM files. (If the reference is the same with the one used to align the BAM files, they would have the same reference contigs.)

.. seealso:: Additionally, the user can refer http://evodify.com/gatk-in-non-model-organism/ for the BAM file processing of non-human organisms. 


