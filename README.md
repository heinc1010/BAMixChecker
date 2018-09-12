MuBaMer:  
Rapid and easy-to-use Multiple Bam files identifier to detect a pair of samples from each individual.
=================================================================================


MuBaMer is a rapid and easy-to-use multiple bam files identifier to detect a pair of samples from each individual.

It’s simple and fast but robust to detect a pair of bam files from same individual with WES, Targeted sequencing dataset. 

And the user can catch the information of swapped or unpaired pairs information as well as the matched pairs information at a glance.




Requied tools
------------------

#### GATK >= 4.0 ( java 8 required )

#### Bedtools

- The tool PATH is needed to be set on “MuBaMer.config” file.

#### Python 2.7 

packages

- scipy.stats

- numpy

- multiprocessing
   

Instruction
-------------
First, clone the repository

```$ git clone https://github.com/heinc1010/MuBaMer```

Or


Download the compressed file.

```$ unzip MuBaMer-master.zip```

And then set the tools PATH on the configuration file, 

```$ cd MuBaMer
  $ vim MuBaMer.config

  GATK=/path/gatk
  BEDTOOLS=/path/bedtools

  *default (If you added the path in $PATH, you don’t need to modify the configuration file )
  GATK=gatk
  BEDTOOLS=bedtools
  ```

Input
-----------
To run MuBaMer, directory path of bam files with –d option or a list of bam files with –l is required.
A pairs of sample list or directory information is needed.


It only gets the pairs of samples. 

The form of list can be two types.


-Tab-divided paired files on each line. MuBaMer will pair files based on the pair information.
  ```	
  /path/Tumor_1.bam /path/Normal_1.bam
  /path/Tumor_2.bam /path/Normal_2.bam
  /path/Tumor_3.bam /path/Normal_3.bam
  /path/Tumor_4.bam /path/Normal_4.bam
```

-One bam file on each line.  MubaMer will check the file names and evaluate whether the files are pair based on the name.
```
  /path/Tumor_1.bam
  /path/Normal_1.bam
  /path/Tumor_2.bam
  /path/Normal_2.bam
  /path/Tumor_3.bam
  /path/Normal_3.bam
  /path/Tumor_4.bam
  /path/Normal_4.bam
```

#### If the number of files is not above 6 or the file names don’t contain detectable regulation, it only pairs by genotype, not by name and skip to make ‘Unmatched_pair.txt’

Usage
------
#### Required arguments
```
-d –-DIR Directory path of the .BAM files 
or
-l --List A file with the list of files

-r --Ref Reference file
```
#### Optional arguments
```
-b --BEDfile  Targeted bed file for Targeted sequencing data mode.
-v –-RefVer ['hg38','hg19']. default is “hg38”.If the reference is hg19, give this option ‘-r hg19’.
-o --OutputDIR Output directory path
-p --MaxProcess The max number of process. Default = 4
--RemoveVCF Use this option to remove called VCF files after running
```

#### MuBaMer runs for Whole Exome data mode as a default without bed file.  If a bed file is given with –b option, it runs as a  targeted mode.


#### 1)	Whole Exome data
```
  $ python MubaMer.py \
  –d /path/aligned/files/ \
  -r /path/reference/HG38/genome.fa \
  –o /path/new/directory 
  ```
 Or
  ```
  $ python MubaMer.py \
  –l /path/aligned/file_list.txt \
  –p 3 \
  -r /path/reference/HG38/genome.fa \
  –o /path/new/directory
```

#### 2)	Targeted sequencing data. 
 ```
  $ python MubaMer.py \
  –d /path/aligned/files/ \
  –p 2 \
  -r /path/reference/HG19/genome.fa \
  –o /path/new/directory \
  –v hg19 \
  -b /path/targeted.bed
```
Or 
```
  $ python MubaMer.py \
  –l /path/aligned/file_list.txt \
  -r /path/reference/HG19/genome.fa \
  –o /path/new/directory \
  –v hg19 
  -b /path/targeted.bed
```

Output
--------
#### “Unmatched_pair.txt
  ```
  #Pair is not matched by name but by genotype.
  /path/Tumor_2.gvcf /path/Normal_3.gvcf 0.98
  /path/Tumor_3.gvcf /path/Normal_2.gvcf 1
  
  #Samples are matched with nothing by genotype.
  /path/Tumor_4.gvcf 
	  -> pair by name with /path/Normal_4.gvcf (score : 0.43 )
  /path/Normal_4.gvcf 
	  -> pair by name with /path/Tumor_4.gvcf (score : 0.43 )
```

#### “Matched_pair.txt”
  ```
  #Matched pair by genotype and name.
  /path/Tumor_1.gvcf /path/Normal_1.gvcf 
  ```
  
#### “Total_result.txt”
  ```
  #Matched pair : 1, Swapped pair : 2, Unpaired sample : 2
  /path/Normal_1.gvcf /path/Normal_2.gvcf 0.3434 unmatched
  /path/Normal_1.gvcf /path/Normal_3.gvcf 0.3123 unmatched 
  /path/Normal_1.gvcf /path/Normal_4.gvcf 0.423 unmatched
  /path/Normal_1.gvcf /path/Tumor_1.gvcf 1 matched
  /path/Normal_1.gvcf /path/Tumor_2.gvcf 0.1244 unmatched
  /path/Normal_1.gvcf /path/Tumor_3.gvcf 0.3223 unmatched 
  …..
  /path/Tumor_3.gvcf /path/Normal_4.gvcf 0.379 unmatched
  /path/Tumor_3.gvcf /path/Tumor_4.gvcf 0.367 unmatched
  /path/Normal_4.gvcf /path/Tumor_4.gvcf 0.43 unmatched 
```

Contact
--------
heinc1010@gmail.com
