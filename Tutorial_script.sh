######################################################################################################
## Tutorial samples were extracted from SRP133424 "Gastric adenoma and carcinoma genomes" from SRA. ##
######################################################################################################

#Command to run the tutorial data 
#Normal
python /path/to/MuBaMer/MuBaMer.py \
-d /path/to/MuBaMer/tutorialData/normal \
-o /output/dir/normal \
-r /hg38/reference/genome.fa \
-b /path/to/MuBaMer/tutorialData/target.bed

#Mismatched sample exist
#'S1254_N.bam' and 'S1345_T.bam ' are swapped, 'S1983_N.bam' and 'S1983_T.bam are unpaired.
python /path/to/MuBaMer/MuBaMer.py \
-d /path/to/MuBaMer/tutorialData/mismatched \
-o /output/dir/mismatched \
-r /hg38/reference/genome.fa \
-b /path/to/MuBaMer/tutorialData/target.bed
