
BASEDIR=/mnt/f/Meng/genotype_snp_calling
cd $BASEDIR

if [ ! -d "${BASEDIR}"/results ]; then 
	mkdir results
fi

DATA=$BASEDIR/bam
REF=$BASEDIR/reference/Ref.fa 
ANC=$BASEDIR/reference/outgrp_ref.fa 
ANGSD=/mnt/f/Meng/git/angsd/angsd
SAMTOOLS=/usr/bin/samtools

$SAMTOOLS faidx $REF 

############################################################
#The workflow is roughly divided into four steps:
#Data filtering and I/O
#Genotype likelihoods
#Genotype calling
#SNP calling
############################################################

$ANGSD --help

cat $BASEDIR/sample_lists/ALL_bams.txt 
echo -e '\n'
wc -l $BASEDIR/sample_lists/ALL_bams.txt 
echo -e '\n'
ls $BASEDIR/sample_lists/*_bams.txt 


######################## step 1 convert DATA to DATA' for SNP calling ############
# with -b we give the file including paths to all BAM files we need to analyse, 
# -ref specifies the reference sequence, 
# -out states the prefix for all output files that will be generated.
# These filters will retain only uniquely mapping reads, not tagged as bad, considering only proper pairs, 
# without trimming, and adjusting for indel/mapping (as in samtools). 
# -C 50 reduces the effect of reads with excessive mismatches, 
# while -baq 1 computes base alignment quality as explained here (BAQ) used to rule out false SNPs close to INDELS.
# Parameter | Meaning |
# --- | --- |
# -minInd 5 | use only sites with data from at least N individuals |
# -setMinDepth 7 | minimum total depth |
# -setMaxDepth 30 | maximum total depth |
####################################################################################

$ANGSD -b $BASEDIR/sample_lists/ALL_bams.txt -ref $REF -out $BASEDIR/results/ALL \
       -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
	   -minMapQ 20 -minQ 20 -minInd 5 -setMinDepth 7 -setMaxDepth 30 -doCounts 1 


######################## step 2 calculate genotype likelihoods for each site at each individual. ############


cd $BASEDIR
echo -e '\n\n genotype likelihood model 2 for PANY bam files '
$ANGSD -b $BASEDIR/sample_lists/PANY_bams.txt -ref $REF -out $BASEDIR/results/PANY \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 \
        -minMapQ 20 -minQ 20 -minInd 5 -setMinDepth 7 -setMaxDepth 30 -doCounts 1 \
        -GL 2 -doGlf 1
		
##########
# Have a look at .glf.gz file. The first two columns are the reference sequence (e.g. chromososome) and position. 
# Then you have 10 likelihoods for all possible genotypes in the order AA,AC,AG,AT,CC,CG,CT,GG,GT,TT. 
# This set of 10 likelihoods is repeated sequentially starting from the left of the file for each individual in the row 
# order of individuals in the BAM file. The values are log-scaled likelihood ratios, all scaled by the most likely genotype.
# Since the likelihoods have been scaled to the most likely and log-transformed, the most likely genotype will have a value of 0.
############


######################## step 3 Genotype calling ############
# In ANGSD, the option to call genotypes is -doGeno:
# $ANGSD -doGeno
# -doGeno 0
#        1: write major and minor
#        2: write the called genotype encoded as -1,0,1,2, -1=not called
#        4: write the called genotype directly: eg AA,AC etc
#        8: write the posterior probability of all possible genotypes
#        16: write the posterior probability of called genotype
#        32: write the posterior probabilities of the 3 gentypes as binary
#        -> A combination of the above can be choosen by summing the values, EG write 0,1,2 types with majorminor as -doGeno 3
#        -postCutoff=0.333333 (Only genotype to missing if below this threshold)
#        -geno_minDepth=-1       (-1 indicates no cutof)
#        -geno_maxDepth=-1       (-1 indicates no cutof)
#        -geno_minMM=-1.000000   (minimum fraction af major-minor bases)
#        -minInd=0       (only keep sites if you call genotypes from this number of individuals)
#        NB When writing the posterior the -postCutoff is not used
#        NB geno_minDepth requires -doCounts
#        NB geno_maxDepth requires -doCounts
# Therefore, if we set -doGeno 2, genotypes are coded as 0,1,2, as the number of alternate/minor alleles. 
# If we want to print the major and minor alleles as well then we set -doGeno 3.

# To calculate the posterior probability of genotypes we need to define a model.
# $ANGSD -doPost
# -doPost 0       (Calculate posterior prob 3xgprob)
#        1: Using frequency as prior
#        2: Using uniform prior
#        3: Using SFS as prior (still in development)
#        4: Using reference panel as prior (still in development), requires a site file with chr pos major minor af ac an
# -doPost 2 uses a uniform prior.

# Furthermore, this calculation requires the specification of how to assign the major and minor alleles (if biallelic).
# $ANGSD -doMajorMinor
#        -doMajorMinor   0
#        1: Infer major and minor from GL
#        2: Infer major and minor from allele counts
#        3: use major and minor from a file (requires -sites file.txt)
#        4: Use reference allele as major (requires -ref)
#        5: Use ancestral allele as major (requires -anc)
#        -rmTrans: remove transitions 0
#        -skipTriallelic 0
##############################################################################################################################

cd $BASEDIR

$ANGSD -glf $BASEDIR/results/PANY.glf.gz -fai $REF.fai -nInd 15 -out $BASEDIR/results/PANY \
    -doMajorMinor 1 -doGeno 3 -doPost 2 -doMaf 1

# The geno.gz output columns are: chromosome, position, major allele, minor allele, genotypes is 0,1,2 format.


# You can control how to set missing genotype when their confidence is low with -postCutoff. 
# For instance, we can set as missing genotypes when their (highest) genotype posterior probability is below 0.95:

$ANGSD -glf $BASEDIR/results/PANY.glf.gz -fai $REF.fai -nInd 15 -out $BASEDIR/results/PANY \
        -doMajorMinor 1 -doGeno 3 -doPost 2 -doMaf 1 -postCutoff 0.95

# how many sites do we have in total? how many sites have at least one missing genotype now?
zcat $BASEDIR/results/PANY.geno.gz | wc -l
zcat $BASEDIR/results/PANY.geno.gz | grep -1 - | wc -l

# Why are there some many sites with missing genotypes?
# The mean depth per sample is around 1-2X, therefore genotypes cannot be assigned with very high confidence. 
# Sites where all genotypes are missing are skipped in the output file.

# In ANGSD we can restrict our analyses on a subset of positions of interest using the -sites option. 
# The file with these positions need to be formatted as (chromosome positions).
# echo Mme_chr24:2558528-4558528 48 > $BASEDIR/results/snp.txt
# echo Mme_chr24:2558528-4558528 61 >> $BASEDIR/results/snp.txt
# We need to index this file in order for ANGSD to process it.
# $ANGSD sites index $BASEDIR/results/snp.txt


####################### allele frequency #############################################################
# In other words, at each site we want to to estimate (or count) how many copies of different alleles 
# (two in case of biallelic variants) we observe in our sample (across all sequenced individuals). 
# However with low depth data direct counting of individually assigned genotypes can lead to biased allele frequencies.

# ANGSD has an option to estimate allele frequencies taking into account data uncertainty from genotype likelihoods:

#$ANGSD -doMaf
# -doMaf  0 (Calculate persite frequencies '.mafs.gz')
#        1: Frequency (fixed major and minor)
#        2: Frequency (fixed major unknown minor)
#        4: Frequency from genotype probabilities
#        8: AlleleCounts based method (known major minor)
#        NB. Filedumping is supressed if value is negative
# -doPost 0       (Calculate posterior prob 3xgprob)
#        1: Using frequency as prior
#        2: Using uniform prior
#        3: Using SFS as prior (still in development)
#        4: Using reference panel as prior (still in development), requires a site file with chr pos major minor af ac an
#Filters:
#        -minMaf         -1.000000       (Remove sites with MAF below)
#        -SNP_pval       1.000000        (Remove sites with a pvalue larger)
#        -rmTriallelic   0.000000        (Remove sites with a pvalue lower)
#Extras:
#        -ref    (null)  (Filename for fasta reference)
#        -anc    (null)  (Filename for fasta ancestral)
#        -eps    0.001000 [Only used for -doMaf &8]
#        -beagleProb     0 (Dump beagle style postprobs)
#        -indFname       (null) (file containing individual inbreedcoeficients)
#        -underFlowProtect       0 (file containing individual inbreedcoeficients)
# NB These frequency estimators requires major/minor -doMajorMinor
# Therefore, the estimation of allele frequencies requires the specification of how to assign the major and minor alleles (if biallelic).

# $ANGSD -doMajorMinor
#        -doMajorMinor   0
#        1: Infer major and minor from GL
#        2: Infer major and minor from allele counts
#        3: use major and minor from a file (requires -sites file.txt)
#        4: Use reference allele as major (requires -ref)
#        5: Use ancestral allele as major (requires -anc)
#        -rmTrans: remove transitions 0
#        -skipTriallelic 0

$ANGSD -b $BASEDIR/sample_lists/PANY_bams.txt -ref $REF -out $BASEDIR/results/PANY \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 \
        -minMapQ 20 -minQ 20 -minInd 5 -setMinDepth 7 -setMaxDepth 30 -doCounts 1 \
        -GL 1 -doGlf 1 -doMajorMinor 1 -doMaf 1





