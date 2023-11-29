BASEDIR=$(pwd)  # set pwd as basedir 
echo "$BASEDIR" #/mnt/f/Meng/data_processing

if [ -d "${BASEDIR}"/adapter_clipped ] ; then
    echo "adapter_clipped exist!"
else 
	mkdir ${BASEDIR}/adapter_clipped
fi

if [ -d "${BASEDIR}/bam" ] ; then
	echo "bam exist!"
else
	mkdir ${BASEDIR}/bam
fi

if [ -d "${BASEDIR}/fastqc" ] ; then
	echo "fastqc exist!"
else
	mkdir ${BASEDIR}/fastqc
fi

########################################################################
# fastq_list.txt contains all he prefix names for wgs fastq files. 
# fastq_table.tsv contains all the sample details for each fastq file. 
#######################################################################
SAMPLELIST=$BASEDIR/sample_lists/fastq_list.txt
SAMPLETABLE=$BASEDIR/sample_lists/fastq_table.tsv 

echo $SAMPLELIST
echo $SAMPLETABLE

for SAMPLEFILE in $(cat $SAMPLELIST) 
do 
	SAMPLE_ID=$(grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 4)
	POPULATION=$(grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 5)
	echo $SAMPLEFILE refers to sampel $SAMPLE_ID from $POPULATION 
done	

########################################################################
# process each fastq file and do fastqc  
#######################################################################

RAWFASTQSUFFIX1=_1.fastq.gz # suffix to raw fastq files. use forward reads with paired end data 
RAWFASTQSUFFIX2=_2.fastq.gz
for SAMPLE in $(cat $SAMPLELIST); do 
    echo -e '\n\nstart to run fastqc on the following fastq gz file ... ...'
	echo $SAMPLE
	# zcat each fastq gz file and run fastqc command 
	##zcat $BASEDIR'/raw_fastq/'$SAMPLE$RAWFASTQSUFFIX1 | fastqc stdin:$SAMPLE$RAWFASTQSUFFIX1 --outdir=$BASEDIR'/raw_fastq/'
	##zcat $BASEDIR'/raw_fastq/'$SAMPLE$RAWFASTQSUFFIX2 | fastqc stdin:$SAMPLE$RAWFASTQSUFFIX2 --outdir=$BASEDIR'/raw_fastq/'	
	echo ' '
done
	
################# All programs we will use to do WGS analysis ####################

FASTQC=/programs/bin/fastqc/fastqc # install fastqc in ubuntu, and run fastqc command directly
TRIMMOMATIC=/mnt/f/Meng/data_processing/trimmomatic-0.39.jar
PICARD=/mnt/f/Meng/data_processing/picard.jar
# https://github.com/broadinstitute/picard/releases   I have to use old version 2.7.5 since old version of java 
# install picard-tools in ubuntu: https://www.cyberithub.com/how-to-install-picard-tools-on-ubuntu-20-04-lts-focal-fossa/
SAMTOOLS=/programs/bin/samtools/samtools # install samtools directly in ubuntu
BOWTIEBUILD=/programs/bin/bowtie2/bowtie2-build
BOWTIE=/programs/bin/bowtie2/bowtie2
# use conda install -c bioconda bowtie2
BAMUTIL=/programs/bamUtil/bam
GATK=/mnt/f/Meng/git/lcwgs-guide-tutorial/tutorial1_data_processing/GenomeAnalysisTK.jar ## Note that GATK-3.7 is needed for indel-realignment, which is no longer suppported in new versions
JAVA=/usr/local/jdk1.8.0_121/bin/java ## An older version of java is required for GATK-3.7 to run on Cornell's BioHPC server

##################################################################################
SAMPLELIST=$BASEDIR/sample_lists/fastq_list.txt
SAMPLETABLE=$BASEDIR/sample_lists/fastq_table.tsv 

RAWFASTQDIR=$BASEDIR/raw_fastq/
RAWFASTQSUFFIX1=_1.fastq.gz # suffix to raw fastq files. use forward reads with paired end data 
RAWFASTQSUFFIX2=_2.fastq.gz
ADAPTERS=$BASEDIR/reference/NexteraPE_NT.fa

FASTQDIR=$BASEDIR/adapter_clipped/ # path to the directory where fastq file are stored after adapter clipped 
FASTQSUFFIX1=_adapter_clipped_f_paired.fastq.gz 
FASTQSUFFIX2=_adapter_clipped_r_paired.fastq.gz 
MAPPINGPRESET=very-sensitive # the pre-set option to use for mapping in bowtie2 

REFERENCE=$BASEDIR/REFERENCE/mme_physalia_testdata_chr24.fa
REFBASENAME="${REFERENCE%.*}"

echo "######################################################################"
echo -e '\n\nstart to adapeter clipping using TRIMMOMATIC'
for SAMPLEFILE in $(cat $SAMPLELIST); do 
	echo $SAMPLEFILE
	# extract relevant values from the table of sample, sequencing, lane ID
	SAMPLE_ID=$(grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 4)
	POP_ID=$(grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 5)
	SEQ_ID=$(grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 3)
	LANE_ID=$(grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 2)
	DATATYPE=$(grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 6)
	
	SAMPLE_UNIQUE_ID=$SAMPLE_ID'_'$POP_ID'_'$SEQ_ID'_'$LANE_ID
	
	# the input and output path and file prefix 
	RAWFASTQ_ID=$RAWFASTQDIR$SAMPLEFILE 
	SAMPLEADAPT=$BASEDIR'/adapter_clipped/'$SAMPLE_UNIQUE_ID
	
	## ADAPTER CLIP THE READS WITH TRIMMOMATIC 
	###################################################################################################
	# When the insert length of a library fragment is shorter than the read length, 
	# the sequencer will read into the adapter sequence. This means that the end of 
	# the read will not be from our actual sample, but will be adapter sequence, 
	# which may lead to lower alignment performance and even biases in the result if not removed.
	# we use Trimmomatic to clip the adapter sequence off the ends of reads where they appear. 
	# This step requires us to input the known adapter sequences that we used when preparing the libraries.
	#################################################################################################
	
	if [ $DATATYPE = pe ]; then
		java -jar $TRIMMOMATIC PE -threads 1 -phred33 $RAWFASTQ_ID$RAWFASTQSUFFIX1 $RAWFASTQ_ID$RAWFASTQSUFFIX2 $SAMPLEADAPT'_adapter_clipped_f_paired.fastq.gz' $SAMPLEADAPT'_adapter_clipped_f_unpaired.fastq.gz' $SAMPLEADAPT'_adapter_clipped_r_paired.fastq.gz' $SAMPLEADAPT'_adapter_clipped_r_unpaired.fastq.gz' 'ILLUMINACLIP:'$ADAPTERS':2:30:10:1:true MINLENGTH:40' 
	elif [ $DATATYPE = se ]; then
		java -jar $TRIMMOMATIC SE -threads 1 -phred33 $RAWFASTQ_ID$RAWFASTQSUFFIX1 $SAMPLEADAPT'_adapter_clipped_se.fastq.gz' 'ILLUMINACLIP:'$ADAPTERS':2:30:10 MINLENGTH:40'
    fi
done
################################################################################################################
# The output from Trimmomatic only shows how many full reads get removed, not how much the reads within 
# the file get truncated. In the library with lots of adapter, many of the reads will now be shorter, 
# but as long as they’re still longer than our threshold of 40bp, they will not get removed. 
# If we wanted to know how much sequence we lost from each fastq, counting the number of bases 
# lost is more informative than the number of sequences. It is also always a good idea to check the 
# adapter clipped fastq files with FastQC to make sure that you did in fact get rid of the adapter sequence. 
# We don’t have time to do this in class, but if you’re interested, you can use the fastqc loop above and 
# just change it to run on the files in your adapter_clipped folder.
################################################################################################################

###############################
echo -e  '\n\nbuild reference index files'
###############################

echo -e '\tuse samtools to index reference fasta file '
samtools faidx $REFERENCE
echo -e '\tuse picard to create dict for fasta'
java -jar $PICARD CreateSequenceDictionary R=$REFERENCE o=$REFBASENAME'.dict'
echo -e '\t use bowtie2-build to create bowtie2 index file\n\n'	
bowtie2-build $REFERENCE $REFBASENAME

##############################################################
# map to the reference, sort, and quality filter
# In this step, we align the short reads within each fastq file to the reference genome using bowtie2. 
# The resulting alignment file, in sam format, will be converted to a binary format bam for more efficient storage. 
# Each mapped read will have a mapping quality, which indicates how confident that mapper is that a read is mapped in 
# the correct position. The bowtie2 manual defines it as “a non-negative integer Q = -10 log10 p, 
# where p is an estimate of the probability that the alignment does not correspond to the read’s true point of origin.
# ” Accordingly, a mapping quality (or MAPQ) of 10 or less indicates that there is at least a 1 in 10 chance that 
# the read truly originated elsewhere, and a MAPQ of 20 indicates at least a 1 in 100 chance.
# Here, to only retain reads for which we are reasonably certain have been mapped in the correct place, 
# we will filter out reads with a mapping quality lower than 20, and 
# after that sort the filtered alignment file for easier computation in the next step.
##############################################################

REFERENCE=$BASEDIR/reference/mme_physalia_testdata_chr24.fa # Path to reference fasta file and file name
REFNAME=mme_physalia_testdata_chr24 # Reference name to add to output files, e.g. gadMor2

echo -e '\n\nuse bowtie2 to map reads after adapter clipped'

for SAMPLEFILE in $(cat $SAMPLELIST); do 

	SAMPLE_ID=$(grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 4)
	POP_ID=$(grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 5)
	SEQ_ID=$(grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 3)
	LANE_ID=$(grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 2)
	DATATYPE=$(grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 6)
	SAMPLE_UNIQUE_ID=$SAMPLE_ID'_'$POP_ID'_'$SEQ_ID'_'$LANE_ID
	
	# the input and output path and file prefix 
	SAMPLETOMAP=$FASTQDIR$SAMPLE_UNIQUE_ID
	SAMPLEBAM=$BASEDIR'/bam/'$SAMPLE_UNIQUE_ID
	
	# define platform unit (PU) which is the lane number 
	PU=$(grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 2)
	
	# define reference base name 
	REFBASENAME="${REFERENCE%.*}"
	
	# MAP reads to the reference 
	echo $SAMPLE_UNIQUE_ID
	
	# map the paired-end reads 
	if [ $DATATYPE = pe ]; then 
		# We ignore the reads that get orphaned during adapter clipping because that is typically 
		# a very small proportion of reads. If a large proportion of reads get orphaned 
		# (loose their mate so they become single-end), these can be mapped in a separate step and the resulting 
		# bam files merged with the paired-end mapped reads.
		bowtie2 -q --phred33 --$MAPPINGPRESET -p 1 -I 0 -X 1500 --fr --rg-id $SAMPLE_UNIQUE_ID  \
		        --rg SM:$SAMPLE_ID --rg LB:$SAMPLE_ID --rg PU:$PU --rg PL:ILLUMINA -x $REFBASENAME \
		        -1 $SAMPLETOMAP$FASTQSUFFIX1 -2 $SAMPLETOMAP$FASTQSUFFIX2 \
				-S $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'.sam' 
    
	# Map the single-end reads
    elif [ $DATATYPE = se ]; then
        bowtie2 -q --phred33 --$MAPPINGPRESET -p 1 --rg-id $SAMPLE_UNIQ_ID --rg SM:$SAMPLE_ID \
    		    --rg LB:$SAMPLE_ID --rg PU:$PU --rg PL:ILLUMINA -x $REFBASENAME -U $SAMPLETOMAP$FASTQSUFFIX1 
				-S $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'.sam'
    
	
	fi
	
	# convert to bam file for storage including all the mapped reads 
	samtools view -bS -F 4 $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'.sam' > $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'.bam'
	rm -f $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'.sam'
	
	
	# filter the mapped reads to only retain reads with high mapping quality 
	# filter bam files to remove poorly mapped reads 
	samtools view -h -q 20 $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'.bam' | samtools view -buS - | samtools sort -o $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'_minq20_sorted.bam'

	echo -e '\n\n'
done
 
echo 'merge bam files'
$BASEDIR/scripts/merge_bam.sh $BASEDIR


################## duplicate and clip overlapping read pairs ##################
# Here, we remove the PCR duplicates and trim the overlapping part of each read pair in pair-end data. 
# It is important to wait to deduplicate until after merging, because PCR duplicates for the same sample 
# may exist in different lanes. We use the Picard Tools MarkDuplicates.
# We also want to clip overlapping reads. We will use the BamUtil clipOverlap
# https://github.com/statgen/bamUtil/issues/61
# https://genome.sph.umich.edu/wiki/BamUtil
# https://genome.sph.umich.edu/wiki/BamUtil
#############################################################################################################

BAMLIST=$BASEDIR/sample_lists/merged_bam_list.txt # path to a list of merged bam files 
REFNAME=mme_physalia_testdata_chr24 # reference name to add to output files 
BAMUTIL=/mnt/f/Meng/git/bamUtil/bin/bam

## loop over each sample 

for SAMPLEBAM in $(cat $BAMLIST); do 

	## remove dplicates and print dupstat file 
	# we used to be able to just specify picard.jar. 
	java -jar $PICARD MarkDuplicates I=$BASEDIR'/bam/'$SAMPLEBAM'_bt2_'$REFNAME'_minq20_sorted.bam' \
	        O=$BASEDIR'/bam/'$SAMPLEBAM'_bt2_'$REFNAME'_minq20_sorted_dedup.bam' \
            M=$BASEDIR'/bam/'$REFNAME'_minq20_sorted_depstat.txt' \
            VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true			
			
	## clip overlapping paired end reads. only necessary for paired-end data 
		$BAMUTIL clipOverlap --in $BASEDIR'/bam/'$SAMPLEBAM'_bt2_'$REFNAME'_minq20_sorted_dedup.bam' --out $BASEDIR'/bam/'$SAMPLEBAM'_bt2_'$REFNAME'_minq20_sorted_dedup_overlapclipped.bam' --stats

done

########################## indel realignment #################################################################
# GATK IndelRealigner 
##############################################################################################################

BAMLIST=$BASEDIR/sample_lists/bam_list_dedup_overlapclipped.list 
REFERENCE=$BASEDIR/reference/mme_physalia_testdata_chr24.fa 
REFNAME=mme_physalia_testdata_chr24

for SAMPLEBAM in $(cat $BASEDIR/sample_lists/merged_bam_list.txt); do 
	echo $BASEDIR'/bam/'$SAMPLEBAM'_bt2_mme_physalia_testdata_chr24_minq20_sorted_dedup_overlapclipped.bam' >> $BAMLIST
done

## LOOP OVER EACH SAMPE 

echo -e '\n\nrun GATK IndelRealigner'

cd $BASEDIR/bam/
for SAMPLEBAM in $(cat $BAMLIST); do 
	samtools index $SAMPLEBAM
done 

## realign around in-dels 
# this is done across all samples at once 

## create list of potential in-dels 
echo 'RealignerTargetCreator'
#java -jar $GATK -T RealignerTargetCreator \
#                -R $REFERENCE -I $BAMLIST -o $BASEDIR'/bam/all_samples_for_indel_realigner.intervals' -drf BadMate 
				
## run the indel realigner tool 
echo 'IndelRealigner'
#java -jar $GATK -T IndelRealigner \
#                -R $REFERENCE \ 
#				-I $BAMLIST \
#				-targetIntervals $BASEDIR'/bam/all_samples_for_indel_realigner.intervals' \
#				--consensusDeterminationModel USE_READS \
#				--nWayOut _realigned.bam 
cd ..          


############################################# estimate read depth in bam files ######################

BAMLIST=$BASEDIR/sample_lists/merged_bam_list.txt 
REFNAME=mme_physalia_testdata_chr24 

for SAMPLEBAM in $(cat $BAMLIST); do 
	samtools depth -aa $BASEDIR'/bam/'$SAMPLEBAM'_bt2_'$REFNAME'_minq20_sorted_dedup_overlapclipped.bam' | cut -f 3 | gzip > $BASEDIR'/bam/'$SAMPLEBAM'_bt2_'$REFNAME'_minq20_sorted_dedup_overlapclipped.bam.depth.gz'
done 

















