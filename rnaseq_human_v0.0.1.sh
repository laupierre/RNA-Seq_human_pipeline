#!/bin/bash

#### This is the RNA-Seq for the mouse species (v.0.0.1dev)

# clean the log.out file
if [ -f "log.out" ] ; then
    rm log.out
fi


# Call getopt to validate the provided input. 
options=$(getopt -o brg --long method: -- "$@")
# Check status code
[ $? -eq 0 ] || { 
    echo "Incorrect options provided"
    exit 1
}

eval set -- "$options"

while true; do
    case "$1" in
    --method)
        shift; # The arg is next in position args
        COLOR=$1
        [[ ! $COLOR =~ ^(star|kallisto|salmon)$ ]] && {
            echo "Incorrect options. Please provide a correct method: star or salmon"
            exit 1
        }
        ;;
    --)
        shift
        break
        ;;
    esac
    shift
done


#################
#### Variables (this folder is not bound by apptainer)
CONTAINER=/projects/ncrrbt_share_la/dev_pipe
INDEX=/projects/ncrrbt_share_la/dev_pipe2
CPUS=12

#### Start message
echo "The RNA-Seq pipeline (v.0.0.1_dev) is for the human species and the $COLOR method was chosen" >> log.out
echo "Starting the RNA-Seq pipeline on `date`" >> log.out
##################



###################
#### BBMap analysis
echo "Starting ribosomal and mitochondrial filtering with BBMap ..." >> log.out

cp /projects/ncrrbt_share_la/dev_pipe2/human_ribosomal.fa .
mkdir projects

var=(`ls *_R1*.fastq.gz`)

	for i in ${var[@]}
	do
	read2=`echo ${i} | sed 's/_R1/_R2/g'`
	prefix=`echo ${i%%_R1*}`
	
	ls -l $prefix\_*filtered.fastq.gz > /dev/null
	
	if [ `echo $?` -eq 0 ]; then
	#echo  "exists"
	continue
	else
	
	apptainer exec $CONTAINER/bbmap.sif /bin/bash -c \
   	"bbduk.sh threads=$CPUS in=$i in2=$read2 out1=$prefix\_R1_001.filtered.fastq out2=$prefix\_R2_001.filtered.fastq ref=human_ribosomal.fa k=31 overwrite=t"
		
	## Put this inside the loop
	if [ $? -eq 0 ]
	then
    	echo "bbmap processed sample ${prefix}" >> log.out
	else
	echo "bbmap failed on sample ${prefix}. Pipeline terminated"  >> log.out
	exit 1
	fi
	fi
	
	done
	

	ls *filtered.fastq | parallel -j 2 pigz -p 6 {}
	mv *filtered.fastq.gz projects

cd projects
mv ../log.out .

AFTER=`date`
echo "BBMap finished on ${AFTER}" >> log.out
####




##################
#### STAR analysis
if [ "$COLOR" = "star" ]; then
echo "Starting alignment with STAR ..." >> log.out

cp -r $INDEX/human_star_index .

var=(`ls *_R1*.fastq.gz`)

	for i in ${var[@]}
	do
	read2=`echo ${i} | sed 's/_R1/_R2/g'`
	prefix=`echo ${i%%_R1*}`
	
	apptainer exec $CONTAINER/star.sif /bin/bash -c \
	"STAR --genomeDir human_star_index --runThreadN $CPUS \
	--readFilesIn $i $read2 \
	--outSAMtype BAM Unsorted \
	--readFilesCommand zcat \
	--outFileNamePrefix star.${prefix}"

	## Put this inside the loop
	if [ $? -eq 0 ]
	then 
    	echo "STAR processed sample ${prefix}" >> log.out
	else
	echo "STAR failed on sample ${prefix}. Pipeline terminated"  >> log.out
	exit 1
	fi
	done
	

mkdir star_results
mv star.* star_results

AFTER=`date`
echo "STAR finished on ${AFTER}" >> log.out



#############################
#### featureCounts after STAR
echo "Starting counting with featureCounts ..." >> log.out
cd star_results
files=`ls -d *bam | xargs -n1000`

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz
gunzip gencode.v43.annotation.gtf.gz


apptainer exec $CONTAINER/featurecounts.sif /bin/bash -c \
"featureCounts -B -C -s 2 -p --countReadPairs -T $CPUS -t exon -g gene_id --extraAttributes gene_name,gene_type \
-a gencode.v43.annotation.gtf \
-o subread.counts.txt $files" || { 
echo "featureCounts has an error. Pipeline terminated" >> log.out
exit 1
}
	
cd ..

AFTER=`date`
echo "featureCounts finished on ${AFTER}" >> log.out



########################################
#### RSeQC analysis and sambamba markdup
echo "Starting RSeQC and sambamba ..." >> log.out

mkdir rseqc_results
cd rseqc_results

wget https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg38_GENCODE_V42_Comprehensive.bed.gz
gunzip hg38_GENCODE_V42_Comprehensive.bed.gz

var=(`ls ../star_results/*bam`)

	for i in ${var[@]}
	do
	prefix=`echo ${i%%_S*}`
	prefix2=`echo ${prefix##*/}`
	apptainer exec $CONTAINER/sambamba.sif /bin/bash -c \
	"sambamba markdup -t $CPUS $i $prefix2.markdup.bam > markdup.$prefix2.log 2>&1" || { 
   	echo "Sambamba has an error. Pipeline terminated" >> log.out
    	exit 1
	}
	done

var=(`ls *.bam`)	
	
	for i in ${var[@]}
	do
	prefix=`echo ${i%%.bam}`
	apptainer exec $CONTAINER/rseqc.sif /bin/bash -c \
	"infer_experiment.py -r hg38_GENCODE_V42_Comprehensive.bed -i $i 1> rseqc.$prefix.infer_experiment.txt" || { 
   	echo "RSeQC has an error. Pipeline terminated" >> log.out
    	exit 1
	}
	done

cd ..

AFTER=`date`
echo "RSeQC and sambamba finished on ${AFTER}" >> log.out

fi




######################################
#### salmon analysis in automatic mode (should fall back on ISR mode)
if [ "$COLOR" = "salmon" ]; then
echo "Starting salmon pseudo-counting ..." >> log.out

cp -r $INDEX/human_salmon_index .

var=(`ls *_R1*.fastq.gz`)

	for i in ${var[@]}
	do
	read2=`echo ${i} | sed 's/R1/R2/g'`
	prefix=`echo ${i%%_R1*}`
	
	apptainer exec $CONTAINER/salmon.sif /bin/bash -c \
   	"salmon quant -i human_salmon_index -p $CPUS -l A --validateMappings -o salmon.${prefix} -1 $i -2 $read2"
	
	## Put this inside the loop
	if [ $? -eq 0 ]
	then
    	echo "salmon processed sample ${prefix}" >> log.out
	else
	echo "salmon failed on sample ${prefix}. Pipeline terminated"  >> log.out
	exit 1
	fi
	done

mkdir salmon_results
mv salmon.IIT* salmon_results

AFTER=`date`
echo "salmon finished on ${AFTER}" >> log.out

fi
####




####################
#### FastQC analysis
echo "Starting FastQC ..." >> log.out

files=`ls *fastq.gz | xargs -n1000`
mkdir fastqc_results

apptainer exec $CONTAINER/fastqc.sif /bin/bash -c \
	"fastqc -t $CPUS -o fastqc_results $files"

AFTER=`date`
echo "FastQC finished on ${AFTER}" >> log.out
####




####################
#### MultiQC analysis
echo "Starting MultiQC ..." >> log.out

if [ "$COLOR" = "star" ]; then
apptainer exec $CONTAINER/multiqc.sif /bin/bash -c \
"multiqc -f -n multiqc_report_rnaseq \
-m featureCounts $PBS_O_WORKDIR/projects/star_results/*summary \
-m star $PBS_O_WORKDIR/projects/star_results/*Log.final.out \
-m sambamba $PBS_O_WORKDIR/projects/rseqc_results/markdup.star.IIT*.log \
-m rseqc $PBS_O_WORKDIR/projects/rseqc_results/*infer_experiment.txt \
-m fastqc $PBS_O_WORKDIR/projects/fastqc_results/*zip"
fi

if [ "$COLOR" = "salmon" ]; then
apptainer exec $CONTAINER/multiqc.sif /bin/bash -c \
"multiqc -f -n multiqc_report_rnaseq \
-m salmon $PBS_O_WORKDIR/projects/salmon_results/* \
-m fastqc $PBS_O_WORKDIR/projects/fastqc_results/*zip"
fi



###########################
### Write results to output

if [ "$COLOR" = "star" ]; then
cd ..
cp $INDEX/star_wrap.R .
cp $INDEX/gencode.v43.annotation.txt .

apptainer exec $CONTAINER/R.sif /bin/bash -c \
"Rscript star_wrap.R"
fi





####
