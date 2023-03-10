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
echo "The RNA-Seq pipeline (v.0.1_dev) is for the human species and the $COLOR method was chosen" >> log.out
echo "Starting the RNA-Seq pipeline on `date`" >> log.out
##################



###################
#### BBMap analysis
echo "Starting ribosomal and mitochondrial filtering with BBMap ..." >> log.out

cp /projects/ncrrbt_share_la/dev_pipe2/human_ribosomal.fa .
mkdir projects

var=(`ls *R1*.fastq.gz`)

	for i in ${var[@]}
	do
	read2=`echo ${i} | sed 's/R1/R2/g'`
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

cp -r $CONTAINER/human_star_index .

var=(`ls *R1*.fastq.gz`)

	for i in ${var[@]}
	do
	read2=`echo ${i} | sed 's/R1/R2/g'`
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
mv star.IIT* star_results

AFTER=`date`
echo "STAR finished on ${AFTER}" >> log.out



#############################
#### featureCounts after STAR
echo "Starting counting with featureCounts ..." >> log.out
cd star_results
files=`ls -d *bam | xargs -n1000`

wget 

apptainer exec $CONTAINER/featurecounts.sif /bin/bash -c \
"featureCounts -B -C -s 2 -p --countReadPairs -T $CPUS -t exon -g gene_id --extraAttributes gene_name,gene_type \
-a /root/gtf/gencode.vM32.annotation.gtf \
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

wget https://sourceforge.net/projects/rseqc/files/BED/Mouse_Mus_musculus/GRCm39_GENCODE_VM27.bed.gz
gunzip GRCm39_GENCODE_VM27.bed.gz

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
	"infer_experiment.py -r GRCm39_GENCODE_VM27.bed -i $i 1> rseqc.$prefix.infer_experiment.txt" || { 
   	echo "RSeQC has an error. Pipeline terminated" >> log.out
    	exit 1
	}
	done

cd ..

AFTER=`date`
echo "RSeQC and sambamba finished on ${AFTER}" >> log.out

fi
####