#!/usr/bin/env bash
#PBS -N Y_bsseeker-002 -S /bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l walltime=100:00:00
#PBS -d /home/xuyubo/Y_hChr1_Feb-002
#PBS -V


#========== Y.X 2016 ==========
#--------------------------------------------------
#pipeline using BSseeker2 to map & analyze simulated reads
#

#--------------------------------------------------
#          runtime variables
#--------------------------------------------------
TAG="Y_hChr1_Feb"                                               #Define the prefix of mapping output file
ID="002"                                                        #ID number of simulating work
readLen=40                                                      #length of read
seqRead=001-RefGen/$TAG-$ID/${ID}e-BIS-SeqReadData-$readLen.fa  #Path of simulated reads file
refGenome="/home/user/genome/hs_ref_GRCh38_chr1.fa"             #Path to reference genome file
conRefGenome="003-bsseeker/converted_genome"                    #Path to converted reference genome
pathToBowtie="/home/user/tools/bowtie2-2.2.6"                   #Path to bowtie2
pathToBsSeeker="/home/user/tools"                               #Prefix of path to BSseeker2


#--------------------------------------------------
#           flow control
#--------------------------------------------------
ALL=1    #Execute all steps
STEP1=0  #Convert reference genome for BSseeker2
STEP2=0  #Map reads
STEP3=0  #Extract methylation from mapping result
STEP4=1  #Extract methylation rate to R table
STEP5=1  #Generate ObsMet VS ExpMet rate figure

if [ $ALL = "1" ]; then
    STEP1=1
    STEP2=1
    STEP3=1
    STEP4=1
    STEP5=1
fi

#----------Clear previous log file and txt file----------
rm -f bsseeker_time.txt
rm -f Y_bsseeker-$i.*
rm -f bsseeker_mapping.txt
rm -f bsseeker_methylation.txt
rm -f bsseeker_prep.txt

#----------Build genome index----------
if [ $STEP1 = "1" ]; then
    echo "- - - - - - - - - - - - - - - - - - - - - - - - - -"
    echo "1. Building index for BSseeker2 "
    SECONDS=0
    if [ ! -d $conRefGenome ]; then
        mkdir -p $conRefGenome 
    fi
    python $pathToBsSeeker/BSseeker2/bs_seeker2-build.py \
        -f $refGenome \
        --aligner=bowtie2 \
        -p $pathToBowtie \
        -d $conRefGenome \
        -r \
        -l 125 \
        -u 250 \
        -c C-CGG
    echo $SECONDS > bsseeker_prep.txt
fi

#----------Mapping reads to genome----------
if [ $STEP2 = "1" ]; then
    echo "- - - - - - - - - - - - - - - - - - - - - - - - - -"
    echo "2. Mapping reads to reference genome"
    python $pathToBsSeeker/BSseeker2/bs_seeker2-align.py \
        -i $seqRead \
        -r \
        -L 125 \
        -U 250 \
        -m 0 \
        -g $refGenome \
        --aligner=bowtie2 \
        -o ${TAG}-${ID}-rrbs.bam \
        -f bam \
        -p $pathToBowtie \
        -d $conRefGenome >& bsseeker_mapping.txt 
    echo $SECONDS > bsseeker_time.txt
fi    

#----------Extract methylation call----------
if [ $STEP3 = "1" ]; then
    echo "- - - - - - - - - - - - - - - - - - - - - - - - - -"
    echo "3. Extracting methylation calls from mapping results"
    python $pathToBsSeeker/BSseeker2/bs_seeker2-call_methylation.py \
        -i ${TAG}-${ID}-rrbs.bam \
        -d $conRefGenome/hs_ref_GRCh38_chr1.fa_rrbs_125_250_bowtie2 >& bsseeker_methylation.txt
fi    
    
#----------Prepare methylation table----------
if [ $STEP4 = "1" ]; then
   echo "- - - - - - - - - - - - - - - - - - - - - - - - - -"
   echo "4. Preparing methylation tables"
   python 05.BSseeker-ScoresSimMethy.py $TAG $ID
fi


#----------Visualize methylation ----------
if [ $STEP5 = "1" ]; then
    echo "- - - - - - - - - - - - - - - - - - - - - - - - - -"
    echo "5. Generating methylation call plots"
    cp 06.BSseeker-SimMetPlots.R 003-bsseeker/
    cd 003-bsseeker/
    R --vanilla <06.BSseeker-SimMetPlots.R --args $TAG $ID
    cd ..
fi


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
echo ''
echo '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
echo '* * * * *   J O B   R U N   D O N E   * * * * * * * * '
echo '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
echo ''
echo ''

## EOF------------------------------------------
