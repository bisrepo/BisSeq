#!/usr/bin/env bash

#========== Y.X 2016 ==========
#--------------------------------------------------
#pipeline using BSMAP to map & analyze simulated reads
#

#--------------------------------------------------
#          runtime variables
#--------------------------------------------------
TAG="Y_hChr1_Feb"                                           #Define the prefix of mapping output file
ID="002"                                                    #ID number of simulating work
readLen=40                                                  #length of read
seqRead="001-RefGen/"$TAG"-"$ID"/"$ID"e-BIS-SeqReadData-"$readLen".fa"  #Path of simulated reads file
refGenome="/home/user/genome/hs_ref_GRCh38_chr1.fa"       #Path to reference genome file
mapResult="004-bsmap/"$TAG"-"$ID"-rrbs.bam"
methResult="004-bsmap/"$TAG"-"$ID"-methratio.txt"

#--------------------------------------------------
#           flow control
#--------------------------------------------------
ALL=0    #Execute all steps
STEP1=0  #Map reads to reference genome
STEP2=0  #Extract methylation rates
STEP3=1  #Extract methylation rate to R table
STEP4=1  #Generate ObsMet VS ExpMet rate figure

if [ $ALL = "1" ]; then
    STEP1=1
    STEP2=1
    STEP3=1
    STEP4=1
fi


#----------remove previous run logs & generate new log----------
#rm -f bsmap_time.txt
#rm -f bsmap_mapping.txt 
#rm -f bsmap_methylation.txt
#rm -f Y_bsmap_$ID.*

#----------Mapping reads to reference genome----------
if [ $STEP1 = "1" ]; then
    echo "- - - - - - - - - - - - - - - - - - - - "
    echo "01. running BSMAP to map to reads to reference genome"
    SECONDS=0
    /home/user/tools/bsmap-2.90/bsmap \
         -a $seqRead \
         -d $refGenome \
         -o $mapResult \
         -D C-CGG \
         -p 1 \
         -v 0 \
         -n 0 \
         -z 64 >& bsmap_mapping.txt
    echo $SECONDS > bsmap_time.txt
fi

#----------Extracting methylation info----------
if [ $STEP2 = "1" ]; then
    echo "- - - - - - - - - - - - - - - - - - - - "
    echo "02. extracting methylation info from bam file"
    python /home/user/tools/bsmap-2.90/methratio.py \
          -d $refGenome \
          -o $methResult \
          -u $mapResult >& bsmap_methylation.txt

fi

#----------Generate table for R plots----------
if [ $STEP3 = "1" ]; then
    echo "- - - - - - - - - - - - - - - - - - - - "
    echo "03. generating table for R plots"
    python 08.BSMAP-ScoresSimMethy.py $TAG $ID $readLen

fi

#----------Produce methylation rate figure----------
if [ $STEP4 = "1" ]; then
    echo "- - - - - - - - - - - - - - - - - - - - "
    echo "04. Visualize mapping results"
    cp 09.BSMAP-SimMetPlots.R 004-bsmap/
    cd 004-bsmap/
    R --vanilla <09.BSMAP-SimMetPlots.R --args $TAG $ID
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
