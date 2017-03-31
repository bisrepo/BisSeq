#!/usr/bin/env bash
# ------------------------------------------------------------------
# Pipeline for Analysis of SIMULATED methyl-DNA seq reads.
# ------------------------------------------------------------------
# AGM-2016
# ------------------------------------------------------------------

# When ReadGenFile is "1" with the goal of processing a standard genome sequence
#   for comparison to a GenPro seqread set, then the runtime variables need to be
#   matched to the comparative set:
#       genCopyNum . . . . epigenetic sequence complexity in sample pool
#       seqCycles  . . . . max potential read representation
#       coinToss . . . . . probability of a frag being sequenced

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#                R U N T I M E   V A R S
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TAG="Adam_Feb_Trial"                             # unique ID str for folders/files
SeqID="001"                                      # file prefix for simulated genome output
gMBsize=252                                      # genome size (MB)
genCopyNum=500                                   # number of genome copies (each with different MET patterns)
seqCycles=5                                      # sequencing cycles (depth of frag sampling to generate seq reads)
readLen=40                                       # sequence read length for each tag
simRefGenFolder="001-RefGen/"                    # use if working from project folder with the 00.0-SimQuantPipe.sh script
ReadGenFile=1                                    # 0 = generate DNA seq; 1 = read DNA seq from an input file
fracBIS=99                                       # percent efficiency of bisulfite conversion chemistry
fastQ=0                                          # 0=fasta; 1 = fastQ format with phred scores
coinToss=2                                       # % probability of a frag in the sample population being sequenced
shearToss=0                                      # % probability of a sheared frag being sequenced
LoadGenome="/your_folder/hs_ref_GRCh38_chr1.fa"  #location for the loaded genome if specified




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#               MAKE FOLDER IF MISSING
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if [ ! -d "001-RefGen/$TAG-$SeqID/01-Genome" ]; then
    mkdir -p 001-RefGen/$TAG-$SeqID/01-Genome
fi


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#                F L O W   C O N T R O L
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 0 = NO, do not execute this step . . .
# 1 = YES, execute this step . . . . . 
STEPALL=1        # Execute all steps from 0-6
STEP0=0          # Generate Simulated SeqRead data files
STEP1=0          # File/folder prep
STEP2=0          # BISMARK genome prep
STEP3=0          # BISMARK
STEP4=0          # BISMARK methylation extractor 
STEP5=0          # Prep data tables for R
STEP6=0          # Run R script for OBS vs EXP analysis

if [ $STEPALL = "1" ]; then
    STEP0=1          # Generate Simulated SeqRead data files
    STEP1=1          # File/folder prep
    STEP2=1          # BISMARK genome prep
    STEP3=1          # BISMARK
    STEP4=1          # BISMARK methylation extractor 
    STEP5=1          # Prep data tables for R
    STEP6=1          # Run R script for OBS vs EXP analysis
fi

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
echo "             ,,                                                        "
echo "\`7MM\"\"\"Yp,   db           .M\"\"\"bgd                              "                     
echo "  MM    Yb               ,MI    \"Y                                    "                    
echo "  MM    dP \`7MM  ,pP\"Ybd \`MMb.      .gP\"Ya    ,dW\"Yvd             "  
echo "  MM\"\"\"bg.   MM  8I   \`\"   \`YMMNq. ,M\`   Yb  ,W\`   MM          "  
echo "  MM    \`Y   MM  \`YMMMa. .     \`MM 8M\"\"\"\"\"\"  8M    MM         "
echo "  MM    ,9   MM  L.   I8 Mb     dM YM.    ,  YA.   MM                  "  
echo ".JMMmmmd9  .JMML.M9mmmP\` P\"Ybmmd\"  \`Mbmmd\`   \`Mbm:dMM            "  
echo "                                                   MM                  "  
echo "                                                 .JMML.                "
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
echo "           Simulated Seq Read Processing and Analysis Pipeline"
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
echo "GenPro/2016"


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 0. Run the Simulation Seq Read script . . . .
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if [ $STEP0 = "1" ]; then
    echo "- - - - - - - - - - - - - - - - - - - - - - - - - -"
    echo "0. Running simulated seq read script to generate .qcfasta file . . . "
    if [ ! -d $simRefGenFolder/$TAG-$SeqID ]; then
        mkdir -p $simRefGenFolder/$TAG-$SeqID
    fi
    head -n 55 00.BIS-SimQuantPipe.sh > $simRefGenFolder/$TAG-$SeqID/00-RunConfig.txt
    python 01.b-GenerateBisulfiteSeqTagData.py $SeqID $TAG $gMBsize $genCopyNum $seqCycles $readLen $simRefGenFolder $ReadGenFile $fracBIS $fastQ $coinToss $shearToss $LoadGenome
fi



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 1. Setup folder and file locations . . .
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if [ $STEP1 = "1" ]; then
    echo "- - - - - - - - - - - - - - - - - - - - - - - - - -"
    echo "1. Setup files/folder Genome File . . . "
    if [ ! -d $simRefGenFolder/$TAG-$SeqID/01-Genome ]; then
        mkdir -p $simRefGenFolder/$TAG-$SeqID/01-Genome
    fi
    if [ $ReadGenFile = "1" ]; then
        cp $LoadGenome 001-RefGen/$TAG-$SeqID/01-Genome/
    else    
        cp 001-RefGen/$TAG-$SeqID/${SeqID}b-GenomeSequence.fa 001-RefGen/$TAG-$SeqID/01-Genome/ 
    fi
fi



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 2. Run Bismark Genome_Preparation Script . . .
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if [ $STEP2 = "1" ]; then
    echo "- - - - - - - - - - - - - - - - - - - - - - - - - -"
    echo "2. Prepare Genome File . . . "
    echo $(date +%c) > bismark_time.txt
        /home/user/tools/bismark_v0.14.5/bismark_genome_preparation \
        --bowtie2 \
        --path_to_bowtie /your_path_to_bowtie/bowtie2-2.2.6/ \
        001-RefGen/$TAG-$SeqID/01-Genome/ > bismark_log.txt
fi


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 3. Run Bismark and map reads to reference target sequence . . . . .
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if [ $STEP3 = "1" ]; then
    echo "- - - - - - - - - - - - - - - - - - - - - - - - - -"
    echo "3. Running BISMARK . . . . "

    if [ ! -d 002-Bismark/$TAG-$SeqID/ ]; then
        mkdir -p 002-Bismark/$TAG-$SeqID/
    fi
    /home/user/tools/bismark_v0.14.5/bismark \
    --fasta \
    --sam \
    --output_dir 002-Bismark/$TAG-$SeqID/ \
    --bowtie2 \
    --path_to_bowtie /your_path_to_bowtie/bowtie2-2.2.6/ \
    001-RefGen/$TAG-$SeqID/01-Genome/ \
    001-RefGen/$TAG-$SeqID/${SeqID}e-BIS-SeqReadData-${readLen}.fa >> bismark_log.txt
    echo $(date +%c) >> bismark_time.txt
fi


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 4. Run Bismark methylation extractor to generate reports . . . . .
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if [ $STEP4 = "1" ]; then
    echo "- - - - - - - - - - - - - - - - - - - - - - - - - -"
    echo "4. Running BISMARK Methylation Extractor . . . . "
    /home/user/tools/bismark_v0.14.5/bismark_methylation_extractor \
    --single-end \
    --output 002-Bismark/$TAG-$SeqID/ \
    --bedGraph \
    --zero_based \
    002-Bismark/$TAG-$SeqID/${SeqID}e-BIS-SeqReadData-${readLen}.fa_bismark_bt2.sam >> bismark_log.txt
fi


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 5. Merge the results into an R table for analysis . . . .
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if [ $STEP5 = "1" ]; then
    echo "- - - - - - - - - - - - - - - - - - - - - - - - - -"
    echo "5. Merge all data sources into R table . . . . . "
    python 02.bis-ScoreSimMethyl.py $TAG $SeqID $ReadGenFile $readLen >> bismark_log.txt
fi


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 6. Generate R plots . . . .
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if [ $STEP6 = "1" ]; then
    echo "- - - - - - - - - - - - - - - - - - - - - - - - - -"
    echo "8. Generate R plots . . . . . "
    cp 03.bis-SimMetPlots.R 002-Bismark/$TAG-$SeqID/
    cd 002-Bismark/$TAG-$SeqID/
    R --vanilla <03.bis-SimMetPlots.R --args $TAG $SeqID  >> bismark_log.txt
    cd ../..
    
fi

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
echo ''
echo '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
echo '* * * * *   J O B   R U N   D O N E   * * * * * * * * '
echo '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
echo ''
echo ''

## EOF------------------------------------------


