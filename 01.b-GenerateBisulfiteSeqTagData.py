#!/usr/bin/python
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
'''
Objective: Generate a model NGS seq read file from a model genome with defined
	cytosine methylation positions that has been treated to RRBS - Reduced Representation
	Bisulfite Sequencing.  The resulting seq files can then be used to
	assess algorithm accuracy.
	
	The key is generating differentially methylated CCGG sites among the frag
	copies such that the final DNA population sampled has %MET at CpG(i) ranging
	from 0% to 100%.
	
	The pool of seq reads generated is dependent upon an enzyme digest and a random
	shearing function to open ends for library preparation.
	
	Seq read output file is structured as a fasta file that has already been QC
	filtered so it has the ".qcfasta".
	
	Based off the RE SeqSim script from 2014. 
	
GenPro/AGM-2014
Modified 2015 YX
'''
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
import os
import re
import sys
import random     # use choice and randint
import string     # use maketrans and translate
from collections import defaultdict as Ddict


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - -   U S E R    V A R I A B L E S  - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TAG       	   	="Yhs_chr1_GRCh38"                     # unique ID str for folders/files
ID 	  			="001"                                 # file prefix for simulated genome output
gMBsize   		=252                        		   # genome size (MB)
genCopyNum		=1                                     # number of genome copies (each with different MET patterns)
seqCycles 		=5                     				   # sequencing cycles (depth of frag sampling to generateseqreads)
readLen 		=51                      			   # sequence read length for each tag
runfolder		="001-RefGen/"   					   # use if working from project folder with the 00.0-SimQuantPipe.sh script
ReadGenFile		=1				    				   # 0 = generate DNA seq; 1 = read DNA seq from an input file
fracBIS			=99                      			   # percent efficiency of bisulfite conversion chemistry
fastQ			=0                         			   # 0=fasta; 1 = fastQ format with phred scores
coinToss 		=10                    				   # % probability of a frag in the sample population being sequenced
shearToss 		=10                   				   # % probability of a sheared frag being sequenced
LoadGenome 		="/Users/user/Thesis/genome/hs_ref_GRCh38_chr1.fa"  #location for the loaded genome if specified

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - -  G L O B A L   V A R I A B L E S  - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Set the NT pool frequency from which random characters will be drawn
#        HUMAN promoter domain composition:
#        A = 0.247; C = 0.251; G = 0.254; T = 0.248
freq = [247, 251, 254, 248]
pNT = ''
for nt in ["A", "C", "G", "T"]:
	fq   = freq.pop(0)
	pNT += ("%s " % nt) * fq

if len(sys.argv) > 1:
	# for ag in sys.argv: print ag
	ID          = sys.argv[1]
	TAG         = sys.argv[2]
	gMBsize     = float(sys.argv[3])            
	genCopyNum  = int(sys.argv[4]) 
	seqCycles   = int(sys.argv[5]) 
	readLen     = int(sys.argv[6])  
	runfolder   = sys.argv[7] 
	ReadGenFile = int(sys.argv[8])
	fracBIS     = int(sys.argv[9])
	fastQ       = int(sys.argv[10])
	coinToss    = int(sys.argv[11]) 
	shearToss   = int(sys.argv[12])
	LoadGenome  = sys.argv[13]
	
runfolder += TAG+"-"+ID+"/"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Structure of genome . . . . . . . . . 
genFragLen = int(gMBsize * 10**6)   # number NT positions in each large frag
CCGG = []
cmet = []
# . . . . . . . . . . . . . 
# Establish distribution of differential methylation values:
#methylstates = [0 for i in range(11)] + [i*5 for i in range(21)]  # Possible set of fractional (%) CpG site methylation states
methylstates = [0 for i in range(5)] + [i*5 for i in range(20)] + [100 for i in range(5)]  # Possible set of fractional (%) CpG site methylation states
ATGC = [ "A", "A", "A", "T", "T", "T", "G", "G", "G", "G", "G", "G", "C" ]  # Limit random CCGG freq by controlling p(C) at second position
ATCG = [ "A", "A", "A", "T", "T", "T", "C", "C", "C", "C", "C", "C", "G" ]  # Limit random CCGG freq by controlling p(G) at fourth position
phred = "456789:;<=>?@ABCDEFGHIJ"
phredz = list(phred)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Defined Files - - - - - - - - - - - - - - - - - - - -
if ReadGenFile == 0:
	GenFile   =  runfolder + "%sb-GenomeSequence.fa" % ID
else:
	GenFile   =  LoadGenome 
#preset methylation rate for single nucleotide in reference genome - - -
refMetTable =  runfolder + "%sd-MetCountDataTable.txt" % ID

if fastQ == 1:
	datafile    =  runfolder + "%se-BIS-SeqReadData-%s.fq" % (ID,str(readLen))
else:
	datafile    =  runfolder + "%se-BIS-SeqReadData-%s.fa" % (ID,str(readLen))
#table of actual methylation rate in synthesised genome or reference genome 	
bisMetTable =  runfolder + "%sf-BIS-MetCountDataTable.txt" % ID

runlog =  runfolder + "%sg-runlog.txt" % ID
WIPE = open(datafile,"w")
WIPE.close()
WIPE = open(bisMetTable,"w")
WIPE.close()
LOG = open(runlog,'w')
LOG.close()
NT = pNT.split()

Hpxy      = string.maketrans('xy','CG')
RC        = string.maketrans("ACGTx", "TGCAy")
Mspxy     = string.maketrans('xy','CC')

refMET   = Ddict( lambda: 0 )
CCGGall  = Ddict( lambda: 0 )
METtable = Ddict( lambda: 0 )

def DUMP(mssg):
	LOG = open(runlog,'a')
	LOG.write(mssg)
	LOG.close()
	print mssg

def GeneID(id):
	n = len(id)
	for i in range (1,(11-n)):
		id = "0"+id
	return id

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#
#      '7MMM.     ,MMF'          db          `7MMF'    `7MN.   `7MF'    
#        MMMb    dPMM           ;MM:           MM        MMN.    M      
#        M YM   ,M MM          ,V^MM.          MM        M YMb   M      
#        M  Mb  M' MM         ,M  `MM          MM        M  `MN. M      
#        M  YM.P'  MM         AbmmmqMA         MM        M   `MM.M      
#        M  `YM'   MM        A'     VML        MM        M     YMM      
#      .JML. `'  .JMML.    .AMA.   .AMMA.    .JMML.    .JML.    YM
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

DUMP("\n\nRunning Genome Seq Tag Construction Set: ID = %s\n\n" % (ID))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# . . . . . Generate the DNA sequence . . . . . . .
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fracGC = 0
if (ReadGenFile == 0):
	fragSeq = 'AAAAAAAAAA'
	lastNT = 'A'
	while(len(fragSeq) < genFragLen-10): 			#genFraglen is the length of whole genome: size times 1000000
		ccgg = random.randint(0,100000)
		if ccgg > 99999:
			fragSeq += 'CCGG'
			fracGC += 4                  			#fracGC is the counter for number of G+C in genome seq
			lastNT = 'G'
		else:
			if lastNT == 'C':
				# Need to reduce random CCGG freq in model sequence; reduce p(C) in second pos
				atgc = random.choice(ATGC)
			elif lastNT == 'G':
				atgc = random.choice(ATCG)
			else:
				# Select from random NT distributions defined at start of script 
				atgc = random.choice(NT)      
			fragSeq += atgc
			if atgc == 'C' or atgc == 'G':
				fracGC += 1
			lastNT = atgc
			
	fragSeq += 'TTTTTTTTTT'
	DUMP("\n\nGenome Seq of %.2f MB.\n" % gMBsize)
	OUT = open(GenFile,'a')
	OUT.write("> Genome Model Sequence %s: %.2f MB\n%s\n" % (ID,gMBsize,fragSeq))
	OUT.close()
else:
	# Load starting sequence to work with . . . . . . 
	
	IN = open(GenFile, "r")
	FILE=IN.readlines()
	IN.close()
	FILE.pop(0)
	for i in xrange(len(FILE)):
		FILE[i] = FILE[i].rstrip()
	fragSeq = ''.join(FILE)
	fracGC = fragSeq.count('C') + fragSeq.count('G')
	genFragLen = len(fragSeq)
	
	# Generate refMettable to work with. . . . . . . .
	index = 0
	while(index > -1):
		index = fragSeq.find('CCGG',index + 3)
		if index > -1 and index < genFragLen-10:
			METtable[index] = random.choice(methylstates)
		else:
			break

	Mettable=open(refMetTable, 'w')
	for pos, pcnt in METtable.iteritems():
		Mettable.write("%s\t%0.2f\n" % (pos+1,pcnt))
	Mettable.close()



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# . . . Process differentially methylated copies . . . .
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if (ReadGenFile==0):
	# Assign target quantitative state to each CCGG site in the sample . . . .
	fracGC = 100*float(fracGC)/float(genFragLen)
	countCCGG = 0
	index = 7
	while(index > -1):
		index = fragSeq.find('CCGG', index+3)
		if index > -1 and index < genFragLen-10:
			CCGGall[index] = random.choice(methylstates)
			countCCGG += 1
	DUMP("Genome Size = %d MB\nGenome G+C = %0.1f%%\nGenome CCGG = %d\n" % (gMBsize,fracGC,countCCGG))


MET   = 0   # uncut
UMT   = 0   # cut
library = 0      # counter for fake gene id numbers in headers
for f in xrange(1,genCopyNum+1):
	print "Genome Copy:", f, " . . . "
	metSeq = list(fragSeq)
	total = 0
	
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	# - - - - - - - - - Methylate - - - - - - - - - - - - -  - -
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	# The list "methylstates" sets probability value for the degree a site will be methylated
	# Use a random selection from 0-100 to translate that probability to a real distribution
	if (ReadGenFile==0):
		for ccggPOS,fracMET in CCGGall.iteritems():
			ccgg = random.randint(0,100)
			if (fracMET > ccgg):
				metSeq[ccggPOS+1] = 'x'
				MET += 1
				refMET[str(ccggPOS+1)] += 1
				# print "\n\n>>> %s <<< \n\n" % metSeq[ccggPOS-14:ccggPOS+10]
			else:
				refMET[str(ccggPOS+1)] += 0    # !! need to establish position key value in dictionary
				UMT += 1
	else:
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	# Use Defined Methylation percentages for each site defined in the Methylation Data Table
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
		
		for ccggPOS,fracMET in METtable.iteritems():
			ccgg = random.randint(0,100)
			if (fracMET > ccgg):
				metSeq[ccggPOS+1] = 'x'
				MET += 1
				refMET[str(ccggPOS+1)] += 1

				# print "\n\npos= %d >>> %s <<< \n\n" % (ccggPOS,metSeq[ccggPOS-14:ccggPOS+10])
			else:
				refMET[str(ccggPOS+1)] += 0    # !! need to establish position key value in dictionary
				UMT += 1

	metSeq = "".join(metSeq)

	# Reverse Complement of the sequence . . . . . . . . . .
	rcseq= metSeq[::-1]
	rcseq = rcseq.translate(RC)    # Cm sites = y in rcseq
	# For bisulfite Msp conversion . . . . need to accurately represent cut-site
	#    The Hpa2 profiling, it was OK to use CxCC/CCyG as teh methyl markers to block RE cuts.
	#    But with Msp/bisulfite, need to shift the methyl mark to the interior cyctosine on the RC copy:
	#				CxGG
	#				GGxC
	# "CxGG" reverse "GGxC" complement "CCyG" shift 5mC marker "CxGG"
	rcseq = re.sub(r'CCyG', 'CxGG', rcseq)
	

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	# - - - - - - - - - Msp I  D i g e s t - - - - - - - - - -
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	Msp1 = []
	recogseq = [ 'CCGG', 'CxGG']    # not methyl-sensitive
	frag = 1
	position = [m.start(0) for m in re.finditer(r'(CCGG)|(CxGG)', metSeq)]
	start = 0
	for i in position:
		fragment = metSeq[start:i+1] + 'CG'
		start = i +1 
		if len(fragment) > readLen:
			Msp1.append(fragment)
	if len(metSeq[position[-1]+1:]) > readLen:
		 Msp1.append(metSeq[position[-1]:])
	# Repeat for RC strand . . . . . . . . . 
	position = [m.start(0) for m in re.finditer(r'(CCGG)|(CxGG)', rcseq)]
	start = 0
	for i in position:
		fragment = rcseq[start:i+1] + 'CG'
		start = i + 1
		if len(fragment) > readLen:
			Msp1.append(fragment)
	if len(rcseq[position[-1]:]) > readLen:
		 Msp1.append(rcseq[position[-1]:])
	frag = 2*(len(position)+1)
	DUMP( "There are %d frags in Msp1.\n" % frag )
	
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	# - - - Sampling Loop to provide coverage depth for Msp cuts 
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	for i in xrange(seqCycles):
		# Shear Msp frags - - - - - - -
		
		Shear = []
		for msp in Msp1:
			if 150 < len(msp) < 225:
				Shear.append(msp)
		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		# Build SeqTags - - - - - - - - -
		OUT = open(datafile,'a')
		for seqrip in Shear:
			odds = random.randint(0,101)
			if (odds <= coinToss and len(seqrip) >= readLen ):
				index = seqrip.find('C',0)
				while(index > -1):
					bisrex = random.randint(0,100)
					if (fracBIS > bisrex):
						# Bisulfite conversion of 5'-OH-Cytosine to Uracil --> Thymidine
						# base x base conversion using reaction efficiency specified at runtime
						seqrip = seqrip[:index] + 'T' + seqrip[index+1:]
					index = seqrip.find('C',index+1)
				seqtrans = seqrip.translate(Mspxy)       # << replace x,y : C,C			
				library += 1
				nid = GeneID("%d" % library)
				fasta = ''
				if fastQ == 1:
					fasta = "@GenProSimSeq-RRBS:%s-%s-%s\n%s\n+\n" % (TAG,ID,nid,seqtrans[0:readLen])
					for i in xrange(0,readLen):
						fasta += "%s" % (random.choice(phredz))
					fasta += "\n"
				else:
					fasta = ">GenProSimSeq-RRBS:%s-%s-%s\n%s\n" % (TAG,ID,nid,seqtrans[0:readLen])
				OUT.write(fasta)
				# print fasta
				# fasta = "> faketag-%s\n%s\n" % (nid,seqtrans[0:readLen])
				# OUT.write(fasta)
		OUT.close()
		
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	# - - - - Sampling loop to provide coverage for random shear sites
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	if shearToss > 0:
		interval = 2
		shift = 10
		for i in xrange(seqCycles):
			
			Shear = []
			# POSTIVE STRAND . . . . . . . 
			# add the ends and tail . . . . . . . .
			lead = re.sub('C', 'T', metSeq[0:readLen+1])
			tail = re.sub('C', 'T', metSeq[-readLen:])
			Shear.append(lead.translate(Mspxy))
			Shear.append(tail.translate(Mspxy))
			# Start shearing . . . . ..  ..  .. . . . 
			L = random.randint(shift,interval*readLen)
			while (L < genFragLen - readLen - 2):
				stop  = L + readLen + 1
				shear = metSeq[L:stop]
				index = shear.find('C',0)
				while(index > -1):
					bisrex = random.randint(0,100)
					if (fracBIS > bisrex):
						# Bisulfite conversion of 5'-OH-Cytosine to Uracil --> Thymidine
						# base x base conversion using reaction efficiency specified at runtime
						shear = shear[:index] + 'T' + shear[index+1:]
					index = shear.find('C',index+1) 
				Shear.append(shear.translate(Mspxy))
				L = L + random.randint(shift,interval*readLen)
		
			# NEGATIVE STRAND . . . . . . .
			# add the ends and tail . . . . . . . . 
			Shear.append(rcseq[0:readLen+1].translate(Hpxy))
			Shear.append(rcseq[-readLen:].translate(Hpxy))
			# Start shearing . . . . ..  ..  .. . . . 
			L = random.randint(shift,interval*readLen)
			while (L < genFragLen - readLen - 2):
				stop  = L + readLen + 1 
				shear = rcseq[L:stop].translate(Hpxy)
				Shear.append(shear)
				L = L + random.randint(shift,interval*readLen)
			# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
			# Build SeqTags - - - - - - - - - - - - - - - - - -
			OUT = open(datafile,'a')
			for s in Shear:
				odds = random.randint(0,101)
				if (odds <= shearToss):
					library += 1
					nid = GeneID("%d" % library)
					fasta = ''
					if fastQ == 1:
						fasta = "@GenProSimSeq-RRBS:%s-%s-%s\n%s\n+\n" % (TAG,ID,nid,s[0:readLen])
						for i in xrange(0,readLen):
							fasta += "%s" % (random.choice(phredz))
						fasta += "\n"
					else:
						fasta = ">GenProSimSeq-RRBS:%s-%s-%s\n%s\n" % (TAG,ID,nid,s[0:readLen])
					OUT.write(fasta)
					# print fasta
			OUT.close()
			

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# OUTPUT Methylation Reference Table . . . . .
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
OUTref=open(bisMetTable, 'w')
for pos, count in refMET.iteritems():
	cov = 100*(float(count)/float(genCopyNum))
	OUTref.write("%s\t%0.2f\n" % (pos,cov))
OUTref.close()


DUMP("\n\n------------------------------------\nThere are %d seq tags in the library." % library)
TOT = MET + UMT
DUMP("\nThere are %d CCGG TOTAL sites.\n" % TOT)
DUMP("There are %d CCGG unmethylated sites.\n" % UMT)
DUMP("There are %d CCGG methylated sites.\n\n" % MET)
DUMP("\n\n * * *   D O N E   * * * \n\n")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  - - - -   E O F  - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


		
