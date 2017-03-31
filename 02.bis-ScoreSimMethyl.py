#!/usr/bin/python
"""
Analyze the Simulated Genome Methylation Profiles.

Compare the algorithm metric scores from the pipe to the defined levels
of 5'mC methylation in the simulated gDNA sample with differential methylation.

Input Files:
	1.  "001-RefGen/141020d-MetCountTable.txt"
		This file contains two columns. The first is the C<C>GG position in the
		ref genome sequence. The second is the differential methylation level as
		a percentage across all DNA copies in the sample.
	2. "002-Bismark/Y_CRBN37k-001/001e-BIS-SeqReadData-51.fa_bismark_bt2.bismark.cov"
		Second column contains the C<C>GG position in the ref genome sequence.
		The score metrics are in the 10th and 11th columns (0 indexing).  
		
YX:2015
"""
import os
import re
import sys
import gzip
from collections import defaultdict as Ddict

# - - - - - - - - - - -   USER RUN TIME VARIABLES - - - - - - - - - - - - - -
TAG = "R2MB"
ID = '004'

# - - - - - - - - - - -   GLOBAL VARIABLES - - - - - - - - - - - - - -
headers = "Pos\tExpMet\tObsMet\tMETbis\tUMTbis\n"

if len(sys.argv) > 1:
	TAG  = sys.argv[1]
	ID  = sys.argv[2]
	ReadGenFile = int(sys.argv[3])
	readLen=sys.argv[4]

if ReadGenFile == 0:
	REFtable   = ID + "d-MetCountDataTable.txt"
else: 
	REFtable   = ID+"f-BIS-MetCountDataTable.txt"

bismarktable = ID + "e-BIS-SeqReadData-%s.fa_bismark_bt2.bismark.cov.gz" % (readLen)
OUTfile    = "SimMethylScoreTable-01-CpG-"+TAG+".txt"
OUTlost    = "SimMethylScoreTable-02-Lost-"+TAG+".txt"
OUTother   = "SimMethylScoreTable-03-Other-"+TAG+".txt"

results = "002-Bismark/%s-%s/" % (TAG,ID)
refdir  = "001-RefGen/%s-%s/" % (TAG,ID)

def GeneID(ID):
	n = len(ID)
	for i in range (1,(7-n)):
		ID = "0"+ID
	return ID

RefMet   = Ddict(lambda: -1.0)
RefSeq   = Ddict(lambda: 'nada')
FoundMet  = Ddict(lambda: -1.0)
ScoreMet = Ddict(lambda: Ddict(lambda: 'na'))
ConCov   = Ddict(lambda: Ddict(lambda: 'na'))     # Initialize with empty list
CovPos   = []

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - -   M A I N  - - - - - - - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

print "\n\nMatching %Methylation to Metric Scores . . . . . "

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 1. Load the RefMet data . . . ..
print "    Loading Reference MET CpG Data . . . . . "
IN=open(refdir+REFtable, 'r')
FILE=IN.readlines()
IN.close()
for line in FILE:
	x = line.rstrip().split('\t')
	RefMet[x[0]] = float(x[1])
	# RefSeq[x[0]] = x[3]
FILE=[]

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 2. Match OBS and EXP %MET Data . . . . . . .
print "    Matching OBS and EXP %MET Data . . . .  "
IN=gzip.open(results + bismarktable, 'rb')
OUT1=open(results + OUTfile, "w")
OUT1.write(headers)
OUT2=open(results + OUTother, "w")
OUT2.write(headers)
found = 0
while(1):
	line = IN.readline().rstrip()
	if not line: break
	# chr1|GenProSim|004-R2MB|2MB	958	959	57.0469798657718	85	64
	# line format: 0=header; 1=5mCpos; 2=1indexpos; 3=%met; 4=metcount; 5=unmetcount
	x = line.split("\t")
	transpos = str(int(x[1]) - 1)   # << 5mC position indexing is +1 indexed in bismark cov 
	if RefMet[transpos] > 0:
		met = float(x[3])
		# headers = "Pos\tExpMet\tObsMet\tMETscore\tUMTscore\tSeq\n"
		OUT1.write("%s\t%0.1f\t%0.1f\t%s\t%s\n" % (transpos, RefMet[transpos], met, x[4], x[5]))
		# OUT1.write("%s\t%s\t%0.1f\t%s\t%s\t%s\n" % (transpos, RefMet[transpos], met, x[4], x[5], RefSeq[transpos]) )
		FoundMet[transpos] = x[3]
		found += 1
	else:
	# Option for dumping all other 5mC calls . . . . . . 
		met = float(x[3])
		# headers = "Pos\tExpMet\tObsMet\tMETscore\tUMTscore\tSeq\n"
		OUT2.write("%s\t%0.1f\t%0.1f\t%s\t%s\tnada\n" % (transpos, RefMet[transpos], met, x[4], x[5]) )
		
IN.close()
OUT1.close()
OUT2.close()

OUT3=open(results + OUTlost,'w')
OUT3.write("POS\tpMET\tLOST\n")
lost = 0
for (pos, count) in RefMet.iteritems():
	if count > 0:
		if FoundMet[pos] < 0:
			OUT3.write("%s\t%s\t0\n" % (pos, RefMet[pos]))
			lost += 1

OUT3.close()
print "    FOUND = ", found," ;  LOST = ", lost
print "\n\n\n * * * * *   D O N E   * * * * * * \n\n\n"

# EOF ------------------------------------------------------------------------
# AGM2010;2011;2012;2013;2014 . .  .   .     .      .       .        .         .          .

	
	



