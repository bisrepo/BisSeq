#!/usr/bin/python
"""
Analyze the Simulated Genome Methylation Profile using BSseeker2

Compare the algorithm metric scores from the GenPro pipe to the defined levels
of 5'mC methylation in the simulated gDNA sample with differential methylation.

Input Files:
	1.  "001-RefGen/Y_hChr1_Feb-001/001f-BIS-MetCountTable.txt"
		This file contains two columns. The first is the C<C>GG position in the
		ref genome sequence. The second is the differential methylation level as
		a percentage across all DNA copies in the sample.
	2. "003-bsseeker/001e-BIS-SeqReadData-76.fa_rrbsse.bam.CGmap.gz"
	    This compressed tab-delimited file contains 8 column. The third column is
	    the <C> position on waston strand, the 6th column is methylattion level
	    for that position.

YX:2016
"""


import sys
import gzip
from collections import defaultdict as Ddict

# - - - - - - - - - - -   USER RUN TIME VARIABLES - - - - - - - - - - - - - -
TAG = "Y_hChr1_Feb"
ID = '003'

# - - - - - - - - - - -   GLOBAL VARIABLES - - - - - - - - - - - - - -
headers = "Pos\tExpMet\tObsMet\tMETbis\tUMTbis\n"

if len(sys.argv) > 1:
	TAG  = sys.argv[1]
	ID  = sys.argv[2]

REFtable   = ID+"f-BIS-MetCountDataTable.txt"
bsseekerMap = TAG  + "-" + ID + "-rrbs.bam.CGmap.gz" 
OUTfile    = ID + "-BSseekerTable-01-CpG-"+TAG+".txt"
OUTlost    = ID + "-BSseekerTable-02-Lost-"+TAG+".txt"
OUTother   = ID + "-BSseekerTable-03-Other-"+TAG+".txt"

results = "003-bsseeker/"
refdir  = "001-RefGen/%s-%s/" % (TAG,ID)

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
IN=gzip.open(bsseekerMap, 'rb')
OUT1=open(results + OUTfile, "w")
OUT1.write(headers)
OUT2=open(results + OUTother, "w")
OUT2.write(headers)
found = 0
while(1):
	line = IN.readline().rstrip()
	if not line: break
	# chr1    C   3001624 CHG CA  0.0 0   9
	# (1) chromosome
    # (2) nucleotide on Watson (+) strand
    # (3) position
    # (4) context (CG/CHG/CHH)
    # (5) dinucleotide-context (CA/CC/CG/CT)
    # (6) methylation-level = #_of_C / (#_of_C + #_of_T).
    # (7) #_of_C (methylated C, the count of reads showing C here)
    # (8) = #_of_C + #_of_T (all Cytosines, the count of reads showing C or T here)

	x = line.split("\t")
        x = map(int, [x[2], x[6], x[7]] ) # position/No.Cs/No.Cs+Ts
	transpos = str(x[0] - 1)   # << 5mC position indexing it's 1-based that BSseeker2 use
	if RefMet[transpos] > 0:
		met = 100 * (float(x[1]) / x[2])
		# headers = "Pos\tExpMet\tObsMet\tMETscore\tUMTscore\tSeq\n"
		OUT1.write("%s\t%0.1f\t%0.1f\t%s\t%s\n" % (transpos, RefMet[transpos], met, "nada", "nada"))
		# OUT1.write("%s\t%s\t%0.1f\t%s\t%s\t%s\n" % (transpos, RefMet[transpos], met, x[4], x[5], RefSeq[transpos]) )
		FoundMet[transpos] = met 
		found += 1
	else:
	# Option for dumping all other 5mC calls . . . . . .
		met =  100 * (x[1] / x[2])
		# headers = "Pos\tExpMet\tObsMet\tMETscore\tUMTscore\tSeq\n"
		OUT2.write("%s\t%0.1f\t%0.1f\t%s\t%s\tnada\n" % (transpos, RefMet[transpos], met, "nada", "nada") )

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



