# ANGES 1.01, reconstructing ANcestral GEnomeS maps
# July 2012.
# Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca

import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))) + '/COMMON')
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))) + '/ACS')
import bm
from acs import *

#######################################################

# ACS_check_for_overlaps_markers.py

#

# Checks if there are overlapping markers in the 

# orthology blocks. If so, finds overlaps and 

# lists species and markers overlapping. 

#######################################################

if len(sys.argv) < 2:
	print "From orthology blocks, checks for overlapping markers."
	print "Usage:"
	print "ACS_check_for_overlaps_markers.py MarkersFile OverlapsFile"
	print '   Input:'
	print '    MarkersFile - the markers file'
	print '   Output:'
	print '    List of overlapping markers, and the species they occur in, with length of each marker and the overlap'
	sys.exit()
#endif

markers, blocks = extract_blocks(sys.argv[1])

o=open(sys.argv[2],'w')
for species in blocks.keys():
#	print species
	for chrom in blocks[species]:
		i=0
		while i<len(chrom)-1:
			if chrom[i]._st< chrom[i+1]._st and chrom[i]._end>chrom[i+1]._st:
				o.write('Overlap detected in species'+ species+'\n')
				if chrom[i+1]._end>chrom[i]._end:
					o.write(str(chrom[i]._id)+' '+str(abs(chrom[i]._end-chrom[i]._st)+1)+' '+str(chrom[i+1]._id)+' '+str(abs(chrom[i+1]._end-chrom[i+1]._st)+1)+' '+str(abs(chrom[i]._end-chrom[i+1]._st))+'\n')
				elif chrom[i+1]._end<chrom[i]._end:
					o.write(str(chrom[i]._id)+' '+str(abs(chrom[i]._end-chrom[i]._st)+1)+' '+str(chrom[i+1]._id)+' '+str(abs(chrom[i+1]._end-chrom[i+1]._st)+1)+' '+str(abs(chrom[i+1]._end-chrom[i+1]._st)+1)+'\n')
			i+=1
			
o.close()
