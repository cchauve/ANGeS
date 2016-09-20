# ANGES 1.01, reconstructing ANcestral GEnomeS maps
# July 2012.
# Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca

import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))) + '/COMMON')

import bm
from acs import *
from acs_nonunique import *

#######################################################

# ACS_compute_common_intervals_sequences.py

#

# From orthology blocks computes the common

# intervals between the genomes of two species, with

# possibly dissimilar gene content.

#######################################################

if len(sys.argv) < 6:
	print "From orthology blocks computes the common intervals between the genomes of two species, taking into account duplicated markers."
	print "Usage:"
	print "ACS_compute_common_intervals_multiple_hits.py MarkersFile Species1 Species2 CommonIntervalsFile [ALL/MAX/STRONG] [0/1]"
	print '   Input:'
	print '    MarkersFile        - the markers file'
	print '    Species1, Species2 - the species pair considered'
	print '   Output:'
	print '    CommonIntervalsFile - the ACS file to add missing markers to'
	print '   Options:'
	print '    ALL    - compute all the common intervals'
	print '    MAX    - compute only the maximal common intervals'
	print '    STRONG - compute only the strong common intervals'
	print '    0      - input genomes are linear'	
	print '    1      - input genomes are circular'	
	sys.exit()
#endif

markers, blocks = extract_blocks(sys.argv[1])

spe1 = sys.argv[2]
spe2 = sys.argv[3]

if sys.argv[5].upper() == 'MAX':
	signal=1
elif sys.argv[5].upper() == 'STRONG':
	signal=2
else:
	signal=0
#endif

if sys.argv[6]=='1':
	circ=True
else:
	circ=False
#endif

#filtered = filter_markers(blocks, [spe1,spe2])
#markers = common_markers(blocks,[spe1,spe2])

sigma, seq1 = join_chromosomes(blocks, spe1, markers, circ)
first_fake=max(sigma)
sigma, seq2 = join_chromosomes(blocks, spe2, sigma, circ)

del markers
del blocks

[POS,NUM]=POSandNUM(sigma,seq1,first_fake)

ACS = common_intervals_sequences(sigma,seq1,seq2,POS,NUM,signal,first_fake,circ)

for row in ACS:
	row._sp = [spe1, spe2]
#endfor

f = open(sys.argv[4], 'w')

ACS.write(f.write)

f.close()
