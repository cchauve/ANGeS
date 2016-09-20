# ANGES 1.01, reconstructing ANcestral GEnomeS maps
# July 2012.
# Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca

import sys

from acs import *
from acs_unique import *

#######################################################
#    ACS_compute_common_intervals.py
#
#    Computes the common intervals of a pair of genomes
#
#######################################################

if len(sys.argv) < 5:
	print "From orthology blocks computes the (maximal) common intervals between the genomes of two species with no duplicated markers."
	print ""
	print "Usage:"
	print "ACS_compute_common_intervals_permutations.py MarkersFile Species1 Species2 CommonIntervalsFile [STRONG|ALL|MAX] [0/1]"
	print '   Input:'
	print '    MarkersFile - the markers file'
	print '    Species1, Species2 - the species pair considered'
	print '   Output:'
	print '    Commonintervals - the common intervals'
	print '   Options:'
	print "    STRONG - computes the strong common intervals"
	print "    ALL    - computes all the common intervals"
	print "    MAX    - computes only maximal common intervals"
	print "    0      - input genomes are linear"
	print "    1      - input genomes are circular"
	
	sys.exit()
#endif

temp, blocks = extract_blocks(sys.argv[1])

spe1 = sys.argv[2]
spe2 = sys.argv[3]

output_intervals=0
if sys.argv[5].upper() == 'ALL':
	output_intervals=1
elif sys.argv[5].upper() == 'MAX':
	output_intervals=2
#bMax = len(sys.argv) <= 5 or sys.argv[5].upper() != 'ALL'

if sys.argv[6]=='1':
	circ=True
else:
	circ=False
#endif

filtered_blocks = filter_markers(blocks, [spe1, spe2])
ci = compute_common_intervals(filtered_blocks, spe1, spe2, output_intervals, circ)

# output common intervals
ci_file = open(sys.argv[4], 'w')

ci_file.write(str(ci))

ci_file.close()
