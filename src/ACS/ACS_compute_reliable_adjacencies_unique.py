# ANGES 1.01, reconstructing ANcestral GEnomeS maps
# July 2012.
# Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca

import sys

from acs import *
from acs_unique import *

#######################################################
#    ACS_compute_reliable_adjacencies.py
#
#    Computing the reliable adjacencies of three species
#    from a set of unique markers
#
#######################################################

if len(sys.argv) < 7:
	print "Computes the reliable adjacencies of three species with no duplicated markers."
	print ""
	print "Usage:"
	print "ACS_compute_reliable_adjacencies.py MarkersFile SupportedAdjacenciesFile"
	print "     Species1 Species2 Species3 ReliableAdjacenciesFile [0/1]"
	print "   Input:"
	print '    MarkersFile                  - the markers file'
	print '    Species1, Species2, Species3 - the species triple considered'
	print '    SupportedAdjacenciesFile     - the supported adjacencies of Species1 and Species2'
	print '   Output:'
	print '    ReliableAdjacenciesFile - the reliable adjacencies of the three species'
	print '   Options:'
	print '    0 - linear chromosomes'
	print '    1 - circular chromosomes'  
	
	sys.exit()
#endif

temp, blocks = extract_blocks(sys.argv[1])

#read supported ajacencies file
sa = bm.BinaryMatrix()

sa.from_file(sys.argv[2])

spe1 = sys.argv[3]
spe2 = sys.argv[4]
spe3 = sys.argv[5]
is_circular = sys.argv[7] == '1'

filtered_blocks = filter_markers(blocks, [spe1, spe2, spe3])

ra = compute_reliable_adjacencies(filtered_blocks, sa, spe1, spe3, is_circular)

# output reliable adjacencies
ra_file = open(sys.argv[6], 'w')

ra_file.write(str(ra))

ra_file.close()
