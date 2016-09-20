# ANGES 1.01, reconstructing ANcestral GEnomeS maps
# July 2012.
# Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca

import sys

from acs import *
from acs_unique import *

#######################################################
#    ACS_compute_supported_adjacencies.py
#
#    Computing supported adjacencies between a pair of species
#
#######################################################

if len(sys.argv) < 5:
	print "Computes the supported adjacencies of two species."
	print ""
	print "Usage:"
	print "ACS_supported_adjacencies_unique.py MarkersFile Species1 Species2 SupportedAdjacenciesFile [0/1]"
	print "   Input:"
	print '    MarkersFile        - the markers file'
	print '    Species1, Species2 - the species pair considered'
	print '   Output:'
	print '    SupportedAdjancenciesFile - the supported adjacencies of the two species'
	print '   Options:'
	print "    0 - input genomes are linear"
	print "    1 - input genomes are circular"
	
	sys.exit()
#endif

if sys.argv[5]=='1':
	circ=True
else:
	circ=False
#endif

temp, blocks = extract_blocks(sys.argv[1])

spe1 = sys.argv[2]
spe2 = sys.argv[3]

filtered_blocks = filter_markers(blocks, [spe1, spe2])

sa = compute_adjacencies(filtered_blocks, spe1, spe2, circ)

# output supported adjacencies
sa_file = open(sys.argv[4], 'w')

sa_file.write(str(sa))

sa_file.close()
