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
# ACS_compute_supported_adjacencies_nonunique.py
#
# Computes the supported adjacencies of a pair of genomes    
#
#######################################################

if len(sys.argv) < 5:
	print "From orthology blocks computes the supported adjacencies between the genomes of two species, taking into account duplicated markers."
	print "Usage:"
	print "ACS_supported_adjacencies_nonunique.py MarkersFile Species1 Species2 SupportedAdjacenciesFile [0/1]"
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

markers, blocks = extract_blocks(sys.argv[1])

spe1 = sys.argv[2]
spe2 = sys.argv[3]

if sys.argv[5]=='1':
	circ=True
else:
	circ=False
#endif

#filtered = filter_markers(blocks, [spe1,spe2])
#markers = common_markers(filtered,[spe1,spe2])

sigma, seq1 = join_chromosomes(blocks, spe1, markers, False)
first_fake=max(sigma)
sigma, seq2 = join_chromosomes(blocks, spe2, sigma, False)

del markers
del blocks

POS=POS(sigma,seq1)

#print 'pre'

sa = supported_adjacencies_sequences(sigma,seq1,seq2,POS,first_fake,circ)

#print 'done'

for row in sa:
	row._sp = [spe1, spe2]
#endfor

f = open(sys.argv[4], 'w')

sa.write(f.write)

f.close()
