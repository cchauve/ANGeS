# ANGES 1.01, reconstructing ANcestral GEnomeS maps
# July 2012.
# Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca

import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))) + '/COMMON')

import bm
from acs import *

#######################################################
#    ACS_add_telomeres.py
#
#    Add conserved telomeric ACS
#
#######################################################

if len(sys.argv) < 5:
	print 'Adds rows with telomeres in a pair of species to an ACS file.'
	print ''
	print "Usage:"
	print 'ACS_add_telomeres.py MarkersFile Species1 Species2 InputACSFile OutputACSFile'
	print '   Input:'
	print '    MarkersFile        - the markers file'
	print '    Species1, Species2 - the species pair considered'
	print '    InputACSFile       - the ACS file to add telomeres to'
	print '   Output:'
	print '    OutputACSFile      - the ACS file with telomere rows'

	sys.exit()	
#endif

temp, ob = extract_blocks(sys.argv[1])		# orthology blocks

acs = bm.BinaryMatrix()		# ancestral consecutive syntenies

acs.from_file(sys.argv[4])

add_telomeres(ob, sys.argv[2], sys.argv[3], acs)

f = file(sys.argv[5], 'w')		# syntenies with telomeres output file

f.write(str(acs))

f.close()
