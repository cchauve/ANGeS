# ANGES 1.01, reconstructing ANcestral GEnomeS maps
# July 2012.
# Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca

import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))) + '/COMMON')

import bm
import weight

def usage():
	print "Computes the weight of a set of ACS based on pattern of occurrence in extant species"
	print ""
	print "Usage:"
	print "ACS_weight_linear_interpolation.py InputACSFile SpeciesTreeFile OutputACSFile [d]"
	print "   Input:"
	print "    InputACSFile - the file of ACS to weight"
	print "    SpeciestreeTreeFile - the species tree of the species"
	print "   Output:"
	print "    OutputACSFile - the weighted mnatrix file"
	print "   Options:"
	print "    d - include if the matrix has doubled markers"
#enddef

if len(sys.argv) < 4:
	usage()
	
	sys.exit()
#endif

double = len(sys.argv) > 4 and sys.argv[4] == 'd'		# True if the matrix has
														# doubled markers to marker
														# single marker weights
														# correctly
														
if len(sys.argv) > 4 and sys.argv[4] != 'd':
	offset = int(sys.argv[4])
elif len(sys.argv) > 5:
	offset = int(sys.argv[5])
else:
	offset = 0
#endif
	
m = bm.BinaryMatrix()		# matrix
		
m.from_file(sys.argv[1])
		
weight.weight_matrix(m, sys.argv[2], double, offset)
		
f = file(sys.argv[3], 'w')		# weighted matrix file
		
f.write(str(m))
		
f.close()
