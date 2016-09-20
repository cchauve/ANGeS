# ANGES 1.01, reconstructing ANcestral GEnomeS maps
# July 2012.
# Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca

import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))) + '/COMMON')

import bm
import c1p

##############################################################
#    C1P_check_C1P.py
#
#    Determines whether a set of ACS is C1P or not.
#    Outputs 'True' if the matrix is C1P and 'False' otherwise.
#
##############################################################

if len(sys.argv) < 2:
	print "Prints \'True\' (0) if a matrix is C1P and \'False\' (1) otherwise"
	print ""
	print "Usage:"
	print "C1P_check_C1P.py ACSFile"
	print '   Input:'
	print '     ACSFile -  the file containing the set of ACS to check'
	print '   Output:'
	print '     \'True\' - the matrix is C1P'
	print '     \'False\' - the matrix is not C1P'
	
	sys.exit(0)
#endif

m = bm.BinaryMatrix()		# matrix
m.from_file(sys.argv[1])
		
if c1p.check_C1P(m):
	print 'True'
	sys.exit(0)
else:
	print 'False'
	sys.exit(1)
#endif

