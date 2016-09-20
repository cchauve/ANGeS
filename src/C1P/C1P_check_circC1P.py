# ANGES 1.01, reconstructing ANcestral GEnomeS maps
# July 2012.
# Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca

import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))) + '/COMMON')

import bm
import cc1p

##############################################################
#     C1P_check_circC1P.py
#
#     Given an ACS file, ouputs whether it is circC1P or not.
#
##############################################################

if len(sys.argv) <= 1:
	print "Checks if a set of ACS is circular C1P"
	print ''
	print 'Usage:'
	print 'C1P_check_circC1P ACSFile'
	print '   Input:'
	print '     ACSFile - the ACS file to check'
	print '   Output:'
	print '     \'True\' - the matrix is circ-C1P'
	print '     \'False\' - the matrix is not circ-C1P'
	
	sys.exit()
#endif

m = bm.BinaryMatrix()

m.from_file(sys.argv[1])

if cc1p.check_circC1P(m):
	print "True"
	
	sys.exit(0)
else:
	print "False"
	
	sys.exit(1)
#endif
	
