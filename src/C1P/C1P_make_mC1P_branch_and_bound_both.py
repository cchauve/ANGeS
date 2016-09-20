# ANGES 1.01, reconstructing ANcestral GEnomeS maps
# July 2012.
# Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca

import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))) + '/COMMON')

import bm
import mc1p
import bab
import tree
import c1p
import sort
from babtester import *

import copy
import time

##############################################################
#    C1P_make_mC1P_branch_and_bound_both.py
#
#    Given a set of weighted ACS with telomeric ACS, computes 
#    a subset that is mC1P using a 2-stage branch-and-bound.
#
##############################################################

if len(sys.argv) < 3:
	print "Use a 2-stage branch-and-bound to compute a set of ACS that satisfies the mC1P"
	print ""
	print "C1P_make_mC1P_branch_and_bound.py InputACSFile OutputACSFile RemovedACSFile"
	print '   Input:'
	print '     InputACSFile - the ACS file to make mC1P'
	print '   Output:'
	print '     OutputACSFile - the mC1P ACS file'
	print '     RemovedACSFile - the ACS removed from the input file to make it mC1P'
	
	sys.exit()
#endif

prop = len(sys.argv) > 3 and sys.argv[3] == 'p'		# True if using weights as probs
mat = bm.BinaryMatrix()		# matrix

mat.from_file(sys.argv[1])

matb, mat_rem = C1P_and_mC1P_bab(mat)

f = file(sys.argv[2], 'w')

f.write(str(matb))

f.close()

f = file(sys.argv[3], 'w')
			
f.write(str(mat_rem))
			
f.close()
