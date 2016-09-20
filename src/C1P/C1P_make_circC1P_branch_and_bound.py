# ANGES 1.01, reconstructing ANcestral GEnomeS maps
# July 2012.
# Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca

import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))) + '/COMMON')

import bm
import c1p
import bab
import sort
import cc1p
from babtester import *

import copy
import math

##############################################################
#     C1P_make_circC1P_branch_nad_bound.py
#
#    Computes a circular C1P maximal subset of ACS using a branch-and-bound
#
##############################################################

if len(sys.argv) < 5:
        print "Use a branch-and-bound to compute a maximal set of ACS that satisfies the circular C1P"
        print ""
	print ""
	print "C1P_make_circC1P_branch_and_bound.py InputACSFile [whole/max] OutputACSFile RemovedACSFile"
	print '   Input:'
	print '     InputACSFile - the ACS file to make circ-C1P'
	print '   Output:'
	print '     OutputACSFile - the circ-C1P ACS file'
	print '     RemovedACSFile - the ACS emoved from the input file to make it circular C1P'
	print '   Options:'
	print "    whole - makes the whole matrix circular circ-C1P"
	print "    max   - makes the maximal connected components of the matrix circ-C1P"
	
	sys.exit()
#endif

mat = bm.BinaryMatrix()		# matrix
			
mat.from_file(sys.argv[1])

mode = sys.argv[2]

if mode == 'whole':
	matb, mat_rem = circC1P_bab(mat)
elif mode == 'max':
	maxs = c1p.make_intersect_components(mat)		# split matrices
	matb = bm.BinaryMatrix()		# C1P matrix
	mat_rem = bm.BinaryMatrix()		# rows removed
	
	j = 1		# iterator for tracing
				
	del mat
	
	for max_comp in maxs:
		print 'Max:' + str(j) + '/' + str(len(maxs)) + ' '
		
		j += 1
		
		circC1P_bab(max_comp, matb, mat_rem)
	#endfor
#endif
	
f = file(sys.argv[3], 'w')
			
f.write(str(matb))
			
f.close()

f = file(sys.argv[4], 'w')
			
f.write(str(mat_rem))
			
f.close()
