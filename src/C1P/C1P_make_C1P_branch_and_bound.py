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
from babtester import *

import copy
import math

##############################################################
#    C1P_make_C1P_branch_and_bound.py
#
#    Computes a maximal subset of ACS that is C1P using a branch-and-bound
#
##############################################################

if len(sys.argv) < 4:
	print "Use a branch-and-bound to compute a set of ACS that satisfies the C1P"
	print ""
	print "C1P_make_C1P_branch_and_bound.py InputACSFile OutputACSFile RemovedACSFile"
	print '   Input:'
	print '     InputACSFile - the ACS file to make C1P'
	print '   Output:'
	print '     OutputACSFile - the C1P ACS file'
	print '     RemovedACSFile - the ACS removed from the input file to make it C1P'
	
	sys.exit()
#endif

prop = len(sys.argv) > 4 and sys.argv[4] == 'p'		# True if using weights as probs
mat = bm.BinaryMatrix()		# matrix
			
mat.from_file(sys.argv[1])
			
ms = c1p.split_matrix(mat)		# split matrices
matb = bm.BinaryMatrix()		# C1P matrix
mat_rem = bm.BinaryMatrix()		# rows removed
			
j = 1		# iterator for tracing
			
del mat

for m in ms:
	print str(j) + '/' + str(len(ms))
				
	j += 1
			
	# sort matrix to get heuristic as first answer
	m.sort()
			
	# branch and bound
	# optimal rows to remove to make compnonent C1P
	rows = bab.branch_and_bound(m, prop, BABTester(m._height))
		
	rows.sort()
		
	for i in xrange(len(rows) - 1, -1, -1):
		row = m.get_row_info(rows[i])		# row to remove
					
		mat_rem.add_row_info(row)
			
		m.remove_row(rows[i])
	#endfor
	
	# collect usable rows into the C1P matrix
	for r in m._rows:
		matb.add_row_info(r)
	#endfor
	
	print ''
#endfor
	
f = file(sys.argv[2], 'w')
			
f.write(str(matb))
			
f.close()

f = file(sys.argv[3], 'w')
			
f.write(str(mat_rem))
			
f.close()
