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
#    C1P_make_mC1P_branch_and_bound.py
#
#    Given a set of weighted ACS with telomeric ACS, computes 
#    the maximum weight subset that is mC1P.
#
##############################################################

def comp_nodes(a, b):
	if len(a._value._support) == len(b._value._support):
		return min(a._value._support) < min(b._value._support)
	#endif
	
	return len(a._value._support) < len(b._value._support)
#endef

if len(sys.argv) < 3:
	print "Use a 1-stage branch-and-bound to compute a set of ACS that satisfies the mC1P"
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
rows_rem = [] 		# rows removed

mat.from_file(sys.argv[1])

ms = c1p.make_intersect_components(mat)		# intersect components
matb = bm.BinaryMatrix()		# mC1P matrix
mat_rem = bm.BinaryMatrix()		# rows removed
			
j = 1		# iterator for tracing

for m in ms:
	print str(j) + '/' + str(len(ms)) + ' ',
	sys.stdout.flush()
				
	j += 1
	
	# sort matrix to get heuristic as first answer
	m.sort()
	
	# branch and bound
	# optimal rows to remove to make compnonent mC1P
	rows = bab.branch_and_bound(m, prop, MC1PTester(m))
		
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
