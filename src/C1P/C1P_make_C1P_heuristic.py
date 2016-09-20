# ANGES 1.01, reconstructing ANcestral GEnomeS maps
# July 2012.
# Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca

import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))) + '/COMMON')

import c1p
import bm

##############################################################
#    C1P_make_C1p_heuristic.py
#
#    Computes a C1P subset of ACS using a greedy heuristic
#
##############################################################

if len(sys.argv) < 3:
	print "Use a greedy heuristic to compute a set of ACS that satisfies the C1P"
	print ""
	print "Usage:"
	print "C1P_make_C1P_heuristic InputACSFile OutputACSFile RemovedACSFile"
	print '   Input:'
	print '     InputACSFile - the ACS file to make C1P'
	print '   Output:'
	print '     OutputACSFile - the C1P ACS file'
	print '     RemovedACSFile - the ACS removed from the input file to make it C1P'
	
	sys.exit(0)
#endif

m = bm.BinaryMatrix()		# matrix
mat_rem = bm.BinaryMatrix()	# rows discarded
						
m.from_file(sys.argv[1])

rows = c1p.make_C1P(m)		# rows to remove to make C1P
			
for j in xrange(len(rows) - 1, -1, -1):
	mat_rem.add_row_info(m.get_row_info(rows[j]))
	
	m.remove_row(rows[j])
#endfor
			
f = file(sys.argv[2], 'w')		# C1P matrix file
			
f.write(str(m))
			
f.close()

f = open(sys.argv[3], 'w')

mat_rem.write(f.write)

f.close()
