# ANGES 1.01, reconstructing ANcestral GEnomeS maps
# July 2012.
# Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca

import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))) + '/COMMON')

import bm
import mc1p
import c1p

#######################################################
#    C1P_check_mC1P.py
#
#    Determines if a binary matrix is mC1P.
#    Outputs the martrix's minumum multiplicity if it is mC1P and -1 otherwise.
#
#######################################################

def usage():
	print 'Determines whether a a matrix is mC1P and prints its minimum multiplicity'
	print ''
	print 'Usage:'
	print 'C1P_check_mC1P.py MatrixFile'
	print '   Input:'
	print '     MatrixFile - the matrix file to check'
	print '   Output:'
	print '     \'-1\' - if for no multiplicity the matrix is mC1P'
	print '     n >= 0 - the minimum multiplicity for the matrix to be mC1P'
#enddef

# exits the program with True or False
# ret - bool
def exit(ret):
	print ret

	if ret >= 0:		
		sys.exit(0)
	else:		
		sys.exit(1)
	#endif
#enddef

if len(sys.argv) < 2:
	usage()
	
	sys.exit()
#endif

# Steps are from Algorithm 1 in "Tractibility results for the consecutive-ones property with multiplicity" by Chauve, Manuch, Patternson and Wittler as of August 2011.
# modification to Algorithm  to fix error in the Algorithm when the LCA's was a Q node and the row was not at the end or beginning of the LCA's children.

# Steps 1 and 2 and initialization
m = bm.BinaryMatrix()		# matrix
multM = bm.BinaryMatrix()		# telomere matrix

m.from_file(sys.argv[1])

m.sort()

# seperate telomere and non telomere rows
for i in xrange(m._height - 1, -1, -1):
	row = m.get_row_info(i)		# current row
	
	if not row._isT:
		break
	#endif
	
	m.remove_row(i)
	
	ad = True		# True if we will add the telomere row to the list of minimal telomere rows
	
	# ensure list of telomere rows is minimal
	for j in xrange(multM._height - 1, -1, -1):
		tRow = multM.get_row_info(j)		# current row in telomere list
		
		if tRow._set < row._set:
			ad = False
			
			break
		#endif
		
		if row._set < tRow._set:
			multM.remove_row(j)
		#endif
	#endfor
	
	if ad:
		multM.add_row_info(row)
	#endif
#endfor

if not c1p.check_C1P(m):		# check matrix is C1P
	exit(-2)
#endif

# make PQ-tree
tre = pqtree.make_PQR_tree(m)

mc1p.fix_tree_node(tre._head)

# Step 3
for r in multM._rows:
	# Steps a and b and c
	if not mc1p.check_LCA_path(r._set, tre):
		exit(-1)
	#endif
#endfor:

# Step 4
if tre._head._value._type == 'Q':
	exit(tre._head._value._used)
# Step 5
else:
	# Step 5a
	K1, K2 = mc1p.count_used(tre)		# K1, K2
	
	# Step 5b
	K = K1 // 2 + K1 % 2 + K2 		# K
	
	if (K1 == 0 and K2 > 0):
		K += 1
	#endif
	
	#Step 5c
	exit(K)
#endif
