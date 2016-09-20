# ANGES 1.01, reconstructing ANcestral GEnomeS maps
# July 2012.
# Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca

import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))) + '/COMMON')

import bm
import c1p
import pqtree

########################################################
#    C1P_computre_PQRTree.py
#
#    Computes the PQR-tree of a given ACS file.
#
########################################################

if len(sys.argv) < 4:
	print 'Computes the PQR-tree of a binary matrix.'
	print ''
	print 'Usage:'
	print 'C1P_compute_PQRtree.py ACSFile PQRTreeFile SpeciesName'
	print '   Input:'
	print '     ACSFile - the ACS file'
	print '     SpeciesName - the species name whose genome the tree represents'
	print '   Output:'
	print '     PQRTreeFile - the PQR-tree of the matrix'
	
	sys.exit()
#endif

m = bm.BinaryMatrix()		# matrix
m.from_file(sys.argv[1])
			
f = file(sys.argv[2], 'w')		# PQR-tree file

f.write(">" + sys.argv[3] + "\n")
pqtree.make_PQR_tree(m).write(f.write)
			
f.close()
