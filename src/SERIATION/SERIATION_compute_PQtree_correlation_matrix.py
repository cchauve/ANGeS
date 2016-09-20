# ANGES 1.01, reconstructing ANcestral GEnomeS maps
# July 2012.
# Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca

import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))) + '/COMMON')

import bm
import numpy
import spectral

########################################################
#    SERIATION_compute_PQTree.py
#
#    Computes the PQ-tree of a binary matrix
#    using the spectral seriation method.
#    Computes correlation matrix as an intermediate
#    step.
########################################################

if len(sys.argv) < 4:
	print 'Computes the PQ-tree of a binary matrix using spectral seriation.'
	print ''
	print 'Usage:'
	print 'SERIATION_compute_PQtree_correlation_matrix.py MatrixFile PQtree SpeciesName'
	print '   Input:'
	print '     MatrixFile - the binary matrix'
	print '     SpeciesName - the species name whose genome the tree represents'
	print '   Output:'
	print '     PQTreeFile - the PQ-tree of the matrix'
	
	sys.exit()
#endif

m = bm.BinaryMatrix()		# matrix
m.from_file(sys.argv[1])

# Create numpy matrix

indices={}
i=0

# indices is a dictionary using the column names as keys for the matrix indices.
for col in m.get_support():
    indices[col]=i
    i+=1 

M=numpy.zeros([m._height,len(indices)])

# Populate binary matrix
i=0
for row in m:
    for col in row._set:
	M[i,indices[col]]=1
    i+=1

del m

# Compute correlation matrix as matrix product
A=numpy.dot(M.T,M)
# Find PQ tree using seriation

T=spectral.seriation(A,indices,0)

o=open(sys.argv[2],'w')
o.write('>'+sys.argv[3].upper()+'\n')
i=1

for t in T:
    o.write('#CAR'+str(i)+'\n')
    o.write(t.printTree()+'\n')
    i+=1

o.close()
