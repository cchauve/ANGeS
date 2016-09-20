# ANGES 1.01, reconstructing ANcestral GEnomeS maps
# July 2012.
# Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca

import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))) + '/COMMON')

import bm
import numpy
import spectral

###############################################################
#    SERIATION_compute_PQTree_dotproduct_correlation_matrix.py
#
#    Computes the PQ-tree of a ternary matrix
#    using the spectral seriation method.
#    Computes correlation matrix as an intermediate
#    step.
###############################################################

if len(sys.argv) < 4:
	print 'Computes the PQ-tree of a ternary matrix using spectral seriation, using parameter alpha for undetermined entries.'
	print ''
	print 'Usage:'
	print 'SERIATION_compute_PQtree_dot_product_correlation_matrix.py MatrixFile alpha PQtree SpeciesName'
	print '   Input:'
	print '     MatrixFile - the ternary matrix'
	print '     alpha      - real valued parameter'
	print '     SpeciesName - the species name whose genome the tree represents'
	print '   Output:'
	print '     PQTreeFile - the PQ-tree of the matrix'
	
	sys.exit()
#endif

m = bm.BinaryMatrix()		# matrix
m.from_file(sys.argv[1])

alpha=float(sys.argv[2])

# Create numpy matrix

indices={}
i=0

# indices is a dictionary using the column names as keys for the matrix indices.
for col in m.get_support():
    indices[col]=i
    i+=1 

# Check if there are columns that only appear as Xs.
for row in m:
	try:
	    row._Xs
	except:
	    row._Xs=[]
	for col in row._Xs:
		if col not in indices.keys():
			indices[col]=i
			i+=1    

A=numpy.zeros([len(indices),len(indices)])

# Compute correlation matrix
for col1 in m.get_support():
    for row in m:
	try:
		row._Xs
	except:
		row._Xs= []
        if col1 in row._set:
            for col2 in row._set:
                A[indices[col1],indices[col2]]+=1
            if row._Xs:
                for col2 in row._Xs:
                    A[indices[col1],indices[col2]]+=alpha
        elif row._Xs and col1 in row._Xs:
            for col2 in row._set:
                A[indices[col1],indices[col2]]+=alpha
            for col2 in row._Xs:
                A[indices[col1],indices[col2]]+=alpha**2
            
# Find PQ tree using seriation
T=spectral.seriation(A,indices,0)

o=open(sys.argv[3],'w')
o.write('>'+sys.argv[4].upper()+'\n')
i=1

for t in T:
    if '_P' in t.printTree().split(' ') or '_Q' in t.printTree().split(' '):
        o.write('#CAR'+str(i)+'\n')
	o.write(t.printTree()+'\n')
        i+=1

o.close()
