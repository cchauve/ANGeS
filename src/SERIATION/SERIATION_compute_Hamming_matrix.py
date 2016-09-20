# ANGES 1.01, reconstructing ANcestral GEnomeS maps
# July 2012.
# Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca

import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))) + '/COMMON')

import bm
import numpy
import hamming_matrix

########################################################
#    SERIATION_compute_Hamming_Matrix.py
#
#    Computes the Hamming distance matrix of a set of ACS
#
########################################################

if len(sys.argv) < 4:
	print 'Computes a parameterized Hamming distance between markers.'
	print ''
	print 'Usage:'
	print 'SERIATION_compute_Hamming_Matrix.py ACSFile Alpha DistancesFile'
	print '   Input:'
	print '     ACSFile - the ACS (possibly with X-markers)'
	print '     Alpha   - the parameter: 0X has distance alpha'
	print '   Format assumption:'
	print '     ACSFile is a ternary 0/1/X matrix'
	print '   Output:'
	print '     DistancesFile - the distances between markers'
	
	sys.exit()
#endif

m = bm.BinaryMatrix()		# matrix
m.from_file(sys.argv[1])

# Read parameter alpha
alpha=float(sys.argv[2])

# Computes distance matrix and stores values in a dictionary
distance_matrix= hamming_matrix.hamming_matrix(m,alpha)

del m

# Write to file
o=open(sys.argv[3],'w')
hamming_matrix.write_to_file(distance_matrix,o)
o.close()
