# ANGES 1.01, reconstructing ANcestral GEnomeS maps
# July 2012.
# Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca

import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))) + '/COMMON')

import bm
import cc1p

##############################################################
#     C1P_make_circC1P_heuristic.py
#
#    Computes a circular C1P subset of ACS using a greedy heuristic
#
##############################################################

if len(sys.argv) < 5:
        print "Use a greedy heuristic to compute a set of ACS that satisfies the circular C1P"
        print ""
	print ""
	print "C1P_make_circC1P_heuristic.py InputACSFile [whole/max] OutputACSFile RemovedACSFile"
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

m = bm.BinaryMatrix()

m.from_file(sys.argv[1])

mat2, mat_rem = cc1p.make_circC1P(m, sys.argv[2])

f = open(sys.argv[3], 'w')

mat2.write(f.write)

f.close()

f = open(sys.argv[4], 'w')

mat_rem.write(f.write)

f.close()
