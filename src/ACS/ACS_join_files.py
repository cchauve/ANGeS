# ANGES 1.01, reconstructing ANcestral GEnomeS maps
# July 2012.
# Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca

import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))) + '/COMMON')

import bm
from acs import *

import gc

#######################################################
#    ACS_join_files.py
#
#    Combines ACS files into a single file
#
#######################################################

if len(sys.argv) < 3:
	print 'Combines the rows of inputs ACS files into one output ACS file'
	print ''
	print 'Usage:'
	print 'ACS_join_files.py [0/1] InputACSFile1 InputACSFile2 ... OutputACSFile'
	print '   Input:'
	print '    InputACSFile1, InputACSFile1 ... - the ACS files to join together'
	print '   Output:'
	print '    OutputACSFile -  the joined ACS files'
	print '   Options:'
	print '    0 - do not merge ACS with the same content and species support'
	print '    1 - merge ACS with the same content and species support'
	print ''
	print 'The default is option 1'
	
	sys.exit()
#enddef

if sys.argv[1] == '0':
	out = file(sys.argv[-1], 'w')		# joined file

	for i in xrange(2, len(sys.argv) - 1):
		in_matrix = open(sys.argv[i])		# matrix file
		
		out.write(in_matrix.read())
		
		in_matrix.close()
	#endfor
	
	out.close()
else:
	if sys.argv[1] == '1':
		st = 2
	else:
		st = 1
	#endif
	
	in_matrices = [bm.BinaryMatrix() for i in xrange(st, len(sys.argv) - 1)]    

	for m in xrange(st, len(sys.argv) - 1):
		in_matrices[m-st].from_file(sys.argv[m])
	#endfor

	out_matrix = join_ACS(in_matrices)		# joined syntenies
	
	out = file(sys.argv[-1], 'w')		# joined file
	
	out_matrix.write(out.write)
	
	out.close()
#endif

