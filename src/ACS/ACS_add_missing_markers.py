# ANGES 1.01, reconstructing ANcestral GEnomeS maps
# July 2012.
# Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca

import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))) + '/COMMON')

from acs import *

#######################################################
#    ACS_add_missing_markers.py
#
#    Adds missing markers to ancestral consecutive syntenies.
#    if markers are not universal
#######################################################

if len(sys.argv) != 7:
	print 'Adds missing markers to an ACS file by either adding 1\'s or X\'s.'
	print ''
	print 'Usage:'
	print 'ACS_add_missing_markers.py MarkersFile Species1 Species2 InputACSFile [1/2/3] OutputACSFile'
	print '   Input:'
	print '    MarkersFile - the markers file'
	print '    Species1, Species2 - the species pair considered'
	print '    InputACSFile - the ACS file to add missing markers to'
	print '   Output:'
	print '    OutputACSFile - the ACS file with missing markers'
	print '   Options:'
	print '    1 - Missing markers spanned in the interval of an ACS are added as 1s to this ACS'
	print '    2 - Missing markers present in either genome are added as Xs'
	print '    3 - Mix option 1 and option 2: markers spanned by intervals are added as 1s, markers missing in both species are added as Xs'
	sys.exit()
#endif

markers, blocks = extract_blocks(sys.argv[1])

spe1 = sys.argv[2]
spe2 = sys.argv[3]

ACS = bm.BinaryMatrix()

ACS.from_file(sys.argv[4])

if sys.argv[5].upper() == '1':
	ACS_missing = add_missing_markers_1(blocks, spe1, spe2, ACS)
elif sys.argv[5].upper() == '2':
	ACS_missing = add_missing_markers_X(blocks, markers, spe1, spe2, ACS)
elif sys.argv[5].upper() == '3':
	ACS_missing = add_missing_markers_1(blocks, spe1, spe2, ACS)
	ACS_missing = add_missing_markers_X(blocks, markers, spe1, spe2, ACS_missing, False)
#endif

f = open(sys.argv[6], 'w')

f.write(str(ACS_missing))

f.close()
