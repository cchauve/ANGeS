# ANGES 1.01, reconstruction ANcestral GEnomeS maps
# June 2012.
# Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca

import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))) + '/COMMON')

import bm
import c1p

##############################################################
#    cc1p
#
#    Circular C1P.
#
##############################################################


# converts a binary matrix to a circular binary matrix using j as reference. 
# if j is None then j is arbitrarily chosen from mat.
#
# mat: the matrix to convert - bm.BinaryMatrix
# j: the reference column(s) - int, set
#
# return: the circular matrix - bm.BinaryMatrix
def convert_to_circular(mat, j = None):
	if j == None:
		j = iter(mat.get_row(0)).next()
	#end

	cmat = bm.BinaryMatrix()		# circular matrix
	s = mat.get_support()		# set of columns of mat
	
	for i in xrange(mat._height):
		row = mat.get_row_info(i)
		
		if type(j) == set:
			if j < row._set:
				cmat.add_row(s - row._set, row._id + "C", row._weight, row._sp, row._isT)
			else:
				cmat.add_row(row._set, row._id, row._weight, row._sp, row._isT)
			#endif
		else:
			if j in row._set:
				cmat.add_row(s - row._set, row._id + "C", row._weight, row._sp, row._isT)
			else:
				cmat.add_row(row._set, row._id, row._weight, row._sp, row._isT)
			#endif
		#endif
				
		cmat._rows[-1]._link = i
	#endfor
	
	return cmat
#endef

# determines whether a matrix's maximal overlap components are circ-C1P or not
#
# mat: the matrix to check - bm.BinaryMatrix
#
# return: True if the matrix's maximal overlap components are circ-C1P - bool
def check_circC1P(mat):
	maxs = c1p.make_intersect_components(mat)		# maximal connected components
	
	for m in maxs:
		cm = convert_to_circular(m) # circular overlap component
		
		if not c1p.check_C1P(cm):
			return False
		#endif
	#endfor
		
	return True
#enddef

# makes a matrix circC1P using a heuristic which considers rows of higher 
# weight first
#
# mat: the matrix to make circC1P - bm.BinaryMatrix
# style: 'whole' to make the matrix as a whole circC1P, 'max' to make the
#    maximal overlap components circC1P - str
#
# return: the circC1P matrix - bm.BinaryMatrix 
def make_circC1P(mat, style):
	mat_rem = bm.BinaryMatrix()
	c1pmat = bm.BinaryMatrix()

	if style == 'whole':
		cmat = convert_to_circular(mat)
		
		mat.sort()
				
		remove = c1p.make_C1P(cmat)
				
		for i in xrange(len(remove) - 1, -1, -1):
			mat_rem.add_row_info(mat.get_row_info(remove[i]))
		
			mat.remove_row(remove[i])
		#endfor
		
		c1pmat = mat
	elif style == 'max':
		maxs = c1p.make_intersect_components(mat)
		
		for m in maxs:
			cm = convert_to_circular(m) # circular overlap component
		
			m.sort()
		
			remove = c1p.make_C1P(cm)
				
			for i in xrange(len(remove) - 1, -1, -1):
				mat_rem.add_row_info(m.get_row_info(remove[i]))
			
				m.remove_row(remove[i])
			#endfor
			
			for row in m:
				c1pmat.add_row_info(row)
			#endfor
		#endfor
	#endif
	
	return c1pmat, mat_rem
#enddef
