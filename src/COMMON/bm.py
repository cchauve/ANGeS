# ANGES 1.01, reconstructing ANcestral GEnomeS maps
# July 2012.
# Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca

import sort
from decimal import *

#######################################################
#    bm.py
#
#    Class to encode a binary/ternary matrix. Uses sets to
#    sparsely record the ones/X's of each row.
#
#######################################################

# a row of a binary matrix
class Row:
	# creates a new Row
	#
	# s: set of columns that have a one - set of int
	# ident: the id of the row - str
	# weight: the weight of the row - decimal.Decimal
	# sp: list of species that realizes the row - list of str
	# isT: True if the row is a telomere - bool
	def __init__(self, s, ident, weight, sp, isT):
		self._set = s		# set of columns that have a one at the row
		self._id = ident		# id of the row
		self._weight = weight		# weight of the row
		self._sp = sp		# list of species that realize the row
		self._isT = isT		# True if it is a telomere row
		self._isX = False		# Not an X row
		self._pos = None		# position of consecutive set in species' genome
	#enddef
	
	# returns a string representation of the Row in the format
	# id|weight;species1,species2,...:column1 column2 ... [T]
	#
	# return: string repesentation of the Row - str
	def __str__(self):
		s = self._id + '|' + str(self._weight) + ';'
					
		for p in xrange(len(self._sp)):
			s = s + self._sp[p]
			
			if p < len(self._sp) - 1:
				s = s + ','
			#endif
		#endfor
		
		s = s + ":"
			
		sorted_set = [x for x in self._set]
		
		sorted_set.sort()
					
		for x in sorted_set:
			s = s + str(x) + ' '
		#endfor
		
		if self._isT:
			s = s + 'T '
		#endif
		
		if self._pos != None:
			s = s + 'P ' + str(pos[0]) + ' ' + str(pos[1])
		#endif
		
		return s
	#enddef
#endclass

# a row of a ternary matrix with X's
class XRow:
	# creates a new Row
	#
	# s: set of columns that have a one - set of int
	# Xs: set of columns that have an X - set of int
	# ident: the id of the row - str
	# weight: the weight of the row - decimal.Decimal
	# sp: list of species that realizes the row - list of str
	# isT: True if the row is a telomere - bool
	def __init__(self, s, Xs, ident, weight, sp, isT):
		self._set = s		# set of columns that have a one at the row
		self._id = ident		# id of the row
		self._weight = weight		# weight of the row
		self._sp = sp		# list of species that realize the row
		self._isT = isT		# True if it is a telomere row
		self._isX = True        # X-row
		self._Xs = Xs		# set of all columns that have an X at the row
		self._pos = None		# position of consecutive set in species' genome
	#enddef
	
	# returns a string representation of the Row in the format
	# id|weight;species1,species2,...:column1 column2 ... X xcolumn1 xcolumn2 ... [T]
	#
	# return: string repesentation of the XRow - str
	def __str__(self):
		s = self._id + '|' + str(self._weight) + ';'
					
		for p in xrange(len(self._sp)):
			s = s + self._sp[p]
			
			if p < len(self._sp) - 1:
				s = s + ','
			#endif
		#endfor
		
		s = s + ":"
			
		sorted_set = [x for x in self._set]
		
		sorted_set.sort()
					
		for x in sorted_set:
			s = s + str(x) + ' '
		#endfor
		
		s =  s + 'X '
		
		sorted_set = [x for x in self._Xs]
		
		sorted_set.sort()
		
		for x in sorted_set:
			s = s+ str(x) + ' '
		#endfor
		
		if self._isT:
			s = s + 'T '
		#endif
		
		if self._pos != None:
			s = s + 'P ' + str(pos[0]) + ' ' + str(pos[1])
		#endif
		
		return s
	#enddef
#endclass

# compares Rows by weight
#
# a - Row, b - Row
#
# return - bool
def weight_comp(a, b):
	if a._isT == b._isT:
		return a._weight > b._weight
	elif a._isT:
		return False
	else:
		return True
	#endif
#enddef

def size_comp(a, b):
	return len(a._set) < len(b._set)
#enddef

# stores a binary matrix in a sparse representation
class BinaryMatrix:
	# creates a new empty BinaryMatrix
	def __init__(self):
		self._rows = []		# list of rows in the matrix
		self._height = 0		# height of the matrix
	#enddef
		
	# returns the string representation of the BinaryMatrix where each row is written on its own line in the format:
	# id|weight;species1,species2,species3,...:column1 column2 column3 ... [X xcolumn1 xcolumn2 ...] [T]
	#
	# return: the the string representation of the BinaryMatrix - str
	def __str__(self):
		s = ''
	
		for r in self._rows:		
			s = s + str(r) + '\n'
		#endfor
		
		return s
	#enddef
	
	def write(self, write_fn):
		for r in self._rows:		
			write_fn(str(r) + '\n')
		#endfor
	#enddef
	
	# reads a Binary Matrix in from file with file name, file_name, where each line is a row in the matrix
	# there are 3 different formats accepted for lines
	#
	# weight\tcolumn1 column2 column3 ...
	# id|species1,species2,species3,...:column1 column2 column3 ...
	# id|weight;species1,species2,species3,...:column1 column2 column3 ...
	#
	# file_name: the file name of the binary matrix - str
	def from_file(self, file_name):
		f = file(file_name, 'r')
		line = f.readline()
		
		while line != '':
			if line != '\n':			
				tok = line.split()
				first = tok[0].split('|')
				isT = False
				pos = None
										
				if len(first) > 1:		# standard				
					ident = first[0]
					first = first[1].split(';')
									
					if len(first) > 1:
						weight = Decimal(first[0])
						
						first = first[1].split(':')
					else:		# weightless
						weight = Decimal('0.0')
												
						first = first[0].split(':')
					# endif
				
					sp = first[0].split(',')
					tok[0] = first[1]	
				else:		# old (testing only)
					ident = '0'
					weight = Decimal(tok[0])
					sp = []
					
					del tok[0]
				#endif
				
				if sp == ['']:
					sp = []
				#endif
			
				X = -1
				i = 0
			
				while i < len(tok):
					try:
						tok[i] = int(tok[i])
					except ValueError:
						if (tok[i] == 'T'):
							isT = True
							
							del tok[i]
							
							i -= 1
						elif(tok[i] == 'X'):
							X = i
							
							del tok[i]
							
							i -= 1
						elif(tok[i] == 'P'):
							pos = [tok[i+1], tok[i+2]]
							
							del tok[i]
							del tok[i+1]
							del tok[i+2]
							
							i -= 1
						#endif
					#endtry
					
					i += 1
				#endwhile
			
				if X == -1:
					self.add_row(set(tok), ident, weight, sp, isT)					
				else:
					self.add_row_info(XRow(set(tok[0:X]), set(tok[X:]), ident, weight, sp, isT))
				#endif
				
				self._rows[-1]._pos = pos
			#endif
			
			line = f.readline()
		#endwhile
	#enddef
	
	# gets the jth row (starting at 0)
	#
	# j: the index of the row to get - int
	#
	# return: the set of ones at the row - set of int
	def get_row(self, j):
		return self._rows[j]._set
	#enddef
	
	# gets the info of the jth row (starting at 0)
	#
	# j: the index of the row to get - int
	#
	# return: the row - Row/XRow
	def get_row_info(self, j):
		return self._rows[j]
	#enddef
		
	# adds row onto the bottom of the matrix
	#
	# row: set of columns with a 1 - set of int
	# ident: id of the row - str
	# weight: weight of the row - decimal.Decimal
	# sp: species that realize the row - list of string
	# isT: True if the row is a telomere - bool
	def add_row(self, row, ident = '-1', weight = Decimal('0'), sp = [], isT = False):
		if ident == '-1':
			ident = str(self._height)
		#endif
	
		self._rows.append(Row(row, ident, weight, sp, isT))
			
		self._height += 1
	#enddef
	
	# adds a row onto the bottom of the matrix 
	#
	# row: the row to add - Row/XRow
	def add_row_info(self, row):
		self._rows.append(row)
		
		self._height += 1
	#enddef
	
	# removes the jth row (starting at 0)
	#
	# j: index of the row to remove - int
	def remove_row(self, j):
		del self._rows[j]
		
		self._height -= 1
	#enddef
	
	# gets the set of columns used by the matrix
	#
	# return: the set of columns of the matrix - set
	def get_support(self):
		s = set([])
		
		for r in self._rows:
			s = s | r._set
		#endfor
		
		return s
	#enddef
	
	# sorts the rows by weight
	def sort(self):
		sort.sort(self._rows, weight_comp)
	#endef
	
	def sort_size(self):
		sort.sort(self._rows, size_comp)
	#enddef
	
	# returns the sum of the weights of the rows of the matrix
	#
	# return - int
	def get_weight(self):
		w = 0
		
		for r in self._rows:
			w += r._weight
		#endfor
		
		return w
	#endef
	
	# iterates over the rows of the matrix
	#
	# return: iterator over the rows - listiterator
	def __iter__(self):
		return iter(self._rows)
	#enddef
#endclass

# set decimal precision
getcontext().prec = 20
