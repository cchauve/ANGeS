# ANGES 1.01, reconstructing ANcestral GEnomeS maps
# July 2012.
# Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca

import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))) + '/COMMON')

import bm
import tree
import c1p

import copy
import sort

###########################################################
#    mc1p.py
#
#    C1P with telomeres
#
###########################################################

# info stored in a PQR-tree node
class PQRTreeValue:
	# creates a new PQRTreeValue
	def __init__(self):
		self._type = None		# type of node (P, Q or R)
		self._support = set([])		# the leaves of the node
		self._allowed = set([])		# allowed telomere adjacencies
		self._children_used = []		# list of the children that have 
										# been used in previous paths
		self._used = 0		# the number of times this node has been used
	#enddef
	
	def __copy__(self):
		cpy = PQRTreeValue()
		
		cpy._type = self._type
		cpy._support = copy.copy(self._support)
		cpy._children_used = copy.copy(self._children_used)
		cpy._used = self._used
		cpy._first = copy.copy(self._first)
		cpy._last = copy.copy(self._last)
		
		return cpy
	#enddef
	
	def __str__(self):
		return "(" + str(self._type) + "," + str(self._allowed) + ", " + str(self._children_used) + ", " + str(self._used) + ")"
	#enddef
#endclass

def populate_node(node):
	if node._value._type == 'Q':
		first = node._child
		
		while node._child != None and type(node._child._value) == set:
			last = node._child
			
			if len(node._child._value) > 1:
				val = PQRTreeValue()
				val._type = 'P'
				val._support = node._child._value
				
				s =  node._child._value
				node._value._allowed = set([])
				chil = node._child._child
		
				while chil != None:
					val._allowed |= chil._value._allowed
					
					s = s - chil._value._support
					chil = chil._brother
				#endwhile
				
				val._allowed |= s
				
				node._child._value = val
			else:
				node._child = node._child._brother
				
				if node._child != None:
					node._child._last = None
				#endif
			#endif
		#endwhile
		
		child = node._child
		
		while child != None:
			last = child
			
			if type(child._value) == set:
				if len(child._value) > 1:
					val = PQRTreeValue()
					val._type = 'P'
					val._support = child._value
				
					s =  child._value
					node._value._allowed = set([])
					chil = child._child
		
					while chil != None:
						val._allowed |= chil._value._allowed
					
						s = s - chil._value._support
						chil = chil._brother
					#endwhile
				
					val._allowed |= s
				
					child._value = val
				else:
					if child._brother != None:
						child._brother._last = child._last
					#endif
					
					if child._last != None:
						child._last._brother = child._brother
					#endif						
				#endif
			#endif
			
			child = child._brother
		#endwhile
		
		if not type(first._value) == set:
			node._value._allowed = first._value._allowed
		else:
			node._value._allowed = first._value
		#endif
		
		if not type(last._value) == set:
			node._value._allowed = node._value._allowed | last._value._allowed
		else:
			node._value._allowed = node._value._allowed | last._value
		#endif
	else:
		s = node._value._support
		node._value._allowed = set([])
		child = node._child
		
		while child != None:
			node._value._allowed |= child._value._allowed
					
			s = s - child._value._support
			child = child._brother
		#endwhile
		
		node._value._allowed |= s
	#endif
#enddef

def make_PQR_tree_from_parts(nodes):
	# find the parents of the nodes
	for c in xrange(len(nodes)):
		S = nodes[c]._value._support		# the support of the current overlap component
						
		for p in xrange(c + 1, len(nodes)):
			if S <= nodes[p]._value._support:
				if nodes[p]._value._type == 'Q':		# Q node
					child = nodes[p]._child		# current child of the current node
				
					while child != None:
#						if type(child._value) == list:
						if type(child._value) == set:
#							if S < child._value[0]:
							if S < child._value:
								nodes[c]._brother = child._child
								
								if child._child != None:
									child._child._last = nodes[c]
								#endif
				 				
				 				child._child = nodes[c]
				 				
				 				break
#							elif S == child._value[0]:
							elif S == child._value:
								child._value = nodes[c]._value
								child._child = nodes[c]._child
								nodes[c] = child
								
								break								
							#endif
						#endif
						
						child = child._brother
					#endwhile
				else:		# P or R node
					nodes[c]._brother = nodes[p]._child
					
					if nodes[p]._child != None:
						nodes[p]._child._last = nodes[c]
					#endif
					
					nodes[p]._child = nodes[c]
				#endif
								
				break
			#endif
		#endfor
		
		populate_node(nodes[c])
	#endfor
	
	return nodes[-1]
#enddef

def fix_tree_node(node):
	# set node value
	val = node._value
	first = None
	last = None
	node._value = PQRTreeValue()
	node._value._type = val
	
	child = node._child
	
	# fix all children
	while child != None:
		# fix child
		fix_tree_node(child)
		
		# if child is leaf
		if len(child._value._support) == 0:
			# add leaf value to support
			node._value._support |= set([child._value._type])
			
			if node._value._type != 'Q':
				node._value._allowed |= set([child._value._type])
			#endif
			
			last = set([child._value._type])
			
			if first == None:
				first = set([child._value._type])
			#endif
			
			# delete leaf
			if child == node._child:
				node._child = child._brother
			else:
				child._last._brother = child._brother
				
				if child._brother != None:
					child._brother._last = child._last
				#endif
			#endif
		else:
			# add child's support
			node._value._support |= child._value._support
			
			if node._value._type != 'Q':
				node._value._allowed |= child._value._allowed
			#endif
			
			last = child._value._allowed
			
			if first == None:
				first = child._value._allowed
			#endif
		#endif
		
		child = child._brother
	#endwhile
	
	if node._value._type == 'Q':
		node._value._allowed = first | last
	#endif
#enddef

# finds the next node in a PQR-tree
# tre - str, i - int
# return - int
def next_node(tre, i):
	try:
		index_p = tre.index('_P', i)		# location of the next P node
	except:
		index_p = -1
	#endtry
	
	try:
		index_q = tre.index('_Q', i)		# location of the next Q node
	except:
		index_q = -1
	#endtry
	
	try:
		index_r = tre.index('_R', i)		# location of the next R node
	except:
		index_r = -1
	#except
		
	if index_p < 0:
		index_p = max(index_q, index_r)
	#endif
		
	if index_q < 0:
		index_q = max(index_p, index_r)
	#endif
		
	if index_r < 0:
		index_r = max(index_p, index_q)
	#endif
		
	return min(index_p, index_q, index_r)
#enddef

# recursive part of parse_tree
# tre - list of str, i - int, current - tree.TreeNode
# return - set of int, int
def parse_tree_rec(tre, i, current):
	current._value._type = tre[i][1]
	child = None		# current child
	first = None
	last = None
	
	i += 1
	
	# find end of node and next node
	end = tre.index(current._value._type + '_', i)		# end of node
	n = next_node(tre, i)		# index of next node
	
	while end > n and n > 0:
		markers = tre[i:n]		# children between next node and current
		
		# add markers to support
		if len(markers) > 0:			
			if first == None:
				first = set([int(markers[0])])
			#endif
			
			current._value._support |= set([int(m) for m in markers])
		#endif
		
		# make child and recurse it
		if child == None:
			current._child = tree.TreeNode(PQRTreeValue())
			child = current._child
		else:
			child._brother = tree.TreeNode(PQRTreeValue())
			child = child._brother
		#endif
		
		a, s, i = parse_tree_rec(tre, n, child)		# child's support, current index
		
		if current._value._first == None:
			first = s
		#endif
		
		last = s
		
		current._value._support |= s
		
		# refind end and next node
		end = tre.index(current._value._type + '_', i)
		n = next_node(tre, i)
	#endwhile
	
	markers = tre[i:end]	# children at end of node
	
	# add markers to support
	if len(markers) > 0:		
		if first == None:
			first = set([int(markers[0])])
		#endif
		
		last = set([int(markers[-1])])
		
		if current._value._type != 'Q':	
			current._value._support |= set([int(m) for m in markers])
		#endif
	#endif
	
	if current._value._type == 'Q':
		current._value._allowed = first | last
	#endif 
	
	return current._value._allowed, current._value._support, end + 1
#enddef

# parses a tree file
# NOT DONE
# tree_file - str
# returns tree.PQRTree
def parse_tree(tree_file):
	f = file(tree_file, 'r')		# tree file

	tre = tree.PQRTree()		# PQR-tree
	
	tre._head._value = PQRTreeValue()
	tre._head._value._type = 'P'	

	first = True		# True if curr is the first node

	for buff in f:
		if buff[0] != '#' and buff[0] != '>':
			if first:
				tre._head._child = tree.TreeNode(PQRTreeValue())
				curr = tre._head._child		# current tree
				
				first = False
			else:
				curr._brother = tree.TreeNode(PQRTreeValue())
				curr = curr._brother		# current tree
			#endif
		
			a, s, i = parse_tree_rec(buff.split(), 0, curr)		# aux, support, aux var
			
			tre._head._value._support |= s
		#endif
	#endfor
	
	f.close()
	
	return tre
#enddef

# finds and check path from head to LCA of row, r
# if fails, run exit(False)
# r - set of int, tre - tree.PQRTree
def check_LCA_path(r, tre):
	current = tre._head		# current node
	
	# Step 3ci* (for head)
#	if tre._head._value._type == 'Q' and not tre._head._value._first <= r and not tre._head._value._first > r and not tre._head._value._last <= r and not tre._head._value._last > r:
	if len(r & tre._head._value._allowed) == 0:
		return False
	#endif

	while True:
		child = current._child		# current child
		current._value._used += 1
		
		while child != None:
			if r <= child._value._support:				
				max_paths = 0		# maximum paths allowed according to steps 3cii-vi
				
				# Step 3cii,iii
				if tre._head._value._type == 'Q' or tre._head != current:
					max_paths = 1
				# Step 3civ
				elif tre._head == current:
					max_paths = 2
				#endif
				
				found = False		# True if child has been used in a path before
				
				for e in current._value._children_used:
					if e[0] == child:
						e[1] += 1
						found = True
						
						# Step 3cii,iii,iv
						if max_paths > 0 and e[1] > max_paths:
							return False
						#endif
					#endif
				#endfor
					
				if not found:
					current._value._children_used.append([child, 1])
				#endif
				
				# fix error in the Algorithm when the LCA was a Q node and the row was not at the end or beginning of the LCA's children
				# Step 3ci*
#				if child._value._type == 'Q' and not child._value._first <= r and not child._value._first > r and not child._value._last <= r and not child._value._last > r:
				if len(r & child._value._allowed) == 0:
					return False
				#endif
				
				current = child
				
				break
			#endif
			
			child = child._brother
		#endwhile
		
		if child == None:
			break
		#endif
	#endwhile
		
	return True
#enddef

# undoes check_LCA_path
# r - set of int, tre - tree.PQRTree
def undo_check_LCA_path(r, tre):
	current = tre._head		# current node
	
	# Step 3ci* (for head)
#	if tre._head._value._type == 'Q' and not tre._head._value._first <= r and not tre._head._value._first > r and not tre._head._value._last <= r and not tre._head._value._last > r:
	if len(r & tre._head._value._allowed) == 0:
		return
	#endif
	
	while True:
		child = current._child		# current child
		current._value._used -= 1
	
		while child != None:
			if r <= child._value._support:				
				max_paths = 0		# maximum paths allowed according to steps 3cii-vi
				
				# Step 3cii,iii
				if tre._head._value._type == 'Q' or tre._head != current:
					max_paths = 1
				# Step 3civ
				elif tre._head == current:
					max_paths = 2
				#endif
				
				found = False		# True if child has been used in a path before
				
				for e in current._value._children_used:
					if e[0] == child:
						e[1] -= 1

						if e[1] == 0:						
							current._value._children_used.remove(e)
						# Step 3cii,iii,iv
						elif max_paths > 0 and e[1] + 1 > max_paths:
							return
						#endif
					#endif
				#endfor
				
				# fix error in the Algorithm when the LCA was a Q node and the row was not at the end or beginning of the LCA's children
				# Step 3ci*
#				if child._value._type == 'Q' and not child._value._first <= r and not child._value._first > r and not child._value._last <= r and not child._value._last > r:
				if len(r & child._value._allowed) == 0:
					return
				#endif
				
				current = child
				
				break
			#endif
			
			child = child._brother
		#endwhile
		
		if child == None:
			break
		#endif
	#endwhile
#enddef

# calculates K1 and K2 according to Step 5a
# tre - tree.PQRTree
# return - int, int
def count_used(tre):
	child = tre._head._child		# current child
	K1 = 0		# K1
	K2 = 0		# K2
	
	while child != None:
		if child._value._used == 1:
			K1 += 1
		elif child._value._used == 2:
			K2 += 1
		#endif
		
		child = child._brother
	#endwhile
	
	return K1, K2
#enddef

def test_telomere_row(row, data):
	minimal = True
	telo_rows = data[0]
	tre = data[1]
	telo_disc = []

	# check minimality
	for r in telo_rows:
		if r < row._set:
			minimal = False
			
			break
		#endif
	#endfor
	
	if minimal:	
		i = 0
	
		# remove non minimal rows
		while i < len(telo_rows):
			if row._set < telo_rows[i]:
				undo_check_LCA_path(telo_rows[i], tre)
				
				telo_disc.append(telo_rows[i])
								
				del telo_rows[i]
			else:
				i = i + 1
			#endif
		#endwhile

		if check_LCA_path(row._set, tre):
			telo_rows.append(row._set)
			
			return True, 0
		else:
			undo_check_LCA_path(row._set, tre)
		
			for r in telo_disc:
				check_LCA_path(r, tre)

				telo_rows.append(r)
			#endfor
		#endif
		
		return False, 0
	else:
		return True, 0
	#endif
#endif

#def test_telomere_row(row, tre):
#	if check_LCA_path(row._set, tre):
#		return True, 0
#	#endif
#		
#	undo_check_LCA_path(row._set, tre)
#	
#	return False, 0
##endif

def remove_row(row, data):
	return False, 0
#enddef

# removes rows until the matrix is C1P
# matrix - bm.BinaryMatrix
# return - list of set
def make_mC1P(matrix):
	p = []		# list of partitions and spanning trees used
	tre = None
	i = 0		# iterator
	rows = []		# rows to remove to make C1P
	telomere = [j for j in xrange(matrix._height)]
	use = [True for j in xrange(matrix._height)]
	support = matrix.get_support()
	telo_rows = []

	# sort the matrix by weight
	matrix.sort()
	
	# process each row
	while i < matrix._height:
		row = matrix.get_row_info(i)
	
		if row._isT:				
			if tre == None:
				nodes = []
				
				# add internal vertices
				for part in p:
					val = PQRTreeValue()
					node = tree.TreeNode(val)
				
					val._support = part._part._support
				
					if part._part._part == part._part._end:
						val._type = 'P'
					else:
						val._type = 'Q'
						
						node._child = part._part._part
					#endif
					
					nodes.append(node)
				#endfor
				
				nodes = sort.sort(nodes, tree.comp_nodes)

				# delete overlap graphs with the same support
				j = 0		# component iterator
				
				while j < len(nodes) - 1:
					if nodes[j]._value._support == nodes[j + 1]._value._support:
						if nodes[j]._value._type == 'P':
							del nodes[j]
						elif nodes[j + 1]._value._type == 'P':
							del nodes[j + 1]
						else:
							j += 1
						#endif			
					else:
						j += 1
					#endif
				#endwhile
							
				# add root
				val = PQRTreeValue()
				
				val._type = 'P'
				val._support = support
				
				nodes.append(tree.TreeNode(val))
								
				tre = tree.PQRTree()
				
				tre._head = make_PQR_tree_from_parts(nodes)
			#endif
			
			if use[i]:
				test_row = test_telomere_row
#				data = tre
				data = [telo_rows, tre]
			else:
				test_row = remove_row
				data = None
			#endif
		else:
			test_row = c1p.test_row
			data = p
			
			for j in xrange(i+1, matrix._height):
				if matrix.get_row_info(j)._set == row._set:
					telomere[i] = j
				#endif
			#endfor
		#endif
		
		yes, temp = test_row(row, data)		# true if test_row succeeds, aux
		
		use[telomere[i]] = yes
		
		if not yes:
			# remove row from matrix
			rows.append(i)			
		#endif
				
		i += 1
	#endwhile
		
	return rows
#enddef
