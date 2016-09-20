# ANGES 1.01, reconstructing ANcestral GEnomeS maps
# July 2012.
# Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca

import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))) + '/COMMON')

import bm
import tree
import part
import sort

import copy

#######################################################
#    c1p.py: linear C1P
#
#    Testing the C1P (linear mode), heuristic for the
#    linear C1P
#
#######################################################

# class for sorting
class PartitionTree:
	# initializes a PartitionTree
	# p - part.Partition, l - list of set
	def __init__(self, p, l):
		self._part = p		# partition
		self._list = l		# list of sets that make up the partition via refinement
	#enddef
#enddef

class OverlapComponent:
	def __init__(self, param):
		if type(param) == list:
			self._head = param[0]
			self._rows = param[1:]
			self._all_rows = param
			self._support = set([i for r in param for i in r._set])
		else:
			self._head = param
			self._rows = []
			self._all_rows = [param]
			self._support = param._set
		#endif
		
		self._nsupport = len(self._support)
	#enddef
	
	def add_row(self, row):
		self._rows.append(row)
		self._all_rows.append(row)
		self._support = self._support | row._set
		
		self._nsupport = len(self._support)
	#enddef
	
	def join_components(self, comp, pos):
		self._rows.extend(comp._all_rows)
		self._all_rows.extend(comp._all_rows)
		self._support = self._support | comp._support
		
		self._nsupport = len(self._support)
	#enddef
	
	def __lt__(self, a):
		if self._nsupport == a._nsupport:
			return min(self._support) < min(a._support)
		#endif	
		
		return self._nsupport < a._nsupport
	#enddef
	
	def __le__(self, a):
		if self._nsupport == a._nsupport:			
			return min(self._support) <= min(a._support)
		#endif
	
		return self._nsupport < a._nsupport
	#enddef
#endclass

# represents a Telomere
class Telomere:
	def __init__(self):
		pass
	#enddif

	def __str__(self):
		return 'T'
	#enddef
#endclass

# compares partition trees
def comp_part(a, b):
	return len(a._part._support) > len(b._part._support)
#enddef

# tests whether two sets overlap
# a - set, b - set
# return - bool
def overlap(a, b):
	if len(a & b) == 0:
		return False
	#endif
	
	if a < b or b < a:
		return False
	#endif
	
	return True
#enddef

# auxiliary 'f' method for traverse
# returns true if v._value and params['v']._value do not overlap
# v - vertex, params = {'og': list of tree.SpanningTree, 'j': int,
# 'v': tree.Vertex, 'found': int}
# return - bool
def test_overlap(v, params):
	if (overlap(v._value._set, params['v']._value._set)):
		if params['found'] >= 0:		# join trees
			params['og'][params['found']].join(params['og'][params['j']], v, params['v'])
			
			del params['og'][params['j']]
		else:		# add to tree
			params['og'][params['j']].add_vertex(params['v'], v)
		#endif
		
		return False
	#endif
		
	return True
#enddef

# auxiliary 'g' method
# flips a tree around so e is the head with child v
# v - vertex, e - vertex, b - bool, params = {'found': int}
def flip_tree_node(v, e, b, params):
	if not b and params['found'] >= 0:
		v.remove_edge(e)
		e.add_edge(v)
	#endif
	
	return b
#enddef

# auxilary 'f' method for traverse
# refines part with v._value
# v - vertex, part - part.Partition
# return - bool
def refine(v, part):
	return part.refine(v._value)
#enddef

# auxiliary 'g' method for traverse
# does nothing but passthrough b
# v - unused, e - unused, b - A, params - unused
# return - A
def do_nothing(v, e, b, params):
	return b
#endef

# auxiliary 'f' method for traverse
# adds the row, v._value, to matrix
# v - vertex, matrix - bm.BinaryMatrix
# return - bool; True
def make_matrix(v, matrix):
	matrix.add_row_info(v._value)

	return True
#endef

# auxiliary 'f' method for traverse
# adds the row, v._value, to the list params
# v - vertex
# params - list 
#
# return - bool; True
def add_to_list(v, params):
	params.append(v._value)
	
	return True
#enddef

# traverses a tree with head v
# calling f(v, params) before traversing v's children
# and g(v, child, child_return_code, params) after traversing each child
# if f and g return False the traversal exits
# params is used as parameters for both f and g
# v - tree.Vertex
# f - boolean(tree.Vertex, params)
# g - boolean(tree.Vertex, tree.Vertex, boolean, params)
# params - anything
#
# return - bool
def traverse(v, f, g, params):
	stack = [[f, None, v]]
	go = True

	while len(stack) > 0:
		job = stack.pop()
	
		if job[0] == f:		# prefix
			if go:
				go = f(job[2], params)
						
				if job[1] != None:
					stack.append([g, job[1], job[2]])
				#endif
				
				if go: 
					for e in job[2]._edges:
						stack.append([f, job[2], e])
					#endfor
				#endif
			#endif
		else:		# postfix
			g(job[1], job[2], go, params)
		#endif
	#endwhile
	
	return go
#enddef

# mem, sort
# tests if the row, row, can be added to the sets of partitions p and
# still be C1P
# row - set, p - list of PartitionTree, mem - list of rows
# return - bool
def test_row(row, p):
	found = False		# True if first intersecting part is found
	failed = False		# True if not C1P
	j = 0		# current index
	mem = [[], []]		# altered rows
		
	# find overlaps
	while j < len(p):
		if len(p[j]._part._support & row._set) > 0:
			overlaps = False		# True if current component overlaps with row
		
			for s in p[j]._list:
				if overlap(row._set, s._set):
					overlaps = True
					
					break
				#endif
			#endfor
			
			if overlaps:
						
				# refine the new found partition
				if not found:
					pNew = PartitionTree(p[j]._part.copy(), copy.copy(p[j]._list))
					found = True
					failed = not pNew._part.refine(row)
						
					if not failed:
						pNew._list.append(row)
						
						mem[0].append(p[j])
						mem[1].append(pNew)
					#endif
				else:		# if partition already found join the two partitions
					# error (this should not happen)
					if len(pNew._part._part._value) == 0:
						print 'error1'
						print row
					
						return
					#endif
					
					pCopy = PartitionTree(p[j]._part.copy(), copy.copy(p[j]._list))
					failed = not pNew._part.join(pCopy._part)
						
					if not failed:
						pNew._list.extend(p[j]._list)
						mem[0].append(p[j])
					#endif
				#endif
				
				# error (this should not happen)
				if len(pNew._part._part._value) == 0:
					print 'error2'
					print row
					
					return
				#endif
				
				if failed:
					break
				#endif
			#endif
		j = j + 1
	#endwhile
		
	# create new partition (new overlap component)
	if not found:
		par =  PartitionTree(part.Partition(row), [row])
	
		sort.insert(p, par, 0, len(p), comp_part)
		mem = [[], [par]]
	elif not failed:					
		# remove refined rows
		for m in mem[0]:
			p.remove(m)
		#endfor

		# add refined row
		sort.insert(p, pNew, 0, found, comp_part)
	else:
		mem = [[], []]
	#endif
		
	return not failed, mem
#endef

# refines an overlapping matrix in the current order of rows
# m - bm.BinaryMatrix
# return - bool
def refine_matrix(m):
	p = part.Partition(m.get_row_info(0))
	C1P = True
		
	for i in xrange(1, m._height):		
		if not p.refine(m.get_row_info(i)):
			return False					
		#endif
	#endfor
	
	return True
#enddef

# creates the overlap graph (as spanning trees) for a given matrix, matrix
#
# ALGORITHM PREMISE:
# OG = set of trees with vertex set a subset of the rows of matrix.
# At any step in the alogithm each tree in OG will have the property that each vertex overlaps with each of its children.
# At the end of the algorithm OG will be a set of spanning trees of the connect components of the overlap graph of matrix
# and each tree's infix order will give a partitio refinement order.
# 
# matrix - bm.BinaryMatrix
#
# return - list of SpanningTree
def create_overlap_graph(matrix):
	og = []		# overlap graph as a list of tree.SpanningTree

	for i in xrange(matrix._height):
		row = matrix.get_row_info(i)		# current row
		f = -1		# index of overlap graph found for row (-1 = undefined)
		v = tree.Vertex(row)		# tree.Vertex for row to
									# be put in spanning tree
		j = 0			# index of overlap graph
			
		# try to add to existing graph
		while j < len(og):		
			params = {'og': og, 'j': j, 'v': v, 'found': f}
			
			if len(og[j]._support & row._set) > 0 and not traverse(og[j]._head, test_overlap, flip_tree_node, params):
				if f < 0:
					f = j
					j += 1		
				#endif
			else:
				j += 1
			#endif
		#endfor
		
		# make new spanning tree
		if f < 0:
			og.append(tree.SpanningTree(v))
		#endif
	#endfor
	
	rows = []
	
	for g in og:
		rows.append([])
	
		traverse(g._head, add_to_list, do_nothing, rows[-1])
	#endfor
	
	return [OverlapComponent(r) for r in rows]
#endfor

# checks whether the binary matrix, matrix, has the consectutive-ones property 
# by using an algorithm akin to the algorithm of Section 5 of (McConell 2004)
# matrix - bm.BinaryMatrix
# return - bool
def check_C1P(matrix):
	# create overlap graphs
	og = create_overlap_graph(matrix)	# overlap graph (as a spanning tree)
				
	# traverse overlap graph and refine partitions
	for j in xrange(len(og)):
		p = part.Partition(og[j]._head)		# current partition to refine
				
		for v in og[j]._rows:
			if not p.refine(v):
				return False
			#endif
		#endfor
	#endfor
	
	return True
#enddef

def make_intersect_components(matrix):
	mats = []
	sups = []
	
	for r in matrix._rows:
		found = False
	
		for i in xrange(len(mats)):
			if len(sups[i] & r._set) > 0:
				mats[i].add_row_info(r)
				
				sups[i] = sups[i] | r._set
				
				found = True
				
				for k in xrange(len(mats) - 1, i, -1):
					if len(sups[k] & r._set) > 0:
						for row in mats[k]._rows:
							mats[i].add_row_info(row)
						#endfor
						
						sups[i] = sups[i] | sups[k]
						
						del sups[k]
						del mats[k]
					#endif
				#endfor
				
				break
			#endif
		#endfor
		
		if not found:
			mats.append(bm.BinaryMatrix())

			mats[-1].add_row_info(r)
			
			sups.append(r._set)
		#endif
	#endfor
	
	return mats
#enddef

# takes a matrix, m, and splits it into matrices that each only have one 
# connected component
# m - BinaryMatrix
# return list of bm.BinaryMatrix
def split_matrix(m):
	# create overlap graph
	g = create_overlap_graph(m)		# overlap graph
	ms = []		# split matrices
	
	for c in g:
		# make a new matrix for each overlap graph
		n = bm.BinaryMatrix()		# current overlap component
		
		for r in c._all_rows:
			n.add_row_info(r)
		#endfor
		
		ms.append(n)
	#endfor
	
	return ms
#enddef

# removes rows until the matrix is C1P
# matrix - bm.BinaryMatrix
# return - list of set
def make_C1P(matrix):
	p = []		# list of partitions and spanning trees used
	i = 0		# iterator
	rows = []		# rows to remove to make C1P

	# sort the matrix by weight
	matrix.sort()
		
	# process each row
	while i < matrix._height:	
		row = matrix.get_row_info(i)		# current row
		yes, temp =  test_row(row, p)		# true if test_row succeeds, aux
		
		if not yes:
			# remove row from matrix
#			print row._set
						
			rows.append(i)
		#endif
		
		i += 1
	#endwhile
		
	return rows
#enddef
