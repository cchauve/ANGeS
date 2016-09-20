# ANGES 1.01, reconstructing ANcestral GEnomeS maps
# July 2012.
# Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca

import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))) + '/COMMON')

import bm
import tree
import part
import c1p
import cc1p

import sort

############################################
#     pqtree.py
#
#     Contains code to compute PQR-trees from
#     binary matrices.
#
############################################

# turns 'P'-nodes with 2 children into Q-nodes
# node - tree.TreeNode
def populate_Pnode(node):
	children = 0		# number of children in node so far
	child = node._child		# current child
			
	while child != None and children <= 2:
		children += 1
				
		child = child._brother
	#endwhile
	
	if children == 2:
		node._value = 'Q'
	elif children == 0:
		node._value = node._value.__iter__().next()
	else:	
		node._value = 'P'
	#endif
#enddef

def make_PQR_tree_from_partitions(support, nodes, og):
	ss = len(support)	# the length of the support

	# find the parents of the nodes
	for c in xrange(len(og) - 1):
		S = og[c]._support		# the support of the current overlap component
		s = max(c + 1, ss)		# index to start looking for parents
		
		nodes[c]._support = og[c]._support
						
		for p in xrange(s, len(og)):
			if S <= og[p]._support:
				if nodes[p]._value == 'Q' or nodes[p]._value == 'C':		# Q node
					child = nodes[p]._child		# current child of the current node
				
					while child != None:
						if type(child._value) == set:
							if S < child._value:
								nodes[c]._brother = child._child
								
								if child._child != None:
									child._child._last = nodes[c]
								#endif
				 				
				 				child._child = nodes[c]
				 				
				 				break
							elif S == child._value:
								child._value = nodes[c]._value
								child._child = nodes[c]._child
								nodes[c] = child
								
								break								
							#endif
						#endif
						
						child = child._brother
					#endwhile
				else:	# P or R node
					nodes[c]._brother = nodes[p]._child
					
					if nodes[p]._child != None:
						nodes[p]._child._last = nodes[c]
					#endif
					
					nodes[p]._child = nodes[c]
				#endif
								
				break
			#endif
		#endfor
		
		# change P node to Q if length is 2
		if nodes[c]._value == 'P':
			populate_Pnode(nodes[c])
		elif nodes[c]._value == 'Q' or nodes[c]._value == 'C':		# fix children of Q nodes
			child = nodes[c]._child		# current child of the current node
			
			while child != None:
				if type(child._value) == set:
					populate_Pnode(child)
				#endif
			
				child = child._brother
			#endwhile
		#endif
	#endfor
	
	return nodes[-1]
#enddef

# makes a PQR tree with support, support, from the set of overlap graphs, og
# support - set, og - list of tree.SpanningTree
# return - tree.TreeNode; root of PQR-tree
def make_PQR_tree_from_graph(support, og):
	# add leaves
	for x in support:
		og.insert(0, c1p.OverlapComponent(bm.Row(set([x]), 'leaf', 0, [], False)))
	#endfor

	# create PQR-tree nodes
	nodes = [tree.TreeNode(g._support) for g in og]		# the nodes of the PQR-tree
	
	# for each overlap graph make a node in the PQR-tree
	for g in xrange(len(og)):
		head = og[g]._head		# the head of the current overlap graph
	
		# if there is only one node in the tree the component is a P node otherwise: 
		if len(og[g]._rows) > 0:		# refine to determine order or R'ness
			is_R = False		# True if node is an R-node
			p = part.Partition(head)		# partition of the
														# current overlap graph
			
			for row in og[g]._rows:
				if not p.refine(row):
					is_R = True
					
					break
				#endif
			#endfor
				
			if is_R:
				nodes[g]._value = 'R'
			else:
				nodes[g]._value = 'Q'
				nodes[g]._child = p._part
			#endif
		else:		# P node
			if len(og[g]._support) == 1:
				nodes[g]._value = og[g]._support.__iter__().next()
			else:
				nodes[g]._value = 'P'
			#endif
		#endif
	#endfor
	
	return make_PQR_tree_from_partitions(support, nodes, og)
#enddef

# creates the PQR Tree of the given binary matrix
#
# matrix: the matrix - bm.Binary Matrix
#
# return: the matrix's PQR-tree - tree.PQRTree
def make_PQR_tree(matrix):
	# seperate telomere rows
	telomere_rows = []
	i = 0
	remove_rows = []
	
	while i < matrix._height:
		row = matrix.get_row_info(i)
	
		if row._isT:			
			ad = True		# True if we will add the telomere row
							# to the list of minimal telomere rows
	
			# ensure list of telomere rows is minimal
			for j in xrange(len(telomere_rows) - 1, -1, -1):
				tRow = telomere_rows[j]		# current row in telomere list
		
				if tRow[0] < row._set:
					ad = False
			
					break
				#endif
		
				if row._set < tRow[0]:
					remove_rows.append(tRow[1])
					
					del telomere_rows[j]
				#endif
			#endfor
	
			if ad:
				# add telomere to row
				telomere_rows.append([row._set, i])
				row._set = row._set | set([c1p.Telomere()])
			else:
				remove_rows.append(i)
			#endif
		#endif
		
		i += 1	
	#endfor
	
	remove_rows.sort()
	
	for i in xrange(len(remove_rows) - 1, -1, -1):
		matrix.remove_row(remove_rows[i])
	#endfor

	# create overlap graphs
	og = c1p.create_overlap_graph(matrix)	 # overlap graph (as an array of spanning trees)
	pq = tree.PQRTree()		# pqr tree
	S = matrix.get_support()		# set of all columns (support of the matrix)
	
	# sort components of the overlap graph by size of support
	og.sort()
	
	# add head node
	og.append(c1p.OverlapComponent(bm.Row(S, 'head', 0, [], False)))
	
	# delete empty lines	(will have to reverse)
	while len(og[0]._support) == 0:
		del og[0]
	#endwhile
			
	# delete overlap graphs with the same support
	i = 0		# component iterator
	
	while i < len(og) - 2:
		if og[i]._support == og[i + 1]._support:
			if len(og[i]._rows) == 0:
				del og[i]
			elif len(og[i + 1]._rows) == 0:
				del og[i + 1]
			else:
				i += 1
			#endif			
		else:
			i += 1
		#endif
	#endwhile
	
	# make S sorted
	S = [x for x in S]
	
	S.sort()

	# create PQR tree
	pq._head = make_PQR_tree_from_graph(S, og)		# head of PQ-tree
		
	return pq
#enddef

# makes a PQCR tree with support, support, from the set of overlap graphs, og
# support - set, og - list of tree.SpanningTree
# return - tree.TreeNode; root of PQR-tree
def make_PQCR_tree_from_graph(support, og):
	# add leaves
	for x in support:
		og.insert(0, c1p.OverlapComponent(bm.Row(set([x]), 'leaf', 0, [], False)))
	#endfor

	# create PQR-tree nodes
	nodes = [tree.TreeNode(g._support) for g in og]		# the nodes of the PQR-tree
	
	# for each overlap graph make a node in the PQR-tree
	for g in xrange(len(og)):
		head = og[g]._head		# the head of the current overlap graph
	
		# if there is only one node in the tree the component is a P node otherwise: 
		if len(og[g]._rows) > 0:		# refine to determine order or R'ness
			is_R = False		# True if node is an R-node
			p = part.Partition(head)		# partition of the current overlap graph
			
#			print str(p)

			for row in og[g]._rows:
				if not p.refine(row):
					is_R = True
					
					break
				#endif
			#endfor
				
			if is_R:
				# check if node is circC1P
				for n in xrange(g+1, len(og) - 1):
					if og[g]._support < og[n]._support:
						if n != len(og) - 1:
							nodes[g]._value = 'R'
						#endif
					#endif
				#endfor
				
				if nodes[g]._value != 'R':
					m = bm.BinaryMatrix()
					
					for row in og[g]._all_rows:
						m.add_row_info(row)
					#endfor
					
					if cc1p.check_circC1P(m):
						p_comp = part.Partition(head)
						non_c1p = []
						
						for row in og[g]._rows:
							if p_comp.test_refine(row._set) < 0:
								p_comp.refine(row)
							else:
								non_c1p.append(row)
							#endif
						#endif
						
						for row in non_c1p:
							p_comp.left_refine(row)
							p_comp.right_refine(row)
							
							if len(row._set - p_comp._support) > 0:
								p_comp.insert_after(tree.TreeNode(row._set - p_comp._support), p_comp._end)
								p_comp._support = p_comp._support | row._set
							#endif
						#endfor
						
						nodes[g]._value = 'C'
						nodes[g]._child = p_comp._part
					else:
						nodes[g]._value = 'R'
					#endfor
				else:
					nodes[g]._value = 'R'
				#endif
			else:
				nodes[g]._value = 'Q'
				nodes[g]._child = p._part
			#endif
		else:		# P node
			if len(og[g]._support) == 1:
				nodes[g]._value = og[g]._support.__iter__().next()
			else:
				nodes[g]._value = 'P'
			#endif
		#endif
	#endfor
		
	return make_PQR_tree_from_partitions(support, nodes, og)
#enddef

# creates the PQCR Tree of the given binary matrix
#
# matrix: the matrix - bm.Binary Matrix
#
# return: the matrix's PQR-tree - tree.PQRTree
def make_PQCR_tree(matrix):
#	print "debug 1"
	# create overlap graphs
	og = c1p.create_overlap_graph(matrix)	 # overlap graph (as an array of spanning trees)
#	print "debug 2"
	pq = tree.PQRTree()		# pqr tree
#	print "debug 3"
	S = matrix.get_support()		# set of all columns (support of the matrix)
#	print "debug 4"
		
	# sort components of the overlap graph by size of support
	og.sort()
#	print len(og)
	
	# add head node
	og.append(c1p.OverlapComponent(bm.Row(S, 'head', 0, [], False)))
	
	# delete empty lines	(will have to reverse)
	while len(og[0]._support) == 0:
		del og[0]
	#endwhile
			
	# delete overlap graphs with the same support
	i = 0		# component iterator
	
	while i < len(og) - 2:
		if og[i]._support == og[i + 1]._support:
			if len(og[i]._rows) == 0:
				del og[i]
			elif len(og[i + 1]._rows) == 0:
				del og[i + 1]
			else:
				i += 1
			#endif			
		else:
			i += 1
		#endif
	#endwhile
	
	# make S sorted
	S = [x for x in S]
	
	S.sort()

	
	# create PQR tree
	pq._head = make_PQCR_tree_from_graph(S, og)		# head of PQ-tree
		
	return pq
#enddef
