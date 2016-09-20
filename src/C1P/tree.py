# ANGES 1.01, reconstructing ANcestral GEnomeS maps
# July 2012.
# Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca

import sys
import copy

#######################################################
#    tree.py
#
#    Classes to encode trees.
#
#######################################################

# vertex of a directed graph / tree / list
class Vertex:
	# creates a new Vertex
	# val - anything
	def __init__(self, val):
		self._value = val		# value of the node
		self._edges = set([])		# the set of neighbours/children of the node
		
	#enddef
	
	# adds an edge to the vertex
	# vert - Vertex
	def add_edge(self, vert):
		self._edges.add(vert)
	#enddef
	
	# adds removes an edge from the vertex
	# vert - Vertex
	def remove_edge(self, vert):
		self._edges.remove(vert)
	#enddef
	
	def __hash__(self):
		return id(self)
	#enddef
#endclass

# represents a graph
class Graph:
	# creates a new Graph with head, vert
	# vert - Vertex
	def __init__(self):
		self._head = vert		# the head of the graph (main entry point)
	#enddef
#endclass

# tree with children stored in sets
class SpanningTree:
	# creates a new SpanningTree with head, vert
	# vert - Vertex
	def __init__(self, vert):
		self._head = vert		# the head of the tree
		self._support = vert._value._set		# the support of the tree
	#enddef
	
	# adds a vertex, vert, to the spanning tree at the vertex, parent
	# parent - Vertex, vert - Vertex
	def add_vertex(self, vert, parent):
		parent.add_edge(vert)
		self._support = self._support | vert._value._set
	#enddef
	
	# joins the spanning tree, tree, at the vertex, v, to the spanning tree
	# at the vertex, parent (does not flip the joining spanning tree)
	# tree - SpanningTree, vertex - Vertex
	def join(self, tree, vertex, parent):
		parent.add_edge(vertex)
		self._support = self._support | tree._support
	#enddef
#endclass

# node in a child-brother tree
class TreeNode:
	# creates a new TreeNode with value, val
	# val - anything
	def __init__(self, val):
		self._value = val		# the node's value
		self._child = None		# the node's child
		self._brother = None		# the node's brother
		self._last = None		# the previous brother (only needed
								# for joining partitions)
	#enddef
	
	# returns a string of the node and all its brothers and children
	# return - str
	def __str__(self):	
		if self._child != None:
			# the string representation of the TreeNode
			s = '_' + str(self._value) + ' ' + self._child.__str__() + str(self._value) + '_ '
		else:
			s = str(self._value) + ' '
		#endif
		
		if self._brother != None:
			s += self._brother.__str__()
		#endif
		
		return s
	#enddef
	
	# returns a string of the node and all its children
	# return - str
	def str_self(self):	
		if self._child != None:
			# the string representation of the TreeNode
			s = '_' + str(self._value) + ' ' + self._child.__str__() + str(self._value) + '_ '
		else:
			s = str(self._value)
		#endif
				
		return s
	#enddef
	
	def write(self, write_fn):	
		if self._child != None:
			write_fn('_' + str(self._value) + ' ')
			self._child.write(write_fn)
			write_fn(str(self._value) + '_ ')
		else:
			write_fn(str(self._value) + ' ')
		#endif
		
		brother = self._brother
		
		while brother != None:
			brother.write_self(write_fn)
			
			brother = brother._brother
		#endif
	#enddef
	
	def write_self(self, write_fn):	
		if self._child != None:
			write_fn('_' + str(self._value) + ' ')
			self._child.write(write_fn)
			write_fn(str(self._value) + '_ ')
		else:
			write_fn(str(self._value)+ ' ')
		#endif
	#enddef
	
	# returns the old string of the node and all its brothers and children
	# return - str
	def oldstr(self):	
		if self._child != None:
			# the string representation of the TreeNode
			s = '_' + str(self._value) + self._child.oldstr() + str(self._value) + '_'
		else:
			s = str(self._value) + ','
		#endif
		
		if self._brother != None:
			s += self._brother.oldstr()
		#endif
		
		return s
	#enddef
	
	def __copy__(self):
		node = TreeNode(copy.copy(self._value))
		node._child = copy.copy(self._child)
		node._brother = copy.copy(self._brother)
		
		if node._brother != None:
			node._brother._last = node
		#endif
		
		return node
	#enddef
#endclass

# a child-brother tree
class PQRTree:
	# creates a new PQRTree
	# support - set
	def __init__(self):
		self._head = TreeNode(None)		# the head of the PQRTree
	#enddef
	
	# returns a string of the PQR Tree in its linear representation
	# return - str
	def __str__(self):
		if self._head._value != 'P':
			s = '#CAR0\n'
			s += self._head.str_self() + '\n'
		else:		
			s = ''		
		
			car = self._head._child
			i = 1
			
			while car != None:
				s += '#CAR' + str(i) + '\n'
				s += car.str_self() + '\n'
				
				car = car._brother
				i += 1		
			#endfor
		#endif
				
		return s
	#enddef
	
	def write(self, write_fn):
		if self._head._value != 'P':
			write_fn('#CAR0\n')
			self._head.write_self(write_fn)
			write_fn('\n')
		else:
			car = self._head._child
			i = 1
			
			while car != None:
				write_fn('#CAR' + str(i) + '\n')
				car.write_self(write_fn)
				write_fn('\n')
				
				car = car._brother
				i += 1		
			#endfor
		#endif
	#enddef
		
	def __copy__(self):
		tre = PQRTree()
		
		tre._head = copy.copy(self._head)
		
		return tre
	#enddef
#endclass

# compares two tree nodes' supports
# a - TreeNode, b - TreeNode
# return - bool
def comp_nodes(a, b):
	if len(a._value._support) == len(b._value._support):
		return min(a._value._support) < min(b._value._support)
	#endif
	
	return len(a._value._support) < len(b._value._support)
#endef

# compares two spanning trees' supports
# a - SpanningTree, b - SpanningTree
# return - bool
def comp_trees(a, b):
	if len(a._support) == len(b._support):
		return min(a._support) < min(b._support)
	#endif
	
	return len(a._support) < len(b._support)
#endef
