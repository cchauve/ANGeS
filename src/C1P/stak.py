# ANGES 1.01, reconstructing ANcestral GEnomeS maps
# July 2012.
# Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca

##################################################
#     stak.py
#
#     Class to encode a linked-list Stack.
#
##################################################

# a linked-list node
class Node:
	# initializes the Node
	def __init__(self, data, next):
		self._value = data		# data
		self._next = next		# next node in list
	#enddef
#endclass

# a stack data structure implemented as a linked-list
class Stack:
	# initailizes an empty Stack object
	def __init__(self):
		self._head = None		# head of stack (next to be popped)
	#enddef
	
	# peeks the value of the next item
	# return - object; None if the stack is empty	
	def peek(self):
		if self._head == None:
			return None
		else:			
			return self._head._value
		#endif
	#enddef
			
	# pops a item off the stack
	# return - object; None if the stack is empty
	def pop(self):
		if self._head == None:
			return None
		else:
			n = self._head
			self._head = self._head._next
			
			return n._value
		#endif
	#enddef

	# pushes data onto the stack
	# data - object
	def push(self, data):
		if data != None:
			self._head = Node(data, self._head)
		#enddef
	#enddef
#endclass
