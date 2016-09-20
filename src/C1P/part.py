# ANGES 1.01, reconstructing ANcestral GEnomeS maps
# July 2012.
# Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca

import tree

import copy

#######################################################
#    part.py
#
#    Partition refinement
#
#######################################################

# tried recording rows used in each class so rows could be removed from the
# partition, but concluded removing rows would be inefficient

class Partition:
	# creates a new Partition with one class, r
	#
	# r: class to add - bm.BinaryMatrix
	def __init__(self, r):
		if r != None:
#			self._part = tree.TreeNode([r._set, set([r])])		# head of partition
			self._part = tree.TreeNode(r._set)		# head of partition
			self._end = self._part		# tail of partition
			self._support = r._set		# support of the partition
		#endif
	#enddef
	
	# string representation of the Partition where each class is written 
	# on a seperate line
	# return - str
	def __str__(self):
		cur = self._part		# current class
		s = ''		# string representation of partition
	
		while cur != None:
#			for x in cur._value[0]:
			for x in cur._value:
				s = s + str(x) + ' '
			#endfor
			
			s = s + '\n'
			cur = cur._brother
		#endwhile
		
		return s
	#enddef
	
	# swap classes a and b
	# a - tree.TreeNode, b -tree.TreeNode
	def swap(self, a, b):
		last = a._last		# class before swapping classes
		brother = b._brother		# class after swapping classes
		
		if b == self._end:
			self._end = a
		else:
			brother._last = a			
		#endif
		
		if a == self._part:
			self._part = b
		else:
			last._brother = b
		#endif
		
		a._last = b
		a._brother = brother
		b._brother = a
		b._last = last
	#enddef
	
	# insert class a before class b
	# a - tree.TreeNode, b -tree.TreeNode
	def insert_before(self, a, b):
		last = b._last		# class before one to be added
		
		if b == self._part:
			self._part = a
		else:
			last._brother = a
		#endif
		
		a._brother = b
		a._last = last
		b._last = a
	#endef
	
	# insert class a after class b
	# a - tree.TreeNode, b -tree.TreeNode
	def insert_after(self, a, b):
		brother = b._brother		# class after one to be added
		
		if b == self._end:
			self._end = a
		else:
			brother._last = a
		#endif
		
		a._last = b
		a._brother = brother
		b._brother = a	
	#endef
	
	# tests if a partition can be refined by a set and if it fails gives the 
	# position of the class that it failed to refine
	# assumes there is an intersection of s and a class of self
	# s - set
	#return - int
	def test_refine(self, s):
		start_full = False		# True if the first interval element is full
		start = -1		# first intersecting item (-1) for unset
		end = -1		# last intersecting item (-1) for unset
		node = self._part	# iterated node
		prevNode = None		# node before iterated node (for insertion before node)
		end_n = 0		# last class in refining set
		
		# find first intersection
		while node != None and len(s) > 0:
#			x = node._value[0]		# class at current node
			x = node._value		# class at current node
			
			y = s & x		# new refined class
			
			if len(y) > 0:		# s intersects x
				start = node
				s = s - y
			
				if len(y) == len(x):		# x = y
					start_full = True
				#endif
				
				end_n += 1
								
				prevNode = node
				node = node._brother
				
				break
			#endif
			
			end_n += 1
			prevNode = node
			node = node._brother
		#endwhile
		
		# look for non full y's
		while node != None and len(s) > 0:
#			x = node._value[0]		# class at current node
			x = node._value		# class at current node
			y = s & x		# new refined class
						
			if len(y) > 0:
				s = s - y
			#endif
			
			if len(y) < len(x):		# y is strictly contained in x	
				s = s - y
						
				end = node
				
				if len(y) > 0:
					end_n += 1
					
					prevNode = node
					node = node._brother
				#endif
				
				break
			#endif
			
			end_n += 1		
			prevNode = node
			node = node._brother
		#endwhile
		
		if end == -1:
			end = node
		#endif
		
		gap = end_n		# position of the end of the gap (if there is one)
		
		# check that part of s doesn't appear later in the partition
		while node != None:
#			if  len(s & node._value[0]) > 0:		# gap found
			if  len(s & node._value) > 0:		# gap found
				return gap
			#endif
						
			node = node._brother
			gap += 1
		#endwhile
		
		# add remainder of s onto the partition
		if len(s) > 0:
			# check if able to add to end
			if end != None:
				# check if able to add to beginning
				if start == self._part and (start_full or prevNode == self._part._brother):
					pass
				else:		# can't add to the begining				
					return end_n		# gap found
				#endif
			#endif
		#endif
		
		return -1
	#enddef
	
	# refines a partition with a set
	# this code is based on Algorithm 1 from (Habib et al. 2000)
	# assumes there is an intersection of s and a class of self
	# r - bm.Row
	# return - bool
	def refine(self, r):
		start_full = False		# True if the first interval element is full
		start = -1		# first intersecting item (-1) for unset
		end = -1		# last intersecting item (-1) for unset
		node = self._part	# iterated node
		prevNode = None		# node before iterated node (for insertion before node)
		s = r._set		# set or ones in the row

		# heuristic to refine the end of the partition
#		if len(self._end._value[0] & s) > 0 and len(self._part._value[0] & s) == 0:
		if len(self._end._value & s) > 0 and len(self._part._value & s) == 0:
			node = self._end
			
#			x = node._value[0]		# class at current node
			x = node._value
			y = s & x		# new refined class
			
			# find intersection of first node
			if x < s:
#				node._value[1] = node._value[1] | set([r])
				start_full = True
			elif len(self._support & (s - x)) == 0:
				if len(s - x) > 0:
#					self.insert_before(tree.TreeNode([x - y, node._value[1]]), self._end)
					self.insert_before(tree.TreeNode(x - y), self._end)
					
#					self._end._value = [y, node._value[1] | set([r])]
					self._end._value = y
				#endif
				
				# add													
				if len(s - y) > 0:
#					self.insert_after(tree.TreeNode([s - y, set([r])]), self._end)
					self.insert_after(tree.TreeNode(s - y), self._end)
				#endif
				
				# update support
				self._support = self._support | s
				
				return True
			else:
#				self._end._value = [x - y, node._value[1]]
				self._end._value = x - y
#				node = tree.TreeNode([y, node._value[1] | set([r])])
				node = tree.TreeNode(y)
				
				self.insert_before(node, self._end)
			#endif
			
#			s = s - node._value[0]
			s = s - node._value
			node = node._last
			
			# look for end of refinement
			while node != self._part and len(s) > 0:
#				x = node._value[0]		# class at current node
				x = node._value		# class at current node
				y = s & x		# new refined class
				
				s = s - y
				
				if len(y) < len(x):
					if len(y) > 0:
#						node._value = [y, node._value[1] | set([r])]
						node._value = y
						
#						self.insert_before(tree.TreeNode([x - y, node._value[1]]), node)
						self.insert_before(tree.TreeNode(x - y), node)
					#endif
					
					break
				#endif
				
				node = node._last
			#endwhile
			
			# add partition to end if possible
			if (len(self._support & s) == 0 and start_full) or len(s) == 0:
				if len(s) > 0:
#					self.insert_after(tree.TreeNode([s, set([r])]), self._end)
					self.insert_after(tree.TreeNode(s), self._end)
				#endif
				
				# update support
				self._support = self._support | s
				
				return True
			else:			
				# gap found
				return False
			#endif
		#endif

		# find first intersection
		while node != None and len(s) > 0:
#			x = node._value[0]		# class at current node
			x = node._value		# class at current node
			
			y = s & x		# new refined class
			
			if len(y) > 0:		# s intersects x
				start = node
				s = s - y
			
				if len(y) == len(x):		# x = y
					start_full = True
#					node._value[1] = node._value[1] | set([r])
				else:
#					node._value[0] = x - y
					node._value = x - y
					
					# breakup partition
#					t = tree.TreeNode([y, node._value[1] | set([r])])
					t = tree.TreeNode(y)
					
					self.insert_after(t, node)

					prevNode = t
					node = t._brother
					
					break
				#endif
				
				prevNode = node
				node = node._brother
				
				break
			#endif
			
			prevNode = node
			node = node._brother
		#endwhile
		
		# look for non full y's
		while node != None and len(s) > 0:
#			x = node._value[0]		# class at current node
			x = node._value		# class at current node
			y = s & x		# new refined class
						
			if len(y) > 0:
				s = s - y
			#endif
			
			if len(y) < len(x):		# y is strictly contained in x
				if len(y) > 0:
#					node._value[0] = x - y
					node._value = x - y
					
					# breakup partition 
#					p = tree.TreeNode([y, node._value[1] | set([r])])
					p = tree.TreeNode(y)
					
					self.insert_before(p, node)
									
					prevNode = p
				#endif
				
				end = node
				
				break
#			else:
#				node._value[1] = node._value[1] | set([r])
			#endif
		
			prevNode = node
			node = node._brother
		#endwhile
		
		if end == -1:
			end = node
		#endif
		
		# add remainder of s onto the partition
		if len(s) > 0:
			# check if able to add to end
			if end != None:
				# check if able to add to beginning
				if start == self._part and (start_full or prevNode == self._part._brother):
					# check that part of s doesn't appear later in the partition
					if len(s & self._support) > 0:		# gap found					
						return False
					#endif
					
					if start_full:
						# add to the begining						
#						self.insert_before(tree.TreeNode([s, set([r])]), self._part)
						self.insert_before(tree.TreeNode(s), self._part)
					elif prevNode == self._part._brother:
						# switch 0th and 1st entries
						self.swap(self._part, self._part._brother)
						
						# add to the begining						
#						self.insert_before(tree.TreeNode([s, set([r])]), self._part)
						self.insert_before(tree.TreeNode(s), self._part)
					#endif
				else:		# can't add to the begining				
					return False		# gap found
				#endif
			else:		# add s to the end
#				self.insert_after(tree.TreeNode([s, set([r])]), self._end)
				self.insert_after(tree.TreeNode(s), self._end)
			#endif
		else:
			# switch the partition to the begining in this case to allow 
			# overlapping partitions to add onto the begining
			if start == self._part and not start_full and prevNode == self._part._brother and self._part._brother._brother != None:
				# switch 0th and 1st entries
				self.swap(self._part, self._part._brother)
			#endif
		#endif
		
		# update support
		self._support = self._support | s
		
		return True
	#enddef
	
	def right_refine(self, row):
		node = self._end
		
		while node != self._part and node._value <= row._set:
			node = node._last
		#endwhile
		
		x = node._value & row._set
		
		if len(x) > 0 and x != node._value:
			y = node._value - x
			node._value = y
			
			self.insert_after(tree.TreeNode(x), node)
		#endif
	#enddef
	
	def left_refine(self, row):
		node = self._part
		
		while node != self._end and node._value <= row._set:
			node = node._brother
		#endwhile
		
		x = node._value & row._set
		
		if len(x) > 0 and x != node._value:
			y = node._value - x
			node._value = y
			
			self.insert_before(tree.TreeNode(x), node)
		#endif
	#enddef
	
	# joins another partition to this one via partition refinement with partitions
	# p - Partition
	# return - bool
	def join(self, p):
		if p._support <= self._support:		# refine inside
			c = self._part		# current class
			
			# find first intersection and direction
			while c != None:			
#				if len(c._value[0] & p._part._value[0]) > 0:
				if len(c._value & p._part._value) > 0:
					# backwards
#					if len(c._value[0] & p._end._value[0]) > 0 and len(p._end._value[0]) == len(c._value[0] & p._end._value[0]):
					if len(c._value & p._end._value) > 0 and len(p._end._value) == len(c._value & p._end._value):
					
						direction = False		# direction to iterate in
						d = p._end		# 'start' of the refining partition
#						s = d._value[0]		# class at the 'start' of the
#											# refining partition
						s = d._value		# class at the 'start' of the
											# refining partition
						end = c		# marks the class that will be the end 
									# of this class once refined
					# forwards
					else: 	
						direction = True
						d = p._part
#						s = d._value[0]
						s = d._value
						end = c
					#endif
					
					break
				#backwards
#				elif len(c._value[0] & p._end._value[0]) > 0:
				elif len(c._value & p._end._value) > 0:
					direction = False
					d = p._end
#					s = d._value[0]
					s = d._value
					end = c
					
					break
				#endif
				
				c = c._brother
			#endwhile
			
			# refine the first intersection
#			while len(s & c._value[0]) > 0:
			while len(s & c._value) > 0:
#				x = tree.TreeNode([s & c._value[0], d._value[1] | c._value[1]])		# class to add
				x = tree.TreeNode(s & c._value)		# class to add
				
				self.insert_after(x, end)
								
				end = x
				
#				s = s - x._value[0]
				s = s - x._value
#				c._value[0] = c._value[0] - x._value[0]
				c._value = c._value - x._value
						
				if len(s) > 0:
					break
				else:
					if direction:
						d = d._brother
					else: 
						d = d._last
					#end
					
					if d != None:
#						s = d._value[0]
						s = d._value
					else:
						break
					#endif
				#endif
			#endwhile
					
			# remove empty partition
#			if len(c._value[0]) == 0:
			if len(c._value) == 0:
				if c._brother != None:
					c._brother._last = c._last
				#endif
				
				if c._last != None:
					c._last._brother = c._brother
				else:
					self._part = c._brother
				#endif
			#endif
					
			c = end._brother
					
			# refine the rest of partition
			while c != None and d != None:		
#				while len(s & c._value[0]) > 0:
				while len(s & c._value) > 0:	
#					x = tree.TreeNode([s & c._value[0], c._value[1] | d._value[1]])		# class to add
					x = tree.TreeNode(s & c._value)		# class to add
					
					self.insert_before(x, c)
											
#					s = s - x._value[0]
					s = s - x._value
#					c._value[0] = c._value[0] - x._value[0]
					c._value = c._value - x._value
											
					if len(s) > 0:
						break
					else:
						if direction:
							d = d._brother
						else: 
							d = d._last
						#end
						
						if d != None:
#							s = d._value[0]
							s = d._value
						else:
							break
						#endif
					#endif
				#endwhile
				
#				if d != None and len(c._value[0]) > 0:		# gap found
				if d != None and len(c._value) > 0:		# gap found
					return False
				else:				
					# remove empty partition
#					if len(c._value[0]) == 0:
					if len(c._value) == 0:
						if c._last != None:
							c._last._brother = c._brother
						#endif
						
						if c._brother != None:
							c._brother._last = c._last
						else:
							self._end = c._last
						#endif
					#endif
				#endwhile
						
				c = c._brother
			#endwhile
			
			if d != None:
				return False
			#endif
#		elif self._support <= p._support:		# refine outside
#			if not p.join(self):
#				return False
#			#endif
#			
#			self._part = p._part
#			self._end = p._end
#		elif len(self._part._value[0] & p._support):		# add to start
		elif len(self._part._value & p._support):		# add to start
			# flip if needed
#			if p._part._value[0] <= self._part._value[0]:
			if p._part._value <= self._part._value:
				p.flip()
#			elif len(self._part._value[0] & p._end._value[0]) == 0:
			elif len(self._part._value & p._end._value) == 0:			
				return False		# gap found
			#endif
			
			# refine the end of the joining partition
			c = p._end		# current node
			
			# find full classes
#			while c != None and c._value[0] <= self._part._value[0]:
			while c != None and c._value <= self._part._value:
#				self._part._value[0] = self._part._value[0] - c._value[0]
				self._part._value = self._part._value - c._value
#				c._value[1] = c._value[1] | self._part._value[1]
								
				c = c._last
			#endif
			
			# split last class
#			if len(c._value[0] & self._part._value[0]) > 0:
			if len(c._value & self._part._value) > 0:
#				x = c._value[0] & self._part._value[0]	# class to add
				x = c._value & self._part._value	# class to add
#				self._part._value[0] = self._part._value[0] - x
				self._part._value = self._part._value - x
#				c._value[0] = c._value[0] - x
				c._value = c._value - x
				
#				p.insert_after(tree.TreeNode([x, c._value[1] | self._part._value[1]]), c)
				p.insert_after(tree.TreeNode(x), c)
			#endif
			
#			if len(p._support & self._part._value[0]) > 0:
			if len(p._support & self._part._value) > 0:
				return False		# gap found
			#endif
			
			# attach partition to start
#			if len(self._part._value[0]) == 0:
			if len(self._part._value) == 0:
				p._end._brother = self._part._brother
				self._part._brother._last = p._end
			else:
				p._end._brother = self._part
				self._part._last = p._end
			#endif
			
			self._part = p._part 			
#		elif len(self._end._value[0] & p._support):		# add to end
		elif len(self._end._value & p._support):		# add to end
			# flip if needed
#			if p._end._value[0] <= self._end._value[0]:
			if p._end._value <= self._end._value:
				p.flip()
#			elif len(self._end._value[0] & p._part._value[0]) == 0:
			elif len(self._end._value & p._part._value) == 0:
				return False		# gap found
			#endif
			
			# refine the end of the joining partition
			c = p._part		# current node
			
			# find full classes			
#			while c != None and c._value[0] <= self._end._value[0]:
			while c != None and c._value <= self._end._value:
#				self._end._value[0] = self._end._value[0] - c._value[0]
				self._end._value = self._end._value - c._value
#				c._value[1] = c._value[1] | self._end._value[1]
				
				c = c._brother
			#endif
			
			# split last class
#			if len(c._value[0] & self._end._value[0]) > 0:
			if len(c._value & self._end._value) > 0:
#				x = c._value[0] & self._end._value[0]		# class to add
				x = c._value & self._end._value		# class to add
#				self._end._value[0] = self._end._value[0] - x
				self._end._value = self._end._value - x
#				c._value[0] = c._value[0] - x
				c._value = c._value - x
				
#				p.insert_before(tree.TreeNode([x, c._value[1] | self._end._value[1]]), c)
				p.insert_before(tree.TreeNode(x), c)
			#endif
						
#			if len(p._support & self._end._value[0]) > 0:
			if len(p._support & self._end._value) > 0:
				return False		# gap found
			#endif
			
			# attach partition to end
#			if len(self._end._value[0]) == 0:
			if len(self._end._value) == 0:
				p._part._last = self._end._last
				self._end._last._brother = p._part
			else:
				p._part._last = self._end
				self._end._brother = p._part
			#endif
			
			self._end = p._end
		else:		# gap found			
			return False
		#endif
		
		self._support = self._support | p._support
		
		return True
	#enddef
	
	# unused/incomplete
	# removes row from the partition
	# row - bm.Row
	# return - list of Partition
	def remove(row):
		st = self._part		# start of row
		r = set([row])		# set contating row
		
		if self._part == self._end and st._val[1] == r:
			return None
		#endif
				
		# find start
		while st != None and not row in st._val[1]:
			st = st._brother
		#endwhile
		
		# row not found
		if st == None:
			return self
		#endif
		
		st._val[1] = st._val[1] - r
		
		curr = st._brother		# current class
		end = st		# end of row
		
		# find end
		while curr != None and row in curr._val[1]:
			end = curr
			curr._val[1] = curr._val[1] - r
		
			curr = curr._brother
		#endwhile
		
		# fix start
		if st != self._part and st._val[1] == st._last._val[1]:
			st._last._val[0] = st._val[0] | st._last._val[0]
			
			self.insert_before(st._last, st._brother)
		#endif
		
		# fix end
		if end != self._end and end._val[1] == end._brother._val[1]:
			end._brother._val[0] = end._val[0] | end._brother._val[0]
			
			self.insert_after(end._brother, end._last)
		#endif
				
		# seperate partitions
	#enddef
	
	# makes a copy of the Partition
	def copy(self):
		p = Partition(None)		# the copy
						
		p._part = tree.TreeNode(copy.copy(self._part._value))
		p._support = self._support
		
		c = self._part._brother		# current node to copy
		d = p._part		# node to copy to
		
		while c != None:
			d._brother = tree.TreeNode(copy.copy(c._value))
			d._brother._last = d
			
			c = c._brother
			d = d._brother
		#endwhile
		
		p._end = d
		
		return p
	#enddef
		
	# flips the partition so _part and _end switch
	def flip(self):
		p = self._part		# current node
		
		self._part = self._end
		self._end = p
		
		while p != None:
			x = p._brother
			p._brother = p._last
			p._last = x
			p = x
		#endwhile
	#enddef
#endclass
