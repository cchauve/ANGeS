# ANGES 1.01, reconstructing ANcestral GEnomeS maps
# July 2012.
# Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca

import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))) + '/COMMON')

import bm
import sort
import c1p
import mc1p
import cc1p
import bab
import pqtree
import tree

class BABTester:
	def __init__(self, h):
		self._part_list = []
		self._mem_list = [[[], []] for l in xrange(h)]
	#enddef
	
	def test(self, row, job):
	
		if not job._rem:
			ret, self._mem_list[job._row] = c1p.test_row(row, self._part_list)

			return ret
		#endif
		
		return True
	#enddef
	
	def clean(self, row, job):
		if not job._rem:
			for m in self._mem_list[job._row][1]:
				self._part_list.remove(m)
			#endfor
			
			for m in self._mem_list[job._row][0]:
				sort.insert(self._part_list, m, 0, len(self._part_list), c1p.comp_part);
			#endfor
		#endif
	#enddef
#endclass

class Support:
	def __init__(self, s):
		self._support = s
	#enddef
#endclass

class MC1PTester:
	def __init__(self, m):
		self._tree = None
#		self._part_list = [[] for l in xrange(h + 1)]
		self._part_list = []
		self._mem_list = [[[], []] for l in xrange(m._height)]
		self._r = 0
		self._used = [True for l in xrange(m._height)]
		self._t_copy = [l for l in xrange(m._height)]
		self._support = m.get_support()
		self._telo_rows = []		# minimal telomere rows added so far
		self._telo_rows_all = []		# all telomere rows added so far
		self._telo_contained_in = [[] for l in xrange(m._height)]
		self._matrix = m
		
		for i in xrange(m._height):
			rowi = m.get_row_info(i)
		
			if rowi._isT:
				break
			#endif
		
			for j in xrange(m._height - 1, -1, -1):
				rowj = m.get_row_info(j)
			
				if not rowj._isT:
					break
				#endif
				
				if rowj._set == rowi._set:
					self._t_copy[i] = j
					
					break
				#endif
				
				if rowi._isT and rowi._set < rowj._set:
					self._telo_contained_in[i].append(j)
				#endif
			#endfor
		#endfor
	#enddef
	
	def test(self, row, job):
		if row._isT:				
			if self._tree == None:
#				mat = bm.BinaryMatrix()
				
#				for l in self._part_list[job._row]:
#				for l in self._part_list:
#					for r in l._list:
#						mat.add_row_info(r)
#					#endif
#				#endfor				
				
#				if mat._height > 0:
#					self._tree = c1p.make_PQR_tree(mat)
				# add leaves
#				og = [Support(set([x])) for x in self._support]
#				nodes = [tree.TreeNode(x) for x in self._support]
				nodes = []
				
				# add internal vertices
				for p in self._part_list:
					pcopy = p._part.copy()
					val = mc1p.PQRTreeValue()
					node = tree.TreeNode(val)
				
					val._support = pcopy._support
				
#					og.append(Support(pcopy._support))
										
					if pcopy._part == pcopy._end:
						val._type = 'P'
					else:
						val._type = 'Q'
						
						node._child = pcopy._part
					#endif
					
					nodes.append(node)
				#endfor
				
				nodes = sort.sort(nodes, tree.comp_nodes)

				# delete overlap graphs with the same support
				i = 0		# component iterator
				
				while i < len(nodes) - 1:
					if nodes[i]._value._support == nodes[i + 1]._value._support:
						if nodes[i]._value._type == 'P':
							del nodes[i]
						elif nodes[i + 1]._value._type == 'P':
							del nodes[i + 1]
						else:
							i += 1
						#endif			
					else:
						i += 1
					#endif
				#endwhile
							
				# add root
				val = mc1p.PQRTreeValue()
				
				val._type = 'P'
				val._support = self._support
				
				nodes.append(tree.TreeNode(val))
								
				self._tree = tree.PQRTree()
				
#				self._tree._head = c1p.make_PQR_tree_from_partitions(self._support, nodes, og)
				self._tree._head = mc1p.make_PQR_tree_from_parts(nodes)
				
#				mc1p.fix_tree_node(self._tree._head)
#				else:
#					self._tree = 0
				#endif
				
				self._r = job._row
			#endif
		
			if not job._rem and self._tree != 0:
#				if not any(row._set in l._list for l in self._part_list[self._r]):
#				if not any(row._set == r._set for r in x._list for x in self._part_list):
				if not self._used[job._row]:
					return False
				#endif
				
				minimal = True
				telo_disc = []
		
				# check minimality
				for r in self._telo_rows:
					if r < row._set:
						minimal = False
						
						break
					#endif
				#endfor
				
				if minimal:
					i = 0
				
					# remove non minimal rows
					while i < len(self._telo_rows):
						if row._set < self._telo_rows[i]:
							mc1p.undo_check_LCA_path(self._telo_rows[i], self._tree)
							
							telo_disc.append(self._telo_rows[i])
							
							del self._telo_rows[i]
						else:
							i = i + 1
						#endif
					#endwhile

					if mc1p.check_LCA_path(row._set, self._tree):
						self._telo_rows.append(row._set)
						self._telo_rows_all.append(job._row)
						
						return True
					else:
						mc1p.undo_check_LCA_path(row._set, self._tree)
					
						for r in telo_disc:
							mc1p.check_LCA_path(r, self._tree)
							
							self._telo_rows.append(r)
						#endfor
					#endif
		
					return False
				else:
					self._telo_rows_all.append(job._row)
				
					return True
				#endif
			
#				if mc1p.check_LCA_path(row._set, self._tree):
#					return True
#				#endif
#		
#				mc1p.undo_check_LCA_path(row._set, self._tree)
#		
#				return False
			#endif
#		else:
		elif not job._rem:
#			self._part_list[job._row + 1] = [c1p.PartitionTree(x._part.copy(), copy.copy(x._list)) for x in self._part_list[job._row]]
	
#			if not job._rem:
#				return c1p.test_row(row._set, self._part_list[job._row + 1])
#			#endif
			ret, self._mem_list[job._row] = c1p.test_row(row, self._part_list)
			
			self._used[self._t_copy[job._row]] = ret
			
			return ret
		else:
			self._used[self._t_copy[job._row]] = False
		#endif
		
		return True
	#enddef
	
	def clean(self, row, job):
		if row._isT:
			if job._rem and self._r == job._row:
				self._tree = None
				self._r = 0
			elif not job._rem:
				if row in self._telo_rows:
					mc1p.undo_check_LCA_path(row._set, self._tree)
					
					self._telo_rows.remove(row)
					self._telo_rows_all.remove(job._row)
					
					for i in self._telo_contained_in[job._row]:
						if i in self._telo_rows_all:
							r = self._matrix.get_row_info(i)
							minimal = True
		
							# check minimality
							for r in self,_telo_rows:
								if r < row._set:
									minimal = False
							
									break
								#endif
							#endfor
							
							if minimal:
								mc1p.check_LCA_path(r, tre)
			
								self._telo_rows.append(r)
							#endif
						#endif
					#endfor
				else:
					self._telo_rows_all.remove(job._row)
				#endif
			#endif
#		else:
#			self._part_list[job._row + 1] = []
		elif not job._rem:			
			for m in self._mem_list[job._row][1]:
				self._part_list.remove(m);
			#endfor
			
			for m in self._mem_list[job._row][0]:
				sort.insert(self._part_list, m, 0, len(self._part_list), c1p.comp_part);
			#endfor
		#endif
	#enddef
#endclass

class PQMC1PTester:
	def __init__(self, m, mult = 0):
		# make PQ-tree
		self._tree = pqtree.make_PQR_tree(m)
		mc1p.fix_tree_node(self._tree._head)
		
		self._mult = mult
		
		self._telo_rows = []		# minimal telomere rows added so far
		self._telo_rows_all = []		# all telomere rows added so far
		self._telo_contained_in = [[] for l in xrange(m._height)]
		self._matrix = m
		
		for i in xrange(m._height):
			rowi = m.get_row_info(i)
		
			if rowi._isT:
				break
			#endif
		
			for j in xrange(m._height - 1, -1, -1):
				rowj = m.get_row_info(j)
			
				if not rowj._isT:
					break
				#endif
				
				if rowi._isT and rowi._set < rowj._set:
					self._telo_contained_in[i].append(j)
				#endif
			#endfor
		#endfor
	#enddef
	
	def test(self, row, job):
		if not job._rem:
			minimal = True
			telo_disc = []
		
			# check minimality
			for r in self._telo_rows:
				if r < row._set:
					minimal = False
					
					break
				#endif
			#endfor
			
			if minimal:
				i = 0
			
				# remove non minimal rows
				while i < len(self._telo_rows):
					if row._set < self._telo_rows[i]:
						mc1p.undo_check_LCA_path(self._telo_rows[i], self._tree)
						
						telo_disc.append(self._telo_rows[i])
						
						del self._telo_rows[i]
					else:
						i = i + 1
					#endif
				#endwhile

				if mc1p.check_LCA_path(row._set, self._tree):
					self._telo_rows.append(row._set)
					self._telo_rows_all.append(job._row)
					
					return True
				else:
					mc1p.undo_check_LCA_path(row._set, self._tree)
				
					for r in telo_disc:
						mc1p.check_LCA_path(r, self._tree)
						
						self._telo_rows.append(r)
					#endfor
				#endif
	
				return False
			else:
				self._telo_rows_all.append(job._row)
			
				return True
			#endif
		#endif
		
		return True
	#enddef
	
	def clean(self, row, job):
		if not job._rem:			
			if row in self._telo_rows:
				mc1p.undo_check_LCA_path(row._set, self._tree)
				
				self._telo_rows.remove(row)
				self._telo_rows_all.remove(job._row)
				
				for i in self._telo_contained_in[job._row]:
					if i in self._telo_rows_all:
						r = self._matrix.get_row_info(i)
						minimal = True
									# check minimality
						for r in self,_telo_rows:
							if r < row._set:
								minimal = False
						
								break
							#endif
						#endfor
						
						if minimal:
							mc1p.check_LCA_path(r, tre)
		
							self._telo_rows.append(r)
						#endif
					#endif
				#endfor
			else:
				self._telo_rows_all.remove(job._row)
			#endif
		#endif
	#enddef
#endclass

def C1P_bab(mat, matb = bm.BinaryMatrix(), mat_rem = bm.BinaryMatrix()):
	ms = c1p.split_matrix(mat)		# split matrices
				
	j = 1		# iterator for tracing
		
	for m in ms:
		print str(j) + '/' + str(len(ms))
		
		j += 1
		
		# sort matrix to get heuristic as first answer
		m.sort()
				
		# branch and bound
		# optimal rows to remove to make compnonent C1P
		rows = bab.branch_and_bound(m, False, BABTester(m._height))
			
		rows.sort()
			
		for i in xrange(len(rows) - 1, -1, -1):
			row = m.get_row_info(rows[i])		# row to remove
						
			mat_rem.add_row_info(row)
				
			m.remove_row(rows[i])
		#endfor
		
		# collect usable rows into the C1P matrix
		for r in m._rows:
			matb.add_row_info(r)
		#endfor
		
		print ''
	#endfor
	
	return matb, mat_rem
#enddef

def mC1P_bab(mat, matb = bm.BinaryMatrix(), mat_rem = bm.BinaryMatrix()):
	ms = c1p.make_intersect_components(mat)		# intersect components
	
	j = 1		# iterator for tracing
	
	for m in ms:
		print str(j) + '/' + str(len(ms)) + ' ',
		sys.stdout.flush()
		
		j += 1
		
		# sort matrix to get heuristic as first answer
		m.sort()
		
		# branch and bound
		# optimal rows to remove to make compnonent mC1P
		rows = bab.branch_and_bound(m, False, MC1PTester(m))
		
		rows.sort()
		
		for i in xrange(len(rows) - 1, -1, -1):
			row = m.get_row_info(rows[i])		# row to remove
						
			mat_rem.add_row_info(row)
				
			m.remove_row(rows[i])
		#endfor
		
		# collect usable rows into the C1P matrix
		for r in m._rows:
			matb.add_row_info(r)
		#endfor
		
		print ''
	#endfor
	
	return matb, mat_rem
#enddef

def C1P_and_mC1P_bab(mat, matb = bm.BinaryMatrix(), mat_rem = bm.BinaryMatrix()):
	ms = c1p.make_intersect_components(mat)		# intersect components
	
	j = 1		# iterator for tracing
	
	for m_full in ms:
		print str(j) + '/' + str(len(ms)) + ' ',
		sys.stdout.flush()
		
		j += 1
		
		m = bm.BinaryMatrix()
		m_tel = bm.BinaryMatrix()
		
		# split telomere rows and non telomere rows
		for row in m_full:
			if row._isT:
				m_tel.add_row_info(row)
			else:
				m.add_row_info(row)
			#endif
		#endfor
		
		# C1P bab
		print 'C1P'
		sys.stdout.flush()
		
		st = mat_rem._height
		
		# sort matrix to get heuristic as first answer
		m.sort()
		
		# branch and bound
		# optimal rows to remove to make compnonent mC1P
		rows = bab.branch_and_bound(m, False, MC1PTester(m))
		
		rows.sort()
		
		for i in xrange(len(rows) - 1, -1, -1):
			row = m.get_row_info(rows[i])		# row to remove
						
			mat_rem.add_row_info(row)
				
			m.remove_row(rows[i])
		#endfor
		
		# collect usable rows into the C1P matrix
		for r in m._rows:
			matb.add_row_info(r)
		#endfor

		# remove unused telomere rows
		for i in xrange(st, mat_rem._height):
			row = mat_rem.get_row(i)
		
			for j in xrange(m_tel._height):
				if m_tel.get_row(j) == row:
					m_tel.remove_row(j)
					
					break
				#endif
			#endfor
		#endfor
	
		# mC1P bab
		print 'telomeres'
		sys.stdout.flush()
			
		# sort matrix to get heuristic as first answer
		m_tel.sort()
			
		# branch and bound
		# optimal rows to remove to make compnonent mC1P
		rows = bab.branch_and_bound(m_tel, False, PQMC1PTester(m))
			
		rows.sort()
			
		for i in xrange(len(rows) - 1, -1, -1):
			row = m_tel.get_row_info(rows[i])		# row to remove
							
			mat_rem.add_row_info(row)
					
			m_tel.remove_row(rows[i])
		#endfor
			
		# collect usable rows into the C1P matrix
		for r in m_tel._rows:
			matb.add_row_info(r)
		#endfor
		
		print ''
	#endfor
			
	return matb, mat_rem
#enddef

def circC1P_bab(mat, matb = bm.BinaryMatrix(), mat_rem = bm.BinaryMatrix()):
	circ_mat = cc1p.convert_to_circular(mat)		# circular matrix
	
	# rows to keep and remove
	circ_matb, circ_mat_rem = C1P_bab(circ_mat, bm.BinaryMatrix(), bm.BinaryMatrix())
	
	rem = [r._link for r in circ_mat_rem]		# indices of rows to remove
	
	rem.sort()
	
	for i in xrange(len(rem) - 1, -1, -1):
		row = mat.get_row_info(rem[i])
		
		mat_rem.add_row_info(row)
		
		mat.remove_row(rem[i])
		
		print ''
	#endfor
	
	for r in mat:
		matb.add_row_info(r)
	#endfor
		
	return matb, mat_rem
#enddef
