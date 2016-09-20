# ANGES 1.01, reconstructing ANcestral GEnomeS maps
# July 2012.
# Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca

import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))) + '/COMMON')

import bm
import stak

import math
import copy
from decimal import *

#######################################################
#    bab.py
#
#    Runs branch and bound on a nonC1P matrix to make it
#    C1P.
#
#######################################################

# represent a node in the branch and bound
class Job:
	# initializes the Job object
	# row - int, rem - bool, clean - bool
	def __init__(self, row, rem, clean):
		self._row = row		# row number to process
		self._rem = rem		# True to remove row False to keep
		self._clean = clean		# True marks cleaning job (fix scores)
	#enddef
#endclass

class Logger:
	def __init__(self, height):
		self._height = height
		self._cell = -1
		self._pos = [0 for i in xrange(self._height // 8 + 1)]
		self._changed = False
	#enddef
	
	def add(self, row):
		self._pos[row // 8] += int(math.pow(2, row % 8))
		
		if self._cell >= 0 and row // 8 <= self._cell:
			self._changed = True
			self._cell = row // 8
		#endif
		
		i = row // 8
		
		while i > 0 and self._pos[i] >= 256:
			self._pos[i] -= 256
			self._pos[i - 1] += 1
			
			if self._cell >= 0 and i - 1 <= self._cell:
				self._changed = True
				self._cell = i - 1
			#endif
			
			i -= 1
		#while
	#enddef
		
	def start(self):
		self._cell = len(self._pos) - 1
		self._changed = True
		
		print '&' + str(self._pos),
	#enddef

	def log(self):
		print '&' + str(self._pos),
	#	if self._changed and self._cell >= 0:
	#		print str(self._cell) + ":" + str(self._pos[self._cell]) + ",",
	#		sys.stdout.flush()
	#		
	#		self._changed = False
	#	#endif
	#enddef
#enclass

# applies the branch and bound methodology to find the minumum conflicting set of rows of a matrix, m,
# with the consecutive-ones property using best_score as
# the best score up till now
# m - bm.Matrix, prop - bool, tester - Tester, best_score - decimal.Decimal
# return - list of set; the rows that are optimal to remove
def branch_and_bound(m, prop, tester, best_score = Decimal('0')):
	unused = []		# list of the current rows to remove
	best_unused = []		# list of the optimal rows to remove so far
	score = Decimal('0')		# current score
	stack = stak.Stack()		# job stack
#	min_C1P_sets = []		# list of minimal rows removed that make C1P sets
	st = 0		# start of removable rows
	pc = 0		# per completed
	leaves = 0		# leaves processed / cut
		
	# find start of removable rows
	while st < m._height and m.get_row_info(st)._weight > Decimal('10'):
		tester.test(m.get_row_info(st), Job(st, False, False))
		
		st += 1
	#endwhile
	
#	pos = Logger(m._height - st)
	
	# current score left down the current branch and bound path
#	if prop:
#		score_left = Decimal('0')
#	else:
	score_left = sum([m.get_row_info(i)._weight for i in xrange(st, m._height)]) 
#	#endif
		
	# put first row on stack
	stack.push(Job(st, True, False))

	r = Job(st, False, False)		# current job
		
	while r != None:	
		# done
		if r._row == m._height:
			if best_score < score: #or best_score == Decimal('0'):		# this should always be true but check ayways
				# check that unused is a minimum conflicting set
#				nC1P = True
			
#				for t in unused:
#					row = m.get_row(t)
#				
#					if p[m._height].test(row):
#						p[m._height].clean(row, False)
#						
#						nC1P = False
#						
#						min_C1P_sets.append(set(unused))
#						
#						break
#					#endif				
#				#endfor
#				
#				if nC1P:
				best_score = score
				best_unused = copy.copy(unused)
								
#				min_C1P_sets.append(set(unused))
				
#				if prop:
#					print math.exp(best_score)
#				else:
#				print ''
				print str(best_score)
				sys.stdout.flush()
#				print ''
#				#endif
#				#endif

#				pos.add(m._height - st - 1)
#				pos.start()
			#endif
		else:
			row = m.get_row_info(r._row)		# current Row
					
			if r._clean:		# clean				
#				if not prop:
				score_left = score_left + row._weight
#				#endif
			
				tester.clean(row, r)
			
				if r._rem:
					del unused[-1]
				else:				
#					if prop:
#						score -= math.log(row._weight)
#					else:
					score = score - row._weight
#					#endif
				#endif
			else:
				if (score_left + score > best_score):
#					if not prop:
					score_left = score_left - row._weight
#					#endif
				
					# test row
					ad = tester.test(row, r)
					
					if r._rem:		# remove r
						# true if path is traversible so far
						ad = ((score + score_left) > best_score)# or best_score == Decimal('0')	

						if ad:
							unused.append(r._row)
							
#							# determine if rows removed are minimal
#							for s in min_C1P_sets:
#								if s < set(unused):
#									ad = False
#									
#									del unused[-1]
#					
#									break
#								#endif
#							#endfor
						#endif				
					else:		# add r					
						if ad:
#							if prop:
#								score += math.log(row._weight)
#							else:
							score = score + row._weight
#							#endif
						#endif
					#endif
				else:
					ad = False
				#endif
			
				if ad:
					# add cleaning job
					r._clean = True
					stack.push(r)

					# add children to job list
					if r._row + 1 < m._height:
						stack.push(Job(r._row + 1, True, False))
					#endif
					
					stack.push(Job(r._row + 1, False, False))
				else:					
#					if not prop:
					score_left = score_left + row._weight
#					#endif

					if r._rem:
						tester.clean(row, r)
					#endif

#					pos.add(r._row - st)
#					pos.log()
				#endif
			#endif
		#endif
	
		r = stack.pop()
	#endwhile
		
	return best_unused
#enddef
