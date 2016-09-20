# ANGES 1.01, reconstructing ANcestral GEnomeS maps
# July 2012.
# Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca

#######################################################
#    sort.py
#
#    Abstract sorting functions to sort lists and do 
#    insertion sorting.
#
#######################################################


# returns true if a < b
# sorts ascending
# a - anything, b - anything
# return - bool
def less(a, b):
	return a < b
#enddef

# returns true if a > b
# sort descending
# a - anything, b - anything
# return - bool
def more(a, b):
	return a > b
#enddef

# unused
# auxialiary method for merge sort, merges to sorted lists into one sorted list
# a - list of A, b - list of A, comp - boolean(A, A)
# return - list of int
def merge(a, b, comp):
	order = range(len(a) + len(b))
	j = 0
			
	for i in xrange(len(a)):
		while j < len(b) and not comp(a[i],b[j]):
			order[i + j] = b[j]
			
			j += 1
		#endwhile
		
		order[i + j] = a[i]
	#endfor
	
	while j < len(b):	
		order[len(a) + j] = b[j]
		
		j += 1
	#endwhile
	
	return order
#endef

#unused
# recursive function that sorts a list using merge sort from start to end
# l - list of A, start - int, end - int, comp - boolean(A, A)
# return list of A
def merge_sort(l, start, end, comp):
	# base case, length of interval <= 2
	if end - start == 1:
		return [l[start]]
	elif end - start == 2:
		if comp(l[start], l[start + 1]):
			return [l[start], l[start + 1]]
		else:
			return [l[start + 1], l[start]]
		#endif
	elif end - start == 0:
		return []
	#endif
	
	# split and merge
	return merge(merge_sort(l, start, (start + end) / 2, comp), merge_sort(l, (start + end) / 2, end, comp), comp)
#enddef

# iterative function that sorts a list using merge sort
# l - list of A, comp - boolean(A, A)
# return - list of A
def merge_sort_iter(l, comp):
	n = len(l)		# length of l
	step = 2		# step number (power of 2)
	
	while step <= n:
		# split array into sections and merge them
		for bound in xrange(n // step):
			a = l[step * bound:step * (2 * bound + 1) // 2]		# left block
			b = l[step * (2 * bound + 1) // 2:step * (bound + 1)]	# right block
			j = 0		# iterator
			
			for i in xrange(len(a)):
				while j < len(b) and not comp(a[i], b[j]):
					l[step * bound + i + j] = b[j]
			
					j += 1
				#endwhile
		
				l[step * bound + i + j] = a[i]
			#endfor
	
			while j < len(b):	
				l[step * bound + len(a) + j] = b[j]
		
				j += 1
			#endwhile
		#endfor
		
		# merge end into the last merged section (if l not divisible by step)
		if n % step != 0:
			a = l[step * (n // step - 1):step * (n // step)]		# left block
			b = l[step * (n // step):n]		# extra block
			j = 0		# iterator
						
			for i in xrange(len(a)):
				while j < len(b) and not comp(a[i],b[j]):
					l[step * (n // step - 1) + i + j] = b[j]
			
					j += 1
				#endwhile
		
				l[step * (n // step - 1) + i + j] = a[i]
			#endfor
	
			while j < len(b):	
				l[step * (n // step - 1) + len(a) + j] = b[j]
		
				j += 1
			#endwhile
		#endif
		
		step = step << 1
	#endfor
	
	return l
#enddef

# sorts a list using comp to compare items
# l - list of A, comp - boolean(A, A)
# return - list of A
def sort(l, comp):
	return merge_sort_iter(l, comp)
#enddef

# insert a into sorted l using a binary insertion sort
# l - list of A, a - A, start - int, end -int, comp - boolean(A, A)
def insert(l, a, start, end, comp):
	s = start		# left bound
	e = max(len(l), end)		# right bound
	p = (s + e) // 2		# midpoint
	
	while e - s > 0 and 0 <= p and p < len(l):	
		if comp(l[p], a):
			s = p + (s + e) % 2
		else:
			e = p
		#endif
		
		p = (s + e) // 2
	#endwhile
	
	l.insert(p, a)
#enddef
