# ANGES 1.01, reconstructing ANcestral GEnomeS maps
# July 2012.
# Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca

import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))) + '/COMMON')

import bm

import copy
import gc

#######################################################
#    acs.py
#
#    Contains code to read orthology blocks files and
#    manipulate ACS files: adding missing markers and 
#    telomeres
#    
#######################################################

## ----------------------------------------------------------
## Markers
## ----------------------------------------------------------

# defines a marker
class Marker:
	def __init__(self, ident, species, ch, st, end):
		self._id = ident
		self._species = species
		self._ch = ch
		self._st = st
		self._end = end
	#enddef
	
	def __lt__(self, a):
		return self._st < a._st
	#enddef
	
	def __le__(self, a):
		return self._st <= a._st
	#enddef
#endclass

# extracts the orthology blocks from an orthology blocks file
#
# blocks_file: orthology blocks file name - str
#
# return: orthology blocks - dict of list of list of Marker
def extract_blocks(blocks_file):
	blocks = {}	# orthology blocks
	blocksfile = open(blocks_file,"r").readlines()		# the file of the orthology blocks
	chm={}		# the chromosomes in the blocks
	markers=[]
	i=0		# file index
	id = ""		# current marker id
	
	while i<len(blocksfile):
		if blocksfile[i][0]=='>':
			id=blocksfile[i][1:-1].split(' ')[0]
			i=i+1
			markers.append(int(id))
			while i<len(blocksfile) and len(blocksfile[i])>1:
				spe=blocksfile[i][:blocksfile[i].find(".")]		# current species
				chr=blocksfile[i][blocksfile[i].find(".")+1:blocksfile[i].find(":")]		# current chromosome
				if not spe in blocks:
					blocks[spe] = []
				#endif
				start=int(blocksfile[i][blocksfile[i].find(":")+1:blocksfile[i].find("-")])		# start pos
				stop=int(blocksfile[i][blocksfile[i].find("-")+1:blocksfile[i].find(" ")])		# end pos
				if(chm.has_key(chr.upper()) == False):
					if chm=={}:
						chm[chr.upper()] = 0
					else:
						chm[chr.upper()] = max(chm.values())+1
				#endif
				if (chm[chr.upper()] >= len(blocks[spe])):
					for k in range(chm[chr.upper()]-len(blocks[spe])+1):
						blocks[spe].append([])
					#endfor
				#endif
				blocks[spe][chm[chr.upper()]].append(Marker(int(id), spe, chr.upper(), int(start), int(stop)))
				i=i+1
			#endwhile
			i=i+1
		else:
			i+=1
		#endif
	#endwhile

	for i in blocks:
		for j in range(len(blocks[i])):
			blocks[i][j].sort()
		#endfor
	#endfor
	return markers,blocks
#enddef

# strips orthology blocks so that each given species has the same set of markers
#
# blocks: orthology blocks - dict of list of list of Marker
# spes: species names - list of str
#
# return: filitered orthology blocks - dict of list of list of Marker
def filter_markers(blocks, spes):
	markers = {}		# flags markers that are used
	
	# find markers used in first species
	for ch in blocks[spes[0]]:	
		for m in ch:
			markers[m._id] = False
		#endfor
	#endfor
	# find markers markers used in the rest of the species
	for sp in xrange(1, len(spes)):
		for ch in blocks[spes[sp]]:
			for m in ch:
				if m._id in markers:
					markers[m._id] = True
				#endfor
			#endfor
		#endfor
		delete = []
		for m in markers:
			if markers[m]:
				markers[m] = False
			else:
				delete.append(m)
			#endif
		#endfor
		for d in delete:
			del markers[d]
		#endfor
	#endfor
	filtered = {}	# the filtered orthology blocks	
	# filter the orthology blocks
	for sp in spes:
		filtered[sp] = []
		for ch in blocks[sp]:
			filtered[sp].append([])			
			for m in ch:
				if m._id in markers:
					filtered[sp][-1].append(m)
				#endif
			#endfor
		#endfor
	#endfor	
	return filtered
#enddef

def common_markers(blocks, spes):
	markers = {}    # flags markers that are used
	common_markers=[]  # list of markers
 # find markers used in first species
	for ch in blocks[spes[0]]:  
		for m in ch:
			markers[m._id] = False
			common_markers.append(m._id)
                #endfor
        #endfor
	return common_markers
#enddef

## ----------------------------------------------------------
## Telomeres
## ----------------------------------------------------------

# adds telomere rows to ancestral consecutive syntenies
# modifies ACS
#
# blocks: (filtered?) orthology blocks -  dict of list of list of Marker
# spe1: name of first species - str
# spe2: name of second species - str
# ACS: (out) ancestral consecutive syntenies - bm.BinaryMatrix
def add_telomeres(blocks, spe1, spe2, ACS):
	markers_to_ch = {spe1: {}, spe2: {}}		# map from marker to its chromosome

	for ch in blocks[spe1]:
		for m in ch:
			markers_to_ch[spe1][m._id] = ch
		#endif
	#endfor
	
	for ch in blocks[spe2]:
		for m in ch:
			markers_to_ch[spe2][m._id] = ch
		#endif
	#endfor
	
	telomere_rows = []		# syntenies with telomeres
	
	for row in ACS:
		s = row._set		# row's set of 1's
		x = iter(s).next()		# arbitrary marker from s
		ch1 = markers_to_ch[spe1][x]		# marker's chromosome in spe1
		ch2 = markers_to_ch[spe2][x]		# marker's chromosome in spe2
		
		if (ch1[0]._id in s or ch1[len(ch1)-1]._id in s) and (ch2[0]._id in s or ch2[len(ch2)-1]._id in s):
			telomere_rows.append(bm.Row(s, row._id + "T", 0, [spe1, spe2], True))
		#endif
	#endfor
	
	for t in telomere_rows:
		ACS.add_row_info(t)
	#endfor
#enddef

## ----------------------------------------------------------
## Missing markers
## ----------------------------------------------------------

def add_missing_markers_X(blocks, markers, spe1, spe2, ACS, addForOneSpecies = True):
	markers_inclusion1 = dict([[m, True] for m in markers])	 # inclusion map of missing markers (True id marker is missing)
	markers_inclusion2 = dict([[m, True] for m in markers])	 # inclusion map of missing markers (True id marker is missing)
	markers_inclusion  = dict([[m, True] for m in markers])	 # inclusion map of missing markers (True id marker is missing)
	markers_missing = set([])		# set of missing markers
	
	# addForOneSpecies - True: set markers_inclusion to True for every marker in spe2 (markers_inclusion = (spe2 - spe1)^C)
	# addForOneSpecies - False: set markers_inclusion to False for every marker in spe1 (markers_inclusion = (spe1 | spe2)^C)

	# set markers_inclusion to False for every marker in spe1 (markers_inclusion = spe1^C)
	for ch in blocks[spe1]:	
		for m in ch:
			markers_inclusion1[m._id] = False
		#endfor
	#endfor
	
	for ch in blocks[spe2]:
		for m in ch:
			markers_inclusion2[m._id] = False
                #endfor
        #endfor

	if addForOneSpecies:
		for m in markers_inclusion.keys():
			if markers_inclusion1[m] == False and  markers_inclusion2[m] == False:
				markers_inclusion[m] = False
	if not  addForOneSpecies:
		for m in markers_inclusion.keys():
			if markers_inclusion1[m] == False or markers_inclusion2[m] == False:
				markers_inclusion[m] = False
		
	
#	# addForOneSpecies - True: negate markers_inclusion for every marker in spe1 (markers_inclusion = (spe1 & spe2)^C)
#	if addForOneSpecies:
#		for ch in blocks[spe1]:	
#			for m in ch:
#				markers_inclusion[m._id] = not markers_inclusion[m._id]
#			#endfor
#		#endfor
#	#endif
	
	for m in markers_inclusion:
		if markers_inclusion[m]:
			markers_missing.add(m)
		#endif
	#endfor
	
	ACS_missing = bm.BinaryMatrix()
	
	for syn in ACS:
		if len(syn._set) == 2:
			it = iter(syn._set)
			a = it.next()
			b = it.next()
			
			if a == b + 1 and a % 2 == 0 or b == a + 1 and b % 2 == 0:
				ACS_missing.add_row_info(bm.XRow(syn._set, set([]), syn._id, syn._weight, syn._sp, syn._isT))
				
				continue
			#endif
		#endif
	
		ACS_missing.add_row_info(bm.XRow(syn._set, markers_missing, syn._id, syn._weight, syn._sp, syn._isT))
	#endfor
	
	return ACS_missing
#enddef

def add_missing_markers_1(blocks, spe1, spe2, ACS):
	for syn in ACS:
		if len(syn._set) == 2:
			it = iter(syn._set)
			a = it.next()
			b = it.next()
			
			if a == b + 1 and a % 2 == 0 or b == a + 1 and b % 2 == 0:
				continue
			#endif
		#endif
			
		s = copy.copy(syn._set)
		t = syn._set
		
		found = False
		left = copy.copy(s)
		
		for ch in blocks[spe1]:	
			for m in ch:
				if m._id in syn._set:
					found = True
					left.remove(m._id)
				elif found:
					t.add(m._id)
				#endif
				
				if len(left) == 0:
					break
				#endif
			#endfor
			
			if found:
				break
			#endif
		#endfor
		
		found = False
		left = copy.copy(s)
		
		for ch in blocks[spe2]:	
			for m in ch:
				if m._id in syn._set:
					found = True
					left.remove(m._id)
				elif found:
					t.add(m._id)
				#endif
				
				if len(left) == 0:
					break
				#endif
			#endfor
			
			if found:
				break
			#endif
		#endfor
	#endfor
	
	return ACS
#enddef

def remove_nonmaximal_intervals(ACS):
	ACS.sort_size()
	
	for i in xrange(ACS._height - 1, -1, -1):	
		for j in xrange(ACS._height - 1, i, -1):
			if ACS.get_row(i) < ACS.get_row(j):
				ACS.remove_row(i)			
				
				break
			#endif
		#endfor
	#endfor
#enddef

## ----------------------------------------------------------
## Misc
## ----------------------------------------------------------

# joins ancestral consecutive symntenies together, combining rows with the same  1's
#
# ACSs: ancestral consecutive syntenies to combine - list of bm.BinaryMatrix
#
# return: combined ancestral consecutive syntenies - bm.BinaryMatrix
def join_ACS(ACSs):
	ACS_sort = []
		
	for acs in ACSs:
		for syn in acs:
			if syn._isX == True:# is bm.XRow:
				S = [x for x in syn._set]
				S.sort()
				val = copy.copy(S)
				val.append('X')
				X = [x for x in syn._Xs]
				X.sort()
				val += X
				val.append(syn._isT)
				val += syn._sp
				val.append(syn._weight)
				val.append(len(X))
				val.append(len(S))
			else:
				S = [x for x in syn._set]
				S.sort()
				val = copy.copy(S)
				val.append(syn._isT)
				val += syn._sp
				val.append(syn._weight)
				val.append(len(S))
			#endif
			ACS_sort.append(val)
		#endfor
		del acs
	#endfor
	
	del ACSs
	ACS_sort.sort()
	joined_ACS = bm.BinaryMatrix()
	r = 0
	while len(ACS_sort) > 0:
		i = 0
		current_row = ACS_sort[i]
		if 'X' in current_row:
			isX = True
			a = current_row[-1]
			b = current_row[-2]
			row = bm.XRow(set(current_row[0:a]), set(current_row[a+1:a+b+1]), str(r), current_row[-3], set(current_row[a + b + 2:-3]), current_row[a + b + 1])
		else:
			isX = False
			row = bm.Row(set(current_row[0:current_row[-1]]), str(r), current_row[-3], set(current_row[current_row[-1] + 1:-2]), current_row[current_row[-1]])
		#endif
		while i < len(ACS_sort):
			if isX:
				if a == ACS_sort[i][-1] and b == ACS_sort[i][-2]:
					if current_row[0:a+b+1] == ACS_sort[i][0:a+b+1] and current_row[a+b+1] == ACS_sort[i][a+b+1]:
						row._sp |= set(ACS_sort[i][a + b + 2:-3])
						row._weight = max(current_row[-3], ACS_sort[i][-3])
					else:
						break
					#endif
				else:
					break
				#endif
			else:
				if current_row[-1] == ACS_sort[i][-1]:
					if current_row[0:current_row[-1]] == ACS_sort[i][0:current_row[-1]] and current_row[current_row[-1]] == ACS_sort[i][current_row[-1]]:
						row._sp |= set(ACS_sort[i][current_row[-1] + 1:-2])
						row._weight = max(current_row[-2], ACS_sort[i][-2])
					else:
						break
					#endif
				else:
					break
				#endif
			#endif
			i = i + 1
		#endwhile
		for j in xrange(i):
	 		del	ACS_sort[0]
	 	#endif
		row._sp = [x for x in row._sp]
		row._sp.sort()
		joined_ACS.add_row_info(row)
		r += 1
	#endwhile
	return joined_ACS
#enddef
