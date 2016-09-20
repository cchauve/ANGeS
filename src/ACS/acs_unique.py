# ANGES 1.01, reconstructing ANcestral GEnomeS maps
# July 2012.
# Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca

import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))) + '/COMMON')

import bm
import copy
from acs import *

#######################################################
#    acs.py
#
#    Contains code to compute common intervals and 
#    adjacencies of permutations (unique markers).
#    
#######################################################

## ----------------------------------------------------------
## Common intervals and adjacencies with unique markers
## ----------------------------------------------------------

# map the markers of spe1 to the identity and spe2 to spe1
#
# blocks:         orthology blocks - dict of list of list of Marker
# spe1:           the first species - str
# spe2:           the second species - str
# gen_to_id:      (out) the map from spe1 to the identity - list of int
# id_to_gen:      (out) the map from the indentity to spe2 - list of int
# seq2:           (out) the map from spe2 to spe1 - list of int
# fictive_blocks: the markers representing a change in chromosome - list of int
def genome_to_identity(blocks, spe1, spe2, gen_to_id, id_to_gen, seq2, fictive_blocks):
	i = 0		# identity iterator
	nb_chr1 = len(blocks[spe1])		# number of chromosomes in spe1
	nb_chr2 = len(blocks[spe2])		# number of chromosomes in spe2

#	fictive_blocks[0].append(i)
#	i+=1
	for j in range(nb_chr1):
		nb_blocks = len(blocks[spe1][j])  # number of blocks in the chromosome
		for k in range(nb_blocks):
			gen_to_id[blocks[spe1][j][k]._id] = i
			id_to_gen[i] = blocks[spe1][j][k]._id
			i+=1
		#endfor		
		if(nb_blocks>0):
			fictive_blocks[0].append(i)
			i+=1
		#endif
	#endfor

#	fictive_blocks[1].append(i)
#	seq2.append(i)
#	i+=1
	for j in range(nb_chr2):
		nb_blocks = len(blocks[spe2][j])
		for k in range(nb_blocks):
			seq2.append(gen_to_id[blocks[spe2][j][k]._id])
		#endfor		
		if(nb_blocks>0):
			seq2.append(i)
			fictive_blocks[1].append(i)
			i+=1
		#endif
	#endfor
	for i in fictive_blocks[0]:
		seq2.append(i)
	#endfor
#enddef

# Inner function for compute_common_intervals.
# Uses the algorithm of "Computing common intervals of K permutations, with
# applications to modular decomposition of graphs" by Bergeron et al
#
# seq2:             the map from seq2 to seq1 - list of int
# fictive_blocks:   the markers representing change in chromosome - list of int
# common_intervals: (out) the common_intervals - list of list of int
# output_intervals: kind of common intervals: 0 (STRONG), 1 (ALL), 2 (MAXIMAL)
def compute_common_intervals_inner(seq2, fictive_blocks, common_intervals, output_intervals):
	IMax=[]
	IMin=[]
	Sup=[]
	Inf=[]
	SupportR=[]
	lenseq2=len(seq2)  #-2

	#initializing
	for i in range(lenseq2):
		IMax.append([0,0])
		IMin.append([0,0])
		Sup.append(0)
		Inf.append(0)
		SupportR.append(0)
	#endfor

	for i in range(lenseq2):	
		#computing IMax
		start = i
		end = i		
		while(start >= 0 and seq2[start] >= seq2[i]):
			start-=1
		#endwhile
		start+=1
		while(end <= (lenseq2-1) and seq2[end] >= seq2[i]):
			end+=1
		#endwhile
		end-=1
		IMax[seq2[i]][0]=start
		IMax[seq2[i]][1]=end
		#computing Sup
		imax = seq2[start:(end+1)]
		imax.sort()
		sup = 0
		while (sup < len(imax)-1 and imax[sup+1] == imax[sup]+1):
			sup+=1
		#endwhile
		Sup[seq2[i]]=imax[sup]
		#computing IMin
		start = i
		end = i
		while(start >= 0 and seq2[start] <= seq2[i]):
			start-=1
		#endwhile
		start+=1
		while(end <= (lenseq2-1) and seq2[end] <= seq2[i]):
			end+=1
		#endwhile
		end-=1
		IMin[seq2[i]][0]=start
		IMin[seq2[i]][1]=end
		#computing Inf
		imin = seq2[start:end+1]
		imin.sort()
		inf = len(imin)-1
		while (inf>0 and imin[inf-1] == imin[inf]-1):
			inf-=1
		#endwhile
		Inf[seq2[i]]=imin[inf]
	#endfor	

	#computing support for R
	SupportR[0] = -1
	S=[]
	S.append(0)
	for i in range(1,lenseq2):
		while(Sup[S[-1]] < Sup[i]):
			S.pop()
		#endwhile
		SupportR[i] = S[-1]
		S.append(i)
	#endfor

	# Computing strong common intervals
	if output_intervals==0: 
		SupportL=[]
		CanR=[]     # canonical generator
		CanL=[]     # canonical generator
		LeftBounds=[]  
		RightBounds=[] 

		#initializing
		for i in range(lenseq2):
			SupportL.append(0)
			CanR.append(0)
			CanL.append(0)
	        #endfor

                #computing support for L
		SupportL[lenseq2-1] = -1
		S=[]
		S.append(lenseq2-1)
		for i in range(0,lenseq2-1)[::-1]:
			while(Inf[S[-1]] > Inf[i]):
				S.pop()
		        #endwhile
			SupportL[i] = S[-1]
			S.append(i)
		#endfor

	        #computing canonical generator R
		CanR[0] = lenseq2-1
		for k in range(1,lenseq2):
			CanR[k] = k
	        #endfor
		for k in range(1,lenseq2)[::-1]:
			i = SupportR[k]
			Ri = Sup[i]
			j = CanR[k]
			Lj = Inf[j]
			if(Ri>=j and Lj<=i):
				CanR[i]=max(CanR[k],CanR[i])
		        #endif
                #endfor
				
	        #computing canonical generator L
		CanL[lenseq2-1] = 0
		for k in range(0,lenseq2-1)[::-1]:
			CanL[k] = k
	        #endfor
		for k in range(0,lenseq2-1):
			j = SupportL[k]
			Lj = Inf[j]
			i = CanL[k]
			Ri = Sup[i]
			if(Ri>=j and Lj<=i):
				CanL[j]=min(CanL[k],CanL[j])
		        #endif
	        #endfor
			
		#setting bounds
		for k in range(lenseq2):
			LeftBounds.append(k)
			LeftBounds.append(CanL[k])
			RightBounds.append(k)
			RightBounds.append(CanR[k])
		#endfor
		LeftBounds.sort()
		RightBounds.sort()

		#computing strong intervals
		S=[]
		i=0
		Leftj=0
		Rightj=0
		while i<4*lenseq2:
			if Leftj<2*lenseq2 and LeftBounds[Leftj]<=RightBounds[Rightj]:
				S.append(LeftBounds[Leftj])
				Leftj=Leftj+1
			else:
				common_intervals.append([S[-1],RightBounds[Rightj]])
				S.pop()
				Rightj=Rightj+1
			#endif
			i=i+1
		#endwhile
		#discarding repeated common intervals
		common_intervals.sort()
		i=1
		while (i< len(common_intervals)):
			if(common_intervals[i] == common_intervals[i-1]):
				del(common_intervals[i])
			else:
				i+=1
        	        #endif
        	#endwhile		
        #endif

	#computing all common intervals or maximal ones
	if output_intervals==1 or output_intervals==2:
		j = lenseq2-1	
		while(j >=0):
			i = j
			while(i>=Inf[j] and i!=-1):
				common_intervals.append([i,j])
				i = SupportR[i]
		        #endwhile
			j-=1
	        #endwhile
	#endif

	#discarding size 1 common intervals
	i=0	
	while (i< len(common_intervals)):
		if(common_intervals[i][1] == common_intervals[i][0]):
			del(common_intervals[i])
		else:
			i+=1
	        #endif
	#endwhile

	# #discarding ajacencies 
	#i=0
	#while (i< len(common_intervals)):
	#	if(common_intervals[i][1] == common_intervals[i][0]+1):
	#		del(common_intervals[i])
	#	else:
	#		i+=1
	#	#endif
	#endwhile

	#discarding common intervals containing fictive blocks
	all_fictive_blocks = fictive_blocks[0]+fictive_blocks[1]
	size = len(all_fictive_blocks)
	common_intervals.sort()
	k=0
	deletion = []
	
	for i in range(len(common_intervals)):
		a=common_intervals[i][0]
		b=common_intervals[i][1]
		while(all_fictive_blocks[k]<a):
			k+=1
		#endwhile
		j=k
		deleted = False
		while(j < size and deleted == False):
			if(a<=all_fictive_blocks[j] and all_fictive_blocks[j]<=b):
				deletion.append(i)
				deleted = True
			#endif
			j+=1
		#endwhile
	#endfor
	i=len(deletion)-1
	while(i >= 0):
		del(common_intervals[deletion[i]])
		i-=1
	#endwhile

	#discarding a common interval containing all blocks
	i=0		
	while (i< len(common_intervals)):
		if common_intervals[i][1] - common_intervals[i][0] +1 == lenseq2 - 2*len(fictive_blocks):
			del common_intervals[i]	
		else:
			i+=1
	#endwhile

	# discarding non-maximum common intervals
#	# also discards a common intervals if it contains all the markers
	if output_intervals==2:
		i=0		
		while (i< len(common_intervals)-1):
#			if common_intervals[i][1] - common_intervals[i][0] >= lenseq2 - 3:
#				del common_intervals[i]
#			elif common_intervals[i+1][1] - common_intervals[i+1][0] >= lenseq2 - 3:
#				del common_intervals[i+1]
			if common_intervals[i][0] == common_intervals[i+1][0]:
				common_intervals[i][1] = common_intervals[i+1][1]				
				del(common_intervals[i+1])
			elif common_intervals[i][1] >= common_intervals[i+1][1]:
				common_intervals[i][1] = common_intervals[i][1]#max(common_intervals[i][1], common_intervals[i+1][1])				
				del(common_intervals[i+1])
			else:
				i+=1
			#endif
		#endwhile
	#endif	
			
	#print common_intervals
#enddef

def compute_common_intervals(blocks, spe1, spe2, output_intervals, is_circular, all_common_intervals = bm.BinaryMatrix()):
	gen_to_id={}
	id_to_gen={}
	seq2=[]
	fictive_blocks=[[],[]]
	
	if len(blocks[spe1]) > 0 and len(blocks[spe2]) > 0:
		genome_to_identity(blocks, spe1, spe2, gen_to_id, id_to_gen, seq2, fictive_blocks)
		common_intervals=[]
		compute_common_intervals_inner(seq2, fictive_blocks, common_intervals, output_intervals)
	else:
		common_intervals=[]
	#endif

	for ci in common_intervals:
		interval = []
		
		for i in xrange(ci[0],ci[1]+1):
			interval.append(id_to_gen[i])
		#endfor
		
		all_common_intervals.add_row(set(interval))
		all_common_intervals.get_row_info(-1)._sp = [spe1, spe2]
		
		# add complement
		if is_circular:
			interval2 = []
		
			for i in xrange(ci[0]):
				interval2.append(id_to_gen[i])
			#endfor
			
			for i in xrange(ci[1] + 1, len(seq2) - len(fictive_blocks)):
				interval2.append(id_to_gen[i])
			#endfor
			
			all_common_intervals.add_row(set(interval2))
			all_common_intervals.get_row_info(-1)._sp = [spe1, spe2]
		#endif
	#endfor
	
	return all_common_intervals
#enddef

def compute_adjacencies_inner(seq2, fictive_blocks, is_circular, adjacencies):
	# find adjacencies	
	lenseq2=len(seq2)  #-2
	for i in xrange(lenseq2 - 1):
		if abs(seq2[i] - seq2[i + 1]) == 1 or (is_circular and abs(seq2[i] - seq2[i + 1]) == lenseq2 - 3):
			adjacencies.append([i, i + 1])
		#endif
	#endfor
	
	if is_circular:
		if abs(seq2[-3] -seq2[0]) == 1 or abs(seq2[-3] - seq2[0]) == lenseq2 - 3:
			adjacencies.append([0, lenseq2 - 3])
		#endif
	#endif
	
	#discarding common intervals containing fictive blocks
	all_fictive_blocks = fictive_blocks[0]+fictive_blocks[1]
	size = len(all_fictive_blocks)
	adjacencies.sort()
	k=0
	deletion = []	
		
	for i in range(len(adjacencies)):	
		a=adjacencies[i][0]
		b=adjacencies[i][1]
		
		while(all_fictive_blocks[k]<a):
			k+=1
		#endwhile		
		j=k
		deleted = False		
		while(j < size and deleted == False):		
			if(seq2[a]==all_fictive_blocks[j] or all_fictive_blocks[j]==seq2[b]):
				deletion.append(i)
				deleted = True
			#endif		
			j+=1
		#endwhile
	#endfor	
	deletion.sort()	
	i=len(deletion)-1	
	while(i >= 0):
		del(adjacencies[deletion[i]])
		i-=1
	#endwhile
	
	return adjacencies
#enddef

def compute_adjacencies(blocks, spe1, spe2, is_circular = False, all_adjacencies = bm.BinaryMatrix()):
	gen_to_id={}
	id_to_gen={}
	seq2=[]
	fictive_blocks=[[],[]]
	
	if len(blocks[spe1]) > 0 and len(blocks[spe2]) > 0:
		genome_to_identity(blocks, spe1, spe2, gen_to_id, id_to_gen, seq2, fictive_blocks)
		adjacencies=[]
		compute_adjacencies_inner(seq2, fictive_blocks, is_circular, adjacencies)
	else:
		adjacencies=[]
	#endif

	for ci in adjacencies:
		interval = [id_to_gen[seq2[ci[0]]], id_to_gen[seq2[ci[1]]]]
		all_adjacencies.add_row(set(interval))
		all_adjacencies.get_row_info(-1)._sp = [spe1, spe2]
	#endfor	
	return all_adjacencies
#enddef

# uv is a reliable adjacency if (1) there exist xy st. xy and uv are realized in
# spe1, xy is realized in spe2 and ux and vy are realized in spe3; and (2) the
# ancestor is on the path between spe1 and spe2 and the path between spe1 and spe3.
# From this it is necessary that xy is a supported adjacency.
#
# If uv is as above then uv is a reliable adjacency of the median of spe2 and spe3, which is an internal node.
#
# This algorithm may output adjacencies that are also supported adjacencies.
#
# blocks: filtered orthology blocks - dict of list of list of Marker
# supported_adjacencies: supported adjacencies of spe1 and spe2 - bm.BinaryMatrix
# spe1: name of species1 - str
# spe3: name of species3 - str
# all_reliable_adjacencies: (out) the set of reliable adjacencies - bm.BinaryMatrix
def compute_reliable_adjacencies(blocks, supported_adjacencies, spe1, spe3, is_circular, all_reliable_adjacencies = bm.BinaryMatrix()):
	gen_to_id={}
	id_to_gen={}
	seq3=[]
	fictive_blocks=[[],[]]

	genome_to_identity(blocks, spe1, spe3, gen_to_id, id_to_gen, seq3, fictive_blocks)
	
	fb = fictive_blocks[0] + fictive_blocks[1]
	
	id_to_id3={}
	
	for i in xrange(len(seq3)):
		id_to_id3[seq3[i]] = i
	#endfor
			
	for adj in supported_adjacencies:
		adj_row = [x for x in adj._set]
	
		# don't try induced adjacencies from doubling
		if (adj_row[0] == adj_row[1] + 1 and adj_row[0] % 2 == 0) or (adj_row[1] == adj_row[0] + 1 and adj_row[1] % 2 == 0):
			continue
		#endif
	
		if adj_row[0] in gen_to_id and adj_row[1] in gen_to_id:
			reliable_adjacencies = []
		
			x = gen_to_id[adj_row[0]]
			y = gen_to_id[adj_row[1]]
		
			# u's and v's in spe3
			if is_circular:
				u1 = seq3[(id_to_id3[x] - 1) % (len(seq3) - 2)]
				u2 = seq3[(id_to_id3[x] + 1) % (len(seq3) - 2)]
				v1 = seq3[(id_to_id3[y] - 1) % (len(seq3) - 2)]
				v2 = seq3[(id_to_id3[y] + 1) % (len(seq3) - 2)]
			else:
				u1 = seq3[id_to_id3[x] - 1]
				u2 = seq3[id_to_id3[x] + 1]
				v1 = seq3[id_to_id3[y] - 1]
				v2 = seq3[id_to_id3[y] + 1]
			#endif
			
			# check if uv is an adjacency in spe1
			if abs(u1 - v1) == 1 or (is_circular and abs(u1-v1) == len(seq3)-3):
				reliable_adjacencies.append([u1, v1, spe1, x, y])
			#endif
			
			if abs(u1 - v2) == 1 or (is_circular and abs(u1-v2) == len(seq3)-3):
				reliable_adjacencies.append([u1, v2, spe1, x, y])
			#endif
			
			if abs(u2 - v1) == 1 or (is_circular and abs(u2-v1) == len(seq3)-3):
				reliable_adjacencies.append([u2, v1, spe1, x, y])
			#endif
			
			if abs(u2 - v2) == 1 or (is_circular and abs(u2-v2) == len(seq3)-3):
				reliable_adjacencies.append([u2, v2, spe1, x, y])
			#endif
			
			for ra in reliable_adjacencies:
				if ra[0] not in fb and ra[1] not in fb:
					a = id_to_gen[ra[0]]
					b = id_to_gen[ra[1]]
			
					# don't include induced adjacencies from doubling
					if (a == b + 1 and a % 2 == 0) or (b == a + 1 and b % 2 == 0):
						continue
					#endif			
					
					all_reliable_adjacencies.add_row(set([a, b]))
					all_reliable_adjacencies.get_row_info(-1)._sp = [ra[2]]
				#endif
			#endfor
		#endif
	#endfor
	
	return all_reliable_adjacencies
#enddef
