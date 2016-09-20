# ANGES 1.01, reconstructing ANcestral GEnomeS maps
# July 2012.
# Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca

import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))) + '/COMMON')

import bm
from acs import *

#######################################################
#    acs_nonunique.py
#
#    Contains code to compute common intervals and
#    adjacencies of sequences (non-unique markers).
#    
#######################################################

## ----------------------------------------------------------
## Common intervals with non-unique markers (sequences)
## ----------------------------------------------------------
def join_chromosomes(blocks, spe, markers,circ):
        extended_markers = copy.copy(markers)
        seq = []
        fake = max(markers) + 1
        extended_markers.append(fake)
        for ch in blocks[spe]:
                if(len(ch) > 0):
                       seq += [m._id for m in ch]
                       if circ:
                           seq += [m._id for m in ch]                
                       seq.append(fake)
        #endfor
        return extended_markers, seq
#enddef

# sigma is the total list of markers, not necessarily present in S1

# Define data structures 

def POSandNUM(sigma, S1, f):
    # POS is a dictionary that stores the positions of different markers
    POS=dict([[x, []] for x in sigma])
    i=0
    while i<len(S1):
        POS[S1[i]].append(i)
        i=i+1
    # NUM is a list that stores the number of distinct markers in an interval of S1
    NUM=[[1 for x in xrange(i,len(S1))] for i in xrange(0,len(S1))]
    i=0
    while i<len(S1)-1:
    	while S1[i]==f:
    		i+=1
        list_of_chars=set([S1[i]])
        j=1
        while j<len(NUM[i]):
            if S1[i+j] not in list_of_chars:
                list_of_chars.add(S1[i+j])
                NUM[i][j]=NUM[i][j]+NUM[i][j-1]
            elif S1[i+j] in list_of_chars:
                NUM[i][j]=NUM[i][j-1]
            j=j+1
        i=i+1
    return [POS,NUM] 

# POS datastructure for finding adjacencies. Does not compute NUM (save quadratic time computation)   
def POS(sigma, S1):
    # POS is a dictionary that stores the positions of different markers
    POS=dict([[x, []] for x in sigma])
    i=0
    while i<len(S1):
        POS[S1[i]].append(i)
        i=i+1
    return POS    

# Function for checking if an interval is left and right maximal    

def fast_maximality(p,intervals):
    intervals[p]=[look_left(p,intervals),look_right(p,intervals)]
    for k in [intervals[p][0],intervals[p][1]]:
   	intervals[k]=intervals[p]
    return intervals[p]


# Function for left maximality

def look_left(p,intervals):
    if p>0 and intervals[p-1]:
        return intervals[p-1][0]
    else:
        return p

# Function for right maximality

def look_right(p,intervals):
    if p<len(intervals)-1 and intervals[p+1]:
        return intervals[p+1][1]
    else:
        return p

# Function to find and store common intervals between two sequences S1 and S2

def common_intervals_sequences(sigma,S1,S2,POS,NUM,MAX,f,is_circ):
    i=0
    acs = bm.BinaryMatrix()
    bounds=[]    
    # List of known common intervals to prevent duplicates
    list_of_common_intervals=[]
    while i<len(S2)-1:
	if S2[i]<f: # If we are not on a fake marker separating chromosomes	    
		# OCC is a binary vector that notes the occurence of a marker in S2
		OCC=dict([[x,0] for x in sigma])
		j=i
    
		# intervals is an array that stores the left and right bounds of the interval containing that position.  
		intervals=[[] for x in S1]
	        
	        if not(is_circ):
			# While j<=|S2| and (i,j) is left-maximal in S2(just make sure that S2[i-1] is not encountered)
	#		while (i>0 and j<len(S2) and S2[i-1]!=S2[j]) or (i==0 and j<len(S2)):
			while (i>0 and S2[j]<f and S2[i-1]!=S2[j]) or (i==0 and S2[j]<f):
				char=S2[j]
				OCC[char]=1
	    
			# While (i,j) is not right-maximal, increase j till it becomes right maximal.
#				while j<len(S2)-1 and OCC[int(S2[j+1])]==1:
				while S2[j]<f and OCC[int(S2[j+1])]==1:
					j=j+1
				for p in POS[char]:
					[start,end]=fast_maximality(p,intervals)

					if end-start > 0 and NUM[start][end-start]==sum([OCC[x] for x in OCC]) and NUM[start][end-start]>1:
		
						if (MAX==0):
							if not(set(S1[start:end+1]) in list_of_common_intervals):
								acs.add_row(set(S1[start:end+1]))
								list_of_common_intervals.append(set(S1[start:end+1]))
						else:
							if (start!=0 or end!=len(S1)-2): # -2 to account for the specific fake marker at the end of 1
								bounds.append([start,end])
								bounds.append([end,start])
				j=j+1
		elif (is_circ):
			
			while (i>0 and S2[j]<f and S2[i-1]!=S2[j] and abs(i-j)<(len(S2)-1)/2 and i<(len(S2)-1)/2) or (i==0 and S2[j]<f and abs(i-j)<(len(S2)-1)/2):
				char=S2[j]
				OCC[char]=1
	    			
			# While (i,j) is not right-maximal, increase j till it becomes right maximal.
#				while j<len(S2)-1 and OCC[int(S2[j+1])]==1:
				while S2[j]<f and OCC[int(S2[j+1])]==1 and abs(i-j)<(len(S2)-1)/2:
					j=j+1

				for p in POS[char]:
					[start,end]=fast_maximality(p,intervals)
					if abs(start-end)>=(len(S1)-1)/2:
					    end=start+(len(S1)-1)/2-1
					check_length=(abs(start-end)<(len(S1)-1)/2 and start<(len(S1)-1)/2)
					    
					# Make sure that the interval covers all markers in S2[i:j+1]
					if check_length and end-start >0 and NUM[start][end-start]==sum([OCC[x] for x in OCC]) and NUM[start][end-start]>1:
						if (MAX==0):
							if not(set(S1[start:end+1]) in list_of_common_intervals):
								acs.add_row(set(S1[start:end+1]))
								list_of_common_intervals.append(set(S1[start:end+1]))
						else:
							bounds.append([start,end])
							bounds.append([end,start])
				j=j+1
	if is_circ and i<(len(S2)-1)/2:
		i+=1
	elif is_circ and i>=(len(S2)-1)/2:
		break
	elif not(is_circ):
		i+=1
    
    if (MAX==0):
	    return acs
    
    elif (MAX==1):
    	bounds.sort()
	# Sorting bounds in increasing order of first index, and subsorting in decreasing order of second index.
	for i in xrange(len(bounds)-1):
	    if bounds[i][0]>bounds[i][1]:
	        j=i+1
		while j<len(bounds) and bounds[j][0]==bounds[i][0] and bounds[j][0]>bounds[j][1]:
		    j=j+1
		    bounds[i:j]=sorted(bounds[i:j],reverse=True,key=lambda b: b[-1])
		i=j
	#Find maximal intervals in S1 using bounds list.
	count=0
	for i in xrange(len(bounds)):
	    if bounds[i][0]<bounds[i][1]:
	        count=count+1
	    elif bounds[i][0]>bounds[i][1]:
		if count==1:
		    acs.add_row(set(S1[bounds[i][1]:bounds[i][0]+1]))
	        count=count-1
        
    elif (MAX==2):
        bounds.sort()
	# Sorting bounds in increasing order of first index, and subsorting in decreasing order of second index.
	for i in xrange(len(bounds)-1):
	    j=i+1
	    while j<len(bounds) and bounds[j][0]==bounds[i][0]:
	        j=j+1
	    bounds[i:j]=sorted(bounds[i:j],reverse=True,key=lambda b: b[-1])
	    i=j

	#Find strong intervals in S1 using bounds list.
        check_intervals=[]
	check_interval_content=[]
	# For all intervals, make sure that no interval strictly overlaps it. If not, add to list of strong intervals
	for i in xrange(len(bounds)):
	    counter=0
	    if bounds[i][0]<bounds[i][1] and not(bounds[i] in check_intervals):
	        j=i+1
		while j<len(bounds) and bounds[j][0]<=bounds[i][1]:
		    if sorted(bounds[j])!=bounds[i] and bounds[j][0]>bounds[j][1]:
		        if bounds[j][1]<bounds[i][0] and bounds[j][0]!=bounds[i][1]:
		    	    counter+=1
		    	    break
		    elif bounds[j][0]<bounds[j][1] and bounds[j]!=bounds[i]:
			if bounds[j][1]>bounds[i][1] and bounds[j][0]!=bounds[i][0]:
			    counter+=1
			    break
		    j+=1
		if counter==0 and not(set(S1[bounds[i][0]:bounds[i][1]+1]) in check_interval_content):
		    check_interval_content.append(set(S1[bounds[i][0]:bounds[i][1]+1]))
		    acs.add_row(set(S1[bounds[i][0]:bounds[i][1]+1]))
		check_intervals.append(bounds[i])
	if [0,len(S1)-2] in bounds and len(set(S1[0:(len(S1)-2)]))>1:
		acs.add_row(set(S1[0:len(S1)-2]))	
	
	if is_circ:
		if [0,(len(S1)-1)/2-1] in bounds and (len(S1)-1)/2-1!=0 and len(check_interval_content)==0 and len(set(S1[0:(len(S1)-1)/2]))>1:
			acs.add_row(set(S1[0:(len(S1)-1)/2]))
    return acs

# Function for computing supported adjacencies given two sequences
    
def supported_adjacencies_sequences(sigma,S1,S2,POS,f,is_circ=False):
    i=0
    acs = bm.BinaryMatrix()
    list_of_adjacencies=[]

    while i<len(S2)-1:
    
        # OCC is a binary vector that notes the occurence of a marker in S2
        OCC=dict([[x,0] for x in sigma])
        j=i
        char1=int(S2[i])
	# While j<=|S2|, move j till exactly two markers are encountered, i.e. an adjacency.
        while (i>0 and j<len(S2) and sum([OCC[x] for x in OCC])<2) or (i==0 and j<len(S2) and sum([OCC[x] for x in OCC])<2):
            char2=int(S2[j])
            OCC[char2]=1
            j=j+1
        # Check for consecutive positions in S1. If circular chromosomes are being handled, check for circular adjacency in S1 (excluding fake marker)    
        for p in POS[char1]:
            for q in POS[char2]:
                
                if not(set([char1,char2]) in list_of_adjacencies):
                    if abs(p-q)==1:
                        list_of_adjacencies.append(set([char1,char2]))
                        acs.add_row(set([char1,char2]))
                    else:
                        if is_circ:
                    	    if abs(p-q)==len(S1)-2:
                                        list_of_adjacencies.append(set([char1,char2]))
                                        acs.add_row(set([char1,char2]))
        # Set pointer to next marker in S2.            
        if j<len(S2) and S2[j]!=f:
        	i=j-1
        else:
        	i=j+1
        
    # Final check: If circular chromosomes are handled, check if the first and last marker of S2 chromosome (excluding fake markers) form adjacencies in S1.
    if is_circ:
        char1=S2[0]
        char2=S2[-2] # Account for fake marker
	if char1!=char2 and not(set([char1,char2]) in list_of_adjacencies):
            for p in POS[char1]:
                for q in POS[char2]:
                    if abs(p-q)==1:
                        list_of_adjacencies.append(set([char1,char2]))
                        acs.add_row(set([char1,char2]))
                    elif abs(p-q)==len(S1)-2:
                    	list_of_adjacencies.append(set([char1,char2]))
                        acs.add_row(set([char1,char2]))
    
    return acs
