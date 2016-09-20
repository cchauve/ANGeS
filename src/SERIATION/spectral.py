# ANGES 1.01, reconstructing ANcestral GEnomeS maps
# July 2012.
# Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca

import sys
import numpy
import pqtrees

sys.setrecursionlimit(10000)

# Seriation algorithm
def seriation(A,column_names,recursive_level):
    n=len(A)

    # List of indices

    indices=sorted(column_names.values())
    if n>1:   
	C=[]
        small=min([A[i,j] for i in xrange(0,n-1) for j in xrange(i+1,n)])
        A=A-small*numpy.ones(n)
	i=0
	# Check for reducibility
        while(len(indices)>0):
            Ci=[]
            components(indices[i],indices,Ci,A)
	    Ci=numpy.sort(Ci)
	    C.append(Ci)

	# If reducible, create P node, and call seriation() on each subcomponent
        if len(C)>1:
            if recursive_level==0:
                C.sort(key= lambda component: component[0])
	        T=[]#pqtrees.Tree(column_names.keys(),'P')
	        for Ci in C:
	    	    subset_columns={}
	    	    i=0
	    	    for index in Ci:
	    	        for pairs in column_names.items():
	    		    if index==pairs[1]:
	    		    	subset_columns[pairs[0]]=i
	    		    	i+=1
		    T.append(seriation(A[Ci,:][:,Ci],subset_columns,recursive_level+1))
	    else:
	        C.sort(key= lambda component: component[0])
	        T=pqtrees.Tree(column_names.keys(),'P')
	        for Ci in C:
	    	    subset_columns={}
	    	    i=0
	    	    for index in Ci:
	    		for pairs in column_names.items():
	    		    if index==pairs[1]:
	    		    	subset_columns[pairs[0]]=i
	    		    	i+=1
		    T.children.append(seriation(A[Ci,:][:,Ci],subset_columns,recursive_level+1))	
	    return T
         
	# If irreducible with smallest off diagonal entry 0, use eigenvector sorting
        elif len(C)==1:

            # Base case: Only two columns. Set as P node. 
            if n==2:
	        subset_columns={}
	        i=0
	        for index in C[0]:
	            for pairs in columns_names.items():
	                if index==pairs[1]:
	    		    subset_columns[pairs[0]]=i
	    		    i+=1
	        T=pqtrees.Tree(subset_columns.keys(),'Q')
	        for c in column_names.items():
		    T.children.append(pqtrees.Tree(c,"leaf"))
		if recursive_level==0:
		    recursive_level+=1
		    T=[T]
 	        return T
	    
 	    # Eigenvector computation
	    elif n>=3:
	    	# Find Laplacian
	        L=numpy.diag(sum(A[:,:]),0)-A
	        evalues,evectors=numpy.linalg.eig(L)

	        # Sort eigenvalues
	        eigsort=numpy.argsort(evalues)

		# For each eigenvector associated to the second smallest eigenvalue, store the permutation that sort the eigenvector
#		for m in multiples:
     	        fsort=numpy.argsort(evectors[:,eigsort[1]].T)
                CC=[]    
                T=pqtrees.Tree(column_names.keys(),'Q')
                i=0
		while i<n:
		   
	           Ci=[fsort[i]]
		   j=i+1
	           # Group like values in an eigenvector. Threshold of .000001 again,
	           while (j<n):
                       if(numpy.absolute(evectors[fsort[i],eigsort[1]]-evectors[fsort[j],eigsort[1]])<.000001):
		           Ci.append(fsort[j])
		           j=j+1
		       elif(numpy.absolute(evectors[fsort[i],eigsort[1]]-evectors[fsort[j],eigsort[1]])>.000001):
		           i=j
			   break
                   i=j
		   CC.append(Ci)
    	        # For each group of eigenvector entries, call seriation() again
    	        for Ci in CC:
    	  	    subset_columns={}
    	  	    i=0
	    	    for index in Ci:
	    	        for pairs in column_names.items():
	    		    if index==pairs[1]:
	    		        subset_columns[pairs[0]]=i
		   		i+=1
		    Temp=seriation(A[Ci,:][:,Ci],subset_columns,recursive_level+1)

		    T.children.append(Temp)

		# Check if this is the first recursive level or not.
    	        if recursive_level==0:
    	            recursive_level+=1
    	            T=[T]
    	        return T
    # Base case of single column: Assign as leaf   
    elif n==1:
	T=pqtrees.Tree(column_names.keys(),"leaf")
	if recursive_level==0:
	    recursive_level+=1
	    T=[T]
        return T	

# Finding connected components
def components(i,indices,C,A):
    if not (i in C):
        C.append(i)
    # Check if i-j edge exists (value in matrix is non-zero). If so, add j to C, and call components on j.    
    if len(indices)>0 and i in indices:
        indices.remove(i)
        neighbours=[]
        for  j in indices:
            if A[i][j]>0 and not (j in C):
                neighbours.append(j)
        for j in neighbours:
            components(j,indices,C,A) 
    return None

# Convert Hamming matrix file to numpy array
def correlation_matrix(hamming,index_dict):
    n=len(index_dict)
    A=numpy.zeros([n,n])
    for row in hamming:
    	for col in hamming[row]:
  	    A[index_dict[row],index_dict[col]]=hamming[row][col]
    	    A[index_dict[col],index_dict[row]]=hamming[row][col]
    return A
    
