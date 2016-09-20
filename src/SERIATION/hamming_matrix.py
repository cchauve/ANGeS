# ANGES 1.01, reconstructing ANcestral GEnomeS maps
# July 2012.
# Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca

# Computes hamming matrix, given a parameter alpha, an XC1P matrix, and the set of indices
import numpy
import bm

def hamming_matrix_numpy(m, alpha):
    support={}
    i=0
#    print m.get_support()
    # Support is a dictionary using the column names as keys for the matrix indices.
    for col in m.get_support():
        support[col]=i
        i+=1 
    
    # List of indices for the matrix
    index_list=numpy.array(sorted(support.values()))
    print index_list
    A=numpy.zeros([len(index_list),len(index_list)])
    print support
    
    # Populate symmetric matrix
    for col in m.get_support():
	for row in m:
		if col not in row._set and col not in row._Xs:
			for second_col in row._Xs:
				A[support[col],support[second_col]]+=alpha
				A[support[second_col],support[col]]+=alpha
			for third_col in row._set:
				A[support[col],support[third_col]]+=1
				A[support[third_col],support[col]]+=1

    return [A,support,index_list]


# Returns a symmetric hamming distance matrix encoded as a dictionary of non-zero marker-marker distance

def hamming_matrix(m, alpha):
    support={}
    i=0
    
    # Support is a dictionary using the column names as keys for the matrix indices.
    for col in m.get_support():
        support[col]=i
        i+=1 
    
    # List of indices for the matrix
    index_list=numpy.array(sorted(support.values()))
    
    full_list=dict([(x,dict([])) for x in m.get_support()])    
    # Populate symmetric matrix
    for col1 in full_list:
	
	for col2 in m.get_support():
		dist=0
		if col1 not in full_list[col2].keys():
		    for row in m:
		        if col1 not in row._set and col1 not in row._Xs:
		            if col2 in row._Xs:
		                dist+=alpha
 		            elif col2 in row._set:
		                dist+=1
		        elif col1 in row._set:
		            if col2 not in row._set and col2 not in row._Xs:
		                dist+=1
		        elif col1 in row._Xs:
		            if col2 not in row._set and col2 not in row._Xs:
		                dist+=alpha
		    if dist!=0:
		        full_list[col1][col2]=dist
		        full_list[col2][col1]=dist        
			
    return full_list

def write_to_file(matrix,f):
    
    for row in matrix:
    	f.write(str(row)+':')
    	for col in matrix[row]:
    	    f.write(str(col)+' '+str(matrix[row][col])+';')
    	f.write('\n')
    return None

