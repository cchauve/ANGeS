# ANGES 1.01, reconstructing ANcestral GEnomeS maps
# July 2012.
# Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca

class Tree:
    children=[]
    operation=None
    character="placeholder"
    def __init__(self,column_names,op):
	self.operation=op
	if op=="leaf":
	    self.children=[]
	    self.character=column_names[0]
	else:
	    self.children=[]
    def printTree(self):
	if not(self==None):
            if self.children==[]:
	        return str(self.character) 
	    else:
		L= "_"+self.operation+" "
	        for child in self.children:
	            if not(child==None):
		        L=L+str(child.printTree())+" "
		L=L+self.operation+"_"
		return L


def compare_trees(T1,T2):
    if T1.operation==T2.operation:
        if T1.operation=='Q':
	    boolcheck1=True
	    boolcheck2=True
	    i=0
	    while i<len(T1.children):
	        boolcheck1=boolcheck1 and compare_trees(T1.children[i],T2.children[i])
	        boolcheck2=boolcheck2 and compare_trees(T1.children[i],T2.children[len(T1.children)-1-i])
		i=i+1
	    return boolcheck1 or boolcheck2
	elif T1.operation=='P':
	    i=0
	    bit=True
            while i<len(T1.children):
	        bit=bit and compare_trees(T1.children[i],T2.children[i])
	        i=i+1
	    return bit
	elif T1.operation=="leaf":
	    if T1.character==T2.character:
	        return True
	    else:
	        return False    
