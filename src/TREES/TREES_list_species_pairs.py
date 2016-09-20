# ANGES 1.01, reconstructing ANcestral GEnomeS maps
# July 2012.
# Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca

#!/usr/bin/python

import sys,string,time,math

########################################################################
#
#  From a species tree where a node is marked as the ancestor 
#  compute all informative species pairs. The tree needs to be resolved.
#
########################################################################


##########################################
##### PARAMETERS #########################
##########################################
def usage():
	print 'Compute all informative species pairs from a species tree.'
	print ''
	print 'Usage:'
	print 'TREES_list_species_pairs.py SpeciesTreeFile SpeciesPairsFile'
	print '   Input:'
	print '    SpeciesTreeFile - the file name of the fully resolved species tree'
	print '   Output:'
	print '    SpecieaPairsFile - the filename of the species pairs file'
#enddef

if len(sys.argv) != 3:
    usage()
    sys.exit(0)

##########################################
##### READING SPECIES TREE FILE ##########
##########################################

ancestor_id = -1
duplication_id = -1

treeseq=open(sys.argv[1],"r").readline()

id=0
tree={}
pile=[]
t=0
while t<len(treeseq):
    if treeseq[t]=="(":
        id=id+1
        tree[id]=["N"+str(id),-1,-1,-1]
        if len(pile)>0:
            tree[id][1]=pile[-1]
            if(tree[pile[-1]][2]==-1):
                tree[pile[-1]][2] = id
            else:
                tree[pile[-1]][3] = id
	
        pile.append(id)
        t=t+1
    elif treeseq[t]==")":
        t=t+1
        if treeseq[t]=="@":
            ancestor_id=pile[-1]
            t=t+1
        if treeseq[t]=="$":
            duplication_id=pile[-1]
            t=t+1
        if treeseq[t]==":":
            t=t+1
            begin=t
            while treeseq[t]!="," and treeseq[t]!=")":
                t=t+1
            length=float(treeseq[begin:t])
        if treeseq[t]==";" or treeseq[t]=="\n" or treeseq == " ":
            t=t+2            
        del pile[-1]
    elif treeseq[t]==",":
        t=t+1
    else:
        id=id+1
        tree[id]=["",-1,-1,-1]
        if len(pile)>0:
            tree[id][1]=pile[-1]
            if(tree[pile[-1]][2]==-1):
                tree[pile[-1]][2] = id
            else:
                tree[pile[-1]][3] = id
        pile.append(id)
        begin=t
        while  t < len(treeseq) and treeseq[t]!="," and treeseq[t]!=")" and treeseq[t]!=":":
            t=t+1
        nom=treeseq[begin:t]
        tree[pile[-1]][0]=nom
        if t < len(treeseq) and treeseq[t]==":":
            t=t+1
            begin=t
            while treeseq[t]!="," and treeseq[t]!=")":
                t=t+1
            length=float(treeseq[begin:t])
        del pile[-1]


def partition_genomes(genome1,genome2,genome3,id,part_id):
	if(tree[id][2]==-1 and tree[id][3]==-1):
		if(part_id==1):
			genome1.append(tree[id][0])
		if(part_id==2):
			genome2.append(tree[id][0])
		if(part_id==3):
			genome3.append(tree[id][0])

	else:
		if(id==ancestor_id):
			partition_genomes (genome1,genome2,genome3,tree[id][2],2)
			partition_genomes (genome1,genome2,genome3,tree[id][3],3)
		elif(id!=duplication_id):
			partition_genomes (genome1,genome2,genome3,tree[id][2],part_id)
			partition_genomes (genome1,genome2,genome3,tree[id][3],part_id)

genomes1=[]
genomes2=[]
genomes3=[]

partition_genomes (genomes1,genomes2,genomes3,1,1)

nb_species = len(genomes1) + len(genomes2) + len(genomes3)

genome1_id=[]
genome2_id=[]
genome3_id=[]

genome_name=[]

i = 0
m={}
for s in genomes1:
	m[s]=i
	genome1_id.append(i)
	genome_name.append(s)
	i+=1
for s in genomes2:
	m[s]=i
	genome2_id.append(i)
	genome_name.append(s)
	i+=1
for s in genomes3:
	m[s]=i
	genome3_id.append(i)
	genome_name.append(s)
	i+=1

output_file = open(sys.argv[2],"w")

# outgroup
output_file.write("#outgroup\n")

for spe1 in genomes1:
    genomes2_3 =genomes2+genomes3
    for spe2  in genomes2_3:
        output_file.write(spe1+" "+spe2+"\n")

# ingroup
output_file.write("#ingroup\n")

for spe1 in genomes2:
    for spe2  in genomes3:
        output_file.write(spe1+" "+spe2+"\n")

output_file.close()
