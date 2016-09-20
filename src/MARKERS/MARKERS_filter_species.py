# ANGES 1.01, reconstructing ANcestral GEnomeS maps
# July 2012.
# Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca

import sys

#######################################################
#   MARKERS_filter_species.py
#
#   Filters an orthology blocks file by a list or species tree 
#   so that only markers with all the species listed are
#   present in the orthology blocks file.
#
#######################################################

def usage():
	print 'Filters a markers file so that only certain species type of blocks are present.'
	print ''
	print 'Usage:'
	print 'MARKERS_filter_species.py InputMarkersFile SpeciesTreeFile [0/1/2] [0/1/2] OutputMarkersFile'
	print '   Input:'
	print '    InputMarkersFile - the file name of the orthology blocks to filter'
	print '    SpeciesTreeFile - a species tree with a marked ancestor'
	print '   Output:'
	print '    OutputMarkersFile - the filename of the filtered orthology blocks'
	print '   Options:'
	print '    [0/1/2]: unique (2) or ingroup unique (1) or repeats allowed everywhere (0)'
	print '    [0/1/2]: universal (2) or ingroup universal (1) or with possible missing markers (0)'
#enddef

if len(sys.argv) != 6:
	usage()
	sys.exit()
#endif

if (sys.argv[3] != '0' and sys.argv[3] != '1' and sys.argv[3] != '2' ):
	usage()
	sys.exit()
#endif
if (sys.argv[4] != '0' and sys.argv[4] != '1'  and sys.argv[4] != '2'):
	usage()
	sys.exit()
#endif

unique=int(sys.argv[3])
universal=int(sys.argv[4])

##########################################
##### READING SPECIES TREE FILE ##########
##########################################

ancestor_id = -1
duplication_id = -1

treeseq=open(sys.argv[2],"r").readline()

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

	elif(id==ancestor_id):
		partition_genomes (genome1,genome2,genome3,tree[id][2],2)
		partition_genomes (genome1,genome2,genome3,tree[id][3],3)
	elif(id!=duplication_id):
		partition_genomes (genome1,genome2,genome3,tree[id][2],part_id)
		partition_genomes (genome1,genome2,genome3,tree[id][3],part_id)

outgroups=[]
ingroups1=[]
ingroups2=[]
partition_genomes (outgroups,ingroups1,ingroups2,1,1)
ingroups=ingroups1+ingroups2
allspecies=ingroups+outgroups
nb_ingroups=len(ingroups)
nb_outgroups=len(outgroups)

##########################################
##### FILTERING ##########################
##########################################
	
def filter_by_list(f):
	out = []
	curr = 0
	l = 0
	species=''
	ingroups_unique=1
	outgroups_unique=1
	ingroups_curr=0
	outgroups_curr=0

	for line in f:
		out.append('')
		
		if line != '\n' and line != '\r\n':
			if line[0] == '>':
				if (unique==2 and (ingroups_unique==0 or outgroups_unique==0)):
					curr = l
				#endif
				if (unique==1 and ingroups_unique==0):
					curr = l
				#endif
				if (universal==2 and ingroups_curr+outgroups_curr < nb_ingroups+nb_outgroups):
					curr = l
				#endif
				if (universal==1 and ingroups_curr < nb_ingroups):
					curr = l
				#endif
			
				l = curr
				ingroups_unique = 1
				outgroups_unique = 1
				ingroups_curr = 0
				outgroups_curr = 0
				species=''
				
				out[curr] = line
				curr += 1
			else:		
				species_prev=species
	 			sp = line.index('.')
				species = line[:sp]
				if species in allspecies:
					if species in ingroups:
						if species==species_prev:
							ingroups_unique=0
						else:
							ingroups_curr += 1						
						#endif
					else:
						if species==species_prev:
							outgroups_unique=0
						else:
							outgroups_curr += 1
						#endif
					#endif
					out[curr] = line
					curr += 1
				#endif
			#endif
		else:		
			out[curr] = line
			curr += 1
		#endif
	#endfor

	if (unique==2 and (ingroups_unique==0 or outgroups_unique==0)):
		curr = l
	#endif
	if (unique==1 and ingroups_unique==0):
		curr = l
	#endif
	if (universal==2 and ingroups_curr+outgroups_curr < nb_ingroups+nb_outgroups):
		curr = l
	#endif
	if (universal==1 and ingroups_curr < nb_ingroups):
		curr = l
	#endif
	
	f.close()
	
	return out[:curr]
#endproc

out = filter_by_list(file(sys.argv[1], 'r'))
f = file(sys.argv[5], 'w')

for l in out:
	f.write(l)
#endfor

f.close()
