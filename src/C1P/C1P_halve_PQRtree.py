# ANGES 1.01, reconstructing ANcestral GEnomeS maps
# July 2012.
# Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca

import sys,string

#######################################################
#	C1P_halve_PQRtree.py
#
#	Halves the leaves of a PQR-tree and adds signs (+/-).
#
#######################################################

if len(sys.argv) < 3:
	print 'Halves the leaves of a PQR-tree or PQRC-tree and adds signs.'
	print ''
	print 'Usage:'
	print 'C1P_halve_PQRtree.py InputPQRTreeFile OutputPQRTreeFile'
	print '   Input:'
	print '     InputPQRTreeFile - the PQR-tree to halve the leaves of'
	print '   Output:'
	print '     OutputPQRTreeFile - the PQR-tree with halved leaves'
	
	sys.exit()
#endif

name_doubled=sys.argv[1]
name_halved=sys.argv[2]

input_file=open(name_doubled,"r").readlines()
output_file=open(name_halved,"w")

for i in range(len(input_file)):
	if input_file[i][0]==">" or input_file[i][0]=="#":
		output_file.write(input_file[i])
	else:
		mots=input_file[i].split()
		m=0
		while m<len(mots):
			if mots[m].find("_")>=0 or mots[m] == "T":
				output_file.write(mots[m]+" ")
				m=m+1
			else:
				m1=int(mots[m])
				try:
					m2=int(mots[m+1])
				except:
					if mots[m-1].find("C")>=0:
						m2=-100
					else:
						m=m+1
						
						continue
					#endif
				#endtry
				
				if abs(m1-m2)!=1:
					if mots[m-1].find("C")>=0:
						if m1%2==0:
							output_file.write("-"+str(m1/2)+" ")
						else:
							output_file.write(str((m1+1)/2)+" ")
						#endif
						
						m=m+1
						
						continue
					else:
						print("ERROR:  C1P_halve_PQRtree " + str(m1) + " " + str(m2))
						sys.exit(0)
					#endif
				#endif
				
				if m1<m2 and m2%2==0:
					output_file.write(str(m2/2)+" ")
				elif m2<m1 and m1%2==0:
					output_file.write("-"+str(m1/2)+" ")
				else:
					if mots[m-1].find("C")>=0:
						if m1%2==0:
							output_file.write("-"+str(m1/2)+" ")
						else:
							output_file.write(str((m1+1)/2)+" ")
						#endif
						
						m=m+1
						
						continue
					else:
						print("ERROR:  C1P_halve_PQRtree " + str(m1) + " " + str(m2))
						sys.exit(0)
					#endif
				#endif
				
				m=m+2
			#endif
		#endwhile
		
		output_file.write("\n")
	#endif
#endwhile

output_file.close()
