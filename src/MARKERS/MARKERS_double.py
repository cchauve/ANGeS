# ANGES 1.01, reconstruction ANcestral GEnomeS maps
# July 2012.
# Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca

import sys,string

#######################################################
#    MARKERS_double.py
#
#    Doubles the markers of an orthology blocks file to
#    allow the preservation of orientation.
#
#######################################################

if len(sys.argv) < 3:
	print 'Doubles the markers of a markers file.'
	print ''
	print 'Usage:'
	print 'MARKERS_double.py InputMarkersFile OutputMarkersFile'
	print '   Input:'
	print '    InputMarkersFile - the original markers file'
	print '   Output:'
	print '    OutputMarkersFile - the markers file with doubled markers'
	
	sys.exit()
#endif

name_file_blocks=sys.argv[1]
name_doubled=sys.argv[2]

file=open(name_file_blocks,"r").readlines()
output=open(name_doubled,"w")

i=0
number=1
while i<len(file):
    if file[i][0]=='>':
        id=file[i][1:file[i].find(" ")]
	comment=file[i][file[i].find(" ")+1:]	
        marker=[]
        j=i
        i=i+1
        while i<len(file) and len(file[i])>1:
            spe=file[i][:file[i].find(".")]
            chr=file[i][file[i].find(".")+1:file[i].find(":")]
            start=int(file[i][file[i].find(":")+1:file[i].find("-")])
            stop=int(file[i][file[i].find("-")+1:file[i].find(" ")])
            dir=file[i][file[i].find(" ")+1:-1]
            marker.append([spe,chr,start,stop,dir])
            i=i+1
        output.write(">"+str(2*int(id)-1)+" "+comment)
        for m in marker:
            if m[4]=="+":
                output.write(m[0]+"."+m[1]+":"+str(m[2])+"-"+str(m[2]+1)+" +\n")
            else:
                output.write(m[0]+"."+m[1]+":"+str(m[3]-1)+"-"+str(m[3])+" +\n")
        output.write("\n")        
        output.write(">"+str(2*int(id))+" "+comment)
        for m in marker:
            if m[4]=="+":
                output.write(m[0]+"."+m[1]+":"+str(m[3]-1)+"-"+str(m[3])+" +\n")
            else:
                output.write(m[0]+"."+m[1]+":"+str(m[2])+"-"+str(m[2]+1)+" +\n")
        output.write("\n")        
        i=i+1
        number=number+2

