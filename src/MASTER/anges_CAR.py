# ANGES 1.01, reconstructing ANcestral GEnomeS maps
# July 2012.
# Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca

import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))) + '/COMMON')

import params

import subprocess
import copy
import shutil
import threading
import time

##########################################################
#   anges.py
#
#   Main file
#
#   Usage:
#   anges_CAR.py ParametersFile
#      Input:
#       ParametersFile - the parameters file of the experiment to run
#
##########################################################

global event
global log_file

def print2(string):
	event.wait(0.05)
	
	if event.isSet():
		log_file.close()
		sys.exit()
	#endif

	try:
		print string
		log_file.write(string + '\n')
#		sys.__stdout__.write(string)
#		sys.__stdout__.flush()
	except:
		pass
	#endtry
	
	if event.isSet():	
		log_file.close()	
		sys.exit()
	#endig
#enddef

def callprocess(params):
	process = subprocess.Popen(params, bufsize=0, executable=None, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=None, close_fds=False, shell=False, cwd=None, env=None, universal_newlines=False, startupinfo=None, creationflags=0)
	
	while process.poll() == None:
		s = process.stdout.read()
	
		if s and s != "":
			print2(s)
		#endif
		
		s = process.stderr.read()
	
		if s and s != "":
			print2(s)
		#endif
	
		time.sleep(.05)
	#endif
	
	s = process.stdout.read()

	if s and s != "":
		print2(s)
	#endif
	
	s = process.stderr.read()

	if s and s != "":
		print2(s)
	#endif
#enddef


def do_c1p(params, c1p_dir, cars_dir, output_prefix, wacs, quiet, suffix, message, make_C1P, compute_PQRtree, m = False):
	acs_c1p	      = c1p_dir + "/" + output_prefix + "ACS_C1P_" + suffix	# C1P matrix file
	acs_discarded = c1p_dir + "/" + output_prefix + "ACS_DISC_" + suffix	# matrix of removed rows file
	pq_tree = cars_dir + "/" + output_prefix + "PQTREE_" + suffix # PQ-tree file	
	if params.markers_doubled:
		pq_tree_doubled = cars_dir + "/" + output_prefix + "PQTREE_DOUBLED_" + suffix # PQ-tree file	
	#endif
	if not quiet:
		print2("----> Making (weighted)ACS matrix C1P (" + message + "): " + wacs + " " + acs_discarded)
	#endif

	if m:
		callprocess(["python", code_dir + make_C1P, wacs, 'max', acs_c1p, acs_discarded])
	else:
		callprocess(["python", code_dir + make_C1P, wacs, acs_c1p, acs_discarded])
	#endif

	if not quiet:
		print2("----> Creating a PQ-tree: " + pq_tree)
	#endif

	if params.markers_doubled:
		callprocess(["python", code_dir + compute_PQRtree, acs_c1p, pq_tree_doubled, params.output_ancestor])
		if not quiet:
			print2("----> Halving PQ-tree columns") 
		#endif
		callprocess(["python", code_dir + "/C1P/C1P_halve_PQRtree.py", pq_tree_doubled, pq_tree])
	else:
		callprocess(["python", code_dir + compute_PQRtree, acs_c1p, pq_tree, params.output_ancestor])
	#endif
#enddef

# Program options

## PARAMETERS -----------------------------------------------------------------------

if len(sys.argv) > 1:
	params = params.Parameters()
	parameters_file= sys.argv[1]
	params.from_file(parameters_file)
	
	code_dir = os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))) # directory where the code is stored
	
	event = threading.Event()
	
	working_dir = os.path.dirname(parameters_file)

	# set the working directory to the directory of the parameters file
	if working_dir and working_dir != '':
		os.chdir(working_dir)
	#endif
else:
	params = app._params
#endif

## CREATING DIRECTORIES ------------------------------------------------------------------------------
output_prefix = params.output_ancestor + "_"	     # Prefix of all filenames
output_dir = params.output_dir
output_tmp = params.output_dir # directory for intermediate files

# create output directory
try:
	os.mkdir(output_dir)
except:
	pass
#endtry

## CREATING LOG FILE
output_log = output_tmp + "/" + output_prefix + "LOG"	# captures code output

log_file = open(output_log, 'w')
	
print2("----> Started")

# TRACKING
quiet = False		# True if in quiet mode

## CHECKING INCOMPATIBILITIES ------------------------------------------------------------------------

## General errors

# Missing input
if not params.markers_provided:
	print2("ERROR: no markers file")
	sys.exit(-1)
#endif
if not params.species_tree_provided:
	print2("ERROR: no species tree file")
	sys.exit(-1)
#endif

# Conflicting models
if params.c1p_circular and params.c1p_linear:
	print2("ERROR: chose either circular or linear C1P")
	sys.exit(-1)
#endif

if (not params.markers_unique) and (params.acs_ra):
	print2("ERROR: reliable adjacencies are not defined for repeated markers")
	sys.exit(-1)
#endif

if params.acs_correction >= 2 and (not params.c1p_spectral):
	print2("ERROR: corrected ACS with Xs require using the spectral seriation algorithm")
	sys.exit(-1)
#endif	
#if params.acs_correction < 2 and params.c1p_sandwich:
#	print2("ERROR: the sandwich c1p requires using corrected ACS with X's ")
#	sys.exit(-1)
#endif	


## Current limitations

if params.acs_ra and (not params.acs_sa):
	print2("ERROR: computing reliable adjacencies require computing supported adjacencies")
	# could be addressed by forcing params.acs_sa to be true if params.acs_ra is true
	sys.exit(-1)
#endif

if params.acs_weight!=1 and not params.acs_file_provided:
	print2("ERROR: weighting ACS by linear interpolation is mandatory currently ")
	sys.exit(-1)
#endif	

#if not params.markers_unique:
#	print2("ERROR: data with non unique markers are currently not handled")
#	sys.exit(-1)
#endif

#if params.c1p_sandwich:
#	print2("ERROR: C1P sandwich not handled currently")
#	sys.exit(-1)
#endif
#if params.c1p_spectral:
#	print2("ERROR: spectral seriation not handled currently")
#	sys.exit(-1)
#endif

if params.c1p_circular:
	if params.acs_aci or params.acs_sci or params.acs_mci:
		print2("ERROR: circular chromosomes can only be computed from adjacencies")
		sys.exit(-1)
        #endif
#endif

if params.c1p_spectral:
	if params.markers_doubled:
		print2("ERROR: spectral seriation can not be used with doubled markers")
		sys.exit(-1)
        #endif
	if params.c1p_circular:
		print2("ERROR: spectral seriation can not be used with circular chromosomes")
		sys.exit(-1)
        #endif
	if params.c1p_telomeres:
		print2("ERROR: spectral seriation can not be used with telomeric ACS")
		sys.exit(-1)
	#endif
#endif


## TRACKING OPTIONS -----------------------------------------------------------------------

#output = sys.stdout

#if quiet:
#	output = open(os.devnull, 'w')
##endif

## COPYING INPUT FILES -------------------------------------------------------------------------
input_dir = output_tmp + "/INPUT"

# create input directory
try:
	os.mkdir(input_dir)
except:
	pass
#endtry

if not quiet:
	print2("----> Copying input")
#endif

copied_params = copy.copy(params)

copied_parameters = input_dir + "/" + output_prefix + "PARAMETERS"
copied_markers = input_dir + "/" + output_prefix + "MARKERS"
copied_tree = input_dir + "/" + output_prefix + "SPECIES_TREE"
copied_acs = input_dir + "/" + output_prefix + "ACS"
copied_pairs = input_dir + "/" + output_prefix + "SPECIES_PAIRS"

copied_params.markers_input = output_prefix + "MARKERS"
copied_params.species_tree = output_prefix + "SPECIES_TREE"
copied_params.output_dir = ".."

if os.path.abspath(params.markers_input) != os.path.abspath(copied_markers):
	shutil.copy(params.markers_input, copied_markers)
#endif

if os.path.abspath(params.species_tree) != os.path.abspath(copied_tree):
	shutil.copy(params.species_tree, copied_tree)
#endif

if params.acs_pairs_provided:
	copied_params.acs_pairs = output_prefix + "SPECIES_PAIRS"
	if os.path.abspath(params.acs_pairs) != os.path.abspath(copied_pairs):
		shutil.copy(params.acs_pairs, copied_pairs)
	#endif
#endif
if params.acs_file_provided:
	copied_params.acs_file = output_prefix + "ACS"
	if os.path.abspath(params.acs_file) != os.path.abspath(copied_acs):
		shutil.copy(params.acs_file, copied_acs)
	#endif
#endif

copied_params.to_file(copied_parameters)

## READING MARKERS ------------------------------------------------------------------------------

markers_dir = output_tmp + "/MARKERS"
#create markers directory
try:
	os.mkdir(markers_dir)
except:
	pass
#endtry

markers_file = markers_dir + "/" + output_prefix + "MARKERS"
if not quiet:
	print2("----> Filtering markers")
#endif
callprocess(["python", code_dir + "/MARKERS/MARKERS_filter_species.py", params.markers_input, params.species_tree, str(params.markers_unique), str(params.markers_universal), markers_file])
if params.markers_doubled:
	markers_file_doubled = markers_file + "_DOUBLED"
 	if not quiet:
 		print2("----> Doubling markers")
 	#endif
 	callprocess(["python", code_dir + "/MARKERS/MARKERS_double.py", markers_file, markers_file_doubled])
	markers_file=markers_file_doubled
#endif
if not quiet:
	print2("----> Markers: " + markers_file)
#endif


## COMPUTING ACS ------------------------------------------------------------------------------
acs_dir = output_tmp + "/ACS"		# intermediate files from ACS code
mci = acs_dir + "/" + output_prefix + "MCI"		# maximal common intervals file
sci = acs_dir + "/" + output_prefix + "SCI"		# strong common intervals file
aci = acs_dir + "/" + output_prefix + "ACI"		# strong common intervals file
sa = acs_dir + "/" + output_prefix + "SA"		# supported adjacencies file
ra = acs_dir + "/" + output_prefix + "RA"		# reliable adjacencies file
acs = acs_dir + "/" + output_prefix + "ACS"		# ancestral contiguous sets file

# make ACS directory
try:
	os.mkdir(acs_dir)
except:
	pass
#endtry

# Clear ACS file
if not quiet:
	print2("----> Cleaning ACS file")
#endif
acsf = file(acs, 'w')
acsf.close()

# Computing species
if not quiet:
	print2("----> Computing species")
#endif

species=acs_dir + "/" + output_prefix + "SPECIES"
if not quiet:
	print2("------> All species: " + species)
#endif
callprocess(["python", code_dir + "/TREES/TREES_list_species.py", params.species_tree, 'ALL', species])

# collect ingroup and outgroup species
in_species = []
out_species = []
ingroup = False

for line in open(species):
	if line.rstrip() == '#ingroup':
		ingroup = True
		
		continue
	elif line.rstrip() == '#outgroup':
		ingroup = False
		
		continue
	#endif
	
	if ingroup:
		in_species.append(line.rstrip())
	else:
		out_species.append(line.rstrip())
	#endif
#endfor

# Computing species pairs
if not quiet:
	print2("----> Computing species pairs to compare")
#endif
if params.acs_pairs_provided:
	if not quiet:
		print2("------> File provided: " + params.acs_pairs)
		species_pairs=params.acs_pairs
	#endif
else:
	species_pairs=acs_dir + "/" + output_prefix + "PAIRS"
	if not quiet:
		print2("------> All informative pairs: " + species_pairs)
	#endif
	callprocess(["python", code_dir + "/TREES/TREES_list_species_pairs.py", params.species_tree, species_pairs])
#endif


ingroup = False

# Computing ACS
if not quiet:
	print2("----> Computing ancestral contiguous sets: " + acs)
#endif
for line in open(species_pairs):
	if line == None or line == '' or line == '\n':		continue
	#endif	
	sp = line.rstrip().split(' ')	
	if line.rstrip() == '#ingroup':
		ingroup = True
		
		continue
	elif line.rstrip() == '#outgroup':
		ingroup = False
		
		continue
	#endif
	if not quiet:
		print2("------> Species: " + sp[0] + " " + sp[1])
	#endif		
	
	acs_sp = acs + "_" + sp[0] + "_" + sp[1]
		
	file_list = []		# list of file to join to make the matrix
	if params.acs_mci:
		if not quiet:
			print2("--------> Computing maximum common intervals: " + mci + "_" + sp[0] + "_" + sp[1])
		#endif
		if params.markers_unique == 2:
			callprocess(["python", code_dir + "/ACS/ACS_compute_common_intervals_unique.py", markers_file, sp[0], sp[1], mci + "_" + sp[0] + "_" + sp[1], "MAX", str(int(params.c1p_circular)) ])
		else:
			callprocess(["python", code_dir + "/ACS/ACS_compute_common_intervals_nonunique.py", markers_file, sp[0], sp[1], mci + "_" + sp[0] + "_" + sp[1], "MAX", str(int(params.c1p_circular)) ])
		file_list.append(mci + "_" + sp[0] + "_" + sp[1]);
		#endif
	#endif
	if params.acs_sci:
		if not quiet:
			print2("--------> Computing strong common intervals: " + sci + "_" + sp[0] + "_" + sp[1])
		#endif	
		if params.markers_unique == 2:
			callprocess(["python", code_dir + "/ACS/ACS_compute_common_intervals_unique.py", markers_file, sp[0], sp[1], sci + "_" + sp[0] + "_" + sp[1], "STRONG", str(int(params.c1p_circular)) ])
		else:
			callprocess(["python", code_dir + "/ACS/ACS_compute_common_intervals_nonunique.py", markers_file, sp[0], sp[1], sci + "_" + sp[0] + "_" + sp[1], "STRONG", str(int(params.c1p_circular)) ])
		#endif
		file_list.append(sci + "_" + sp[0] + "_" + sp[1]);
	#endif
	if params.acs_aci:
		if not quiet:
			print2("--------> Computing all common intervals: " + aci + "_" + sp[0] + "_" + sp[1])
		#endif	
		if params.markers_unique == 2:
			callprocess(["python", code_dir + "/ACS/ACS_compute_common_intervals_unique.py", markers_file, sp[0], sp[1], aci + "_" + sp[0] + "_" + sp[1], "ALL", str(int(params.c1p_circular)) ])
		else:
			callprocess(["python", code_dir + "/ACS/ACS_compute_common_intervals_nonunique.py", markers_file, sp[0], sp[1], aci + "_" + sp[0] + "_" + sp[1], "ALL", str(int(params.c1p_circular)) ])
		#endif
 		file_list.append(aci + "_" + sp[0] + "_" + sp[1]);
 	#endif	
 	if params.acs_sa:
 		if not quiet:
 			print2("--------> Computing supported adjacencies: " + sa + "_" + sp[0] + "_" + sp[1])
 		#endif
		if params.markers_unique == 2:
			callprocess(["python", code_dir + "/ACS/ACS_compute_supported_adjacencies_unique.py", markers_file, sp[0], sp[1], sa + "_" + sp[0] + "_" + sp[1], str(int(params.c1p_circular))])		
		else:
			callprocess(["python", code_dir + "/ACS/ACS_compute_supported_adjacencies_nonunique.py", markers_file, sp[0], sp[1], sa + "_" + sp[0] + "_" + sp[1], str(int(params.c1p_circular))])		
		#endif
		file_list.append(sa + "_" + sp[0] + "_" + sp[1])
 	#endif	
 	if params.acs_ra:
 		if ingroup:
 			for sp3 in out_species:
 				if sp3 == sp[0] or sp3 == sp[1]:
 					continue
 				#endif
 			
 				if not quiet:
 					print2("--------> Computing reliable adjacencies: " + ra + "_" + sp[0] + "_" + sp[1] + "_" + sp3)
		 		#endif
 			
 				callprocess(["python", code_dir + "/ACS/ACS_compute_reliable_adjacencies_unique.py", markers_file, sa + "_" + sp[0] + "_" + sp[1], sp[0], sp[1], sp3, ra + "_" + sp[0] + "_" + sp[1] + "_" + sp3, str(int(params.c1p_circular))])
				file_list.append(ra + "_" + sp[0] + "_" + sp[1] + "_" + sp3)
				
				if not quiet:
 					print2("--------> Computing reliable adjacencies: " + ra + "_" + sp[1] + "_" + sp[0] + "_" + sp3)
		 		#endif
				
				callprocess(["python", code_dir + "/ACS/ACS_compute_reliable_adjacencies_unique.py", markers_file, sa + "_" + sp[0] + "_" + sp[1], sp[1], sp[0], sp3, ra + "_" + sp[1] + "_" + sp[0] + "_" + sp3, str(int(params.c1p_circular))])
				file_list.append(ra + "_" + sp[1] + "_" + sp[0] + "_" + sp3)
			#endfor
 		else:
 			if sp[0] in out_species:
	 			for sp3 in in_species:
	 				if sp3 == sp[0] or sp3 == sp[1]:
 						continue
 					#endif
	 			
	 				if not quiet:
	 					print2("--------> Computing reliable adjacencies: " + ra + "_" + sp[0] + "_" + sp[1] + "_" + sp3)
			 		#endif
	 				
 					callprocess(["python", code_dir + "/ACS/ACS_compute_reliable_adjacencies_unique.py", markers_file, sa + "_" + sp[0] + "_" + sp[1], sp[0], sp[1], sp3, ra + "_" + sp[0] + "_" + sp[1] + "_" + sp3, str(int(params.c1p_circular))])
					file_list.append(ra + "_" + sp[0] + "_" + sp[1] + "_" + sp3)
				#endfor
			else:
	 			for sp3 in in_species:
	 				if sp3 == sp[0] or sp3 == sp[1]:
 						continue
 					#endif
	 			
	 				if not quiet:
	 					print2("--------> Computing reliable adjacencies: " + ra + "_" + sp[1] + "_" + sp[0] + "_" + sp3)
			 		#endif
	 			
					callprocess(["python", code_dir + "/ACS/ACS_compute_reliable_adjacencies_unique.py", markers_file, sa + "_" + sp[0] + "_" + sp[1], sp[1], sp[0], sp3, ra + "_" + sp[1] + "_" + sp[0] + "_" + sp3, str(int(params.c1p_circular))])
					file_list.append(ra + "_" + sp[1] + "_" + sp[0] + "_" + sp3)
				#endfor
			#endif
 		#endif
 	#endif
 	
 	if len(file_list) == 0:
 		if not quiet:
 			print2("------> Skipping ACS computation")
 		#endif
 		
 		break
 	#endif
		
	prog = ["python", code_dir + "/ACS/ACS_join_files.py", "0"]		# join files program commandline
	prog.extend(file_list)
	prog.append(acs_sp)
	if not quiet:
		print2("------> Joining ACS files: " + acs_sp)
	#endif
 	
 	callprocess(prog)
 	 	
 	# adding telomeres
	if params.c1p_telomeres > 0:
		acs_telomeres = acs_sp + "_TEL"
		if not quiet:
			print2("------> Adding telomeres")
		#endif
	
		callprocess(["python", code_dir + "/ACS/ACS_add_telomeres.py", markers_file, sp[0], sp[1], acs_sp, acs_telomeres])
		
		acs_sp = acs_telomeres
	#endif
	
	# Correcting computed ACS
	if (not ingroup and params.markers_universal == 1) or (params.markers_universal == 0) and params.acs_correction > 0:
		if params.acs_correction == 1:
			acs_corrected = acs_sp + "_MM1"
		elif params.acs_correction == 2:
			acs_corrected = acs_sp + "_MMX"
		elif params.acs_correction == 3:
			acs_corrected = acs_sp + "_MMX1"
		#endif
		
		if not quiet:
			print2("------> Correcting computed ancestral contiguous sets: " + acs_sp)
   	    #endif		
   	     
		if params.acs_correction > 0:
			callprocess(["python", code_dir + "/ACS/ACS_add_missing_markers.py", markers_file, sp[0], sp[1], acs_sp, str(params.acs_correction), acs_corrected])
		
			acs_sp = acs_corrected
		#endif
	#endif
	
 	if not quiet:
 		print2("------> Joining ACS files: " + acs_sp + ", " + acs)
 	#endif	
	
	callprocess(["python", code_dir + "/ACS/ACS_join_files.py", acs, acs_sp, acs])
#endfor

# Reading provided ACS (if any)
if params.acs_file_provided:
 	if not quiet:
 		print2("----> Adding provided ACS files: " + params.acs_file)
 	#endif

	callprocess(["python", code_dir + "/ACS/ACS_join_files.py", acs, params.acs_file, acs])
#endif

## WEIGHTING ACS ------------------------------------------------------------------------------
weight_dir = output_tmp + "/WEIGHT"
weight_input=acs
wacs = weight_dir + "/" + output_prefix + "WACS"	# weighted ancestral contiguous sets
# make WEIGHT directory
try:
	os.mkdir(weight_dir)
except:
	pass
#endtry
if params.acs_weight==1:
	if not quiet:
		print2("----> Weighting ACS: " + wacs)
        #endif
	if params.markers_doubled:
		callprocess(["python", code_dir +"/WEIGHT/WEIGHT_weight_linear_interpolation.py", weight_input, params.species_tree, wacs, "d"])
	else:
		callprocess(["python", code_dir +"/WEIGHT/WEIGHT_weight_linear_interpolation.py", weight_input, params.species_tree, wacs])
	#endif
elif params.acs_weight==0:		# weight already given
	wacs = acs
#endif

## COMPUTING A C1P MATRIX ------------------------------------------------------------------------------
c1p_dir = output_tmp + "/C1P"		# intermediate files from C1P code
cars_dir = output_tmp + "/CARS"
pqr_tree = cars_dir + "/" + output_prefix + "PQRTREE"		# PQR-tree file
if params.markers_doubled:
	pqr_tree_doubled = cars_dir + "/" + output_prefix + "PQRTREE_DOUBLED"
#endif
# make C1P directory
try:
 	os.mkdir(c1p_dir)
except:
 	pass
#endtry
# make CARS directory
try:
 	os.mkdir(cars_dir)
except:
 	pass
#endtry
if params.c1p_linear or params.c1p_telomeres > 0:
	if not quiet:
		print2("----> Creating PQR-tree: " + pqr_tree)
	#endif

	if params.markers_doubled:
		callprocess(["python", code_dir +"/C1P/C1P_compute_PQRtree.py", wacs, pqr_tree_doubled, params.output_ancestor])
 		if not quiet:
 			print2("----> Halving PQR-tree columns") 
                #endif
		callprocess(["python", code_dir +"/C1P/C1P_halve_PQRtree.py", pqr_tree_doubled, pqr_tree])
	else:
		callprocess(["python", code_dir +"/C1P/C1P_compute_PQRtree.py", wacs, pqr_tree, params.output_ancestor])
	#endif
elif params.c1p_circular:
	if not quiet:
		print2("----> Creating PQR-tree: " + pqr_tree)
	#endif

	if params.markers_doubled:
		callprocess(["python", code_dir +"/C1P/C1P_compute_PQCRtree.py", wacs, pqr_tree_doubled, params.output_ancestor])
 		if not quiet:
 			print2("----> Halving PQR-tree columns") 
                #endif
		callprocess(["python", code_dir +"/C1P/C1P_halve_PQRtree.py", pqr_tree_doubled, pqr_tree])
	else:
		callprocess(["python", code_dir +"/C1P/C1P_compute_PQCRtree.py", wacs, pqr_tree, params.output_ancestor])
	#endif
#endif

# heuristic
if params.c1p_heuristic and params.c1p_telomeres == 0:
	if params.c1p_linear:
		do_c1p(params, c1p_dir, cars_dir, output_prefix, wacs, quiet, "HEUR", "heuristic", "/C1P/C1P_make_C1P_heuristic.py", "/C1P/C1P_compute_PQRtree.py")
	elif params.c1p_circular:
		do_c1p(params, c1p_dir, cars_dir, output_prefix, wacs, quiet, "HEUR", "heuristic", "/C1P/C1P_make_circC1P_heuristic.py", "/C1P/C1P_compute_PQCRtree.py", True)
	#endif
#endif

# branch and bound
if params.c1p_bab and params.c1p_telomeres == 0:
	if params.c1p_linear:
		do_c1p(params, c1p_dir, cars_dir, output_prefix, wacs, quiet, "BAB", "branch and bound", "/C1P/C1P_make_C1P_branch_and_bound.py", "/C1P/C1P_compute_PQRtree.py")
	elif params.c1p_circular:
		do_c1p(params, c1p_dir, cars_dir, output_prefix, wacs, quiet, "BAB", "branch and bound", "/C1P/C1P_make_circC1P_branch_and_bound.py", "/C1P/C1P_compute_PQCRtree.py", True)
	#endif
#endif

# seriation
if params.c1p_spectral:
	pq_tree = cars_dir + "/" + output_prefix + "PQTREE_SERIATION"	

	if params.acs_correction >= 2:
 		if not quiet:
 			print2("----> Computing PQ-tree from correlation matrix of a ternary matrix using parameter alpha") 
                #endif
		callprocess(["python", code_dir +"/SERIATION/SERIATION_compute_PQtree_dotproduct_correlation_matrix.py", wacs, str(params.c1p_spectral_alpha), pq_tree, params.output_ancestor])
	else:
 		if not quiet:
 			print2("----> Computing PQ-tree using spectral seriation on correlation matrix") 
                #endif
		callprocess(["python", code_dir +"/SERIATION/SERIATION_compute_PQtree_correlation_matrix.py", wacs, pq_tree, params.output_ancestor])

        #endif
#endif

# telomere heuristic
if params.c1p_telomeres == 1:
	do_c1p(params, c1p_dir, cars_dir, output_prefix, wacs, quiet, "TEL_HEUR", "heuristic", "/C1P/C1P_make_mC1P_heuristic.py", "/C1P/C1P_compute_PQRtree.py")
#endif

# telomere branch and bound, add after
if params.c1p_telomeres == 2:
	do_c1p(params, c1p_dir, cars_dir, output_prefix, wacs, quiet, "TEL_BAB1", "branch and bound, added after", "/C1P/C1P_make_mC1P_branch_and_bound_both.py", "/C1P/C1P_compute_PQRtree.py")
#endif

#telomere branch and bound, add during
if params.c1p_telomeres == 3:
	do_c1p(params, c1p_dir, cars_dir, output_prefix, wacs, quiet, "TEL_BAB2", "branch and bound, added during", "/C1P/C1P_make_mC1P_branch_and_bound.py", "/C1P/C1P_compute_PQRtree.py")
#endif

if not quiet:
 	print2("----> Done")
#endif

log_file.close()
