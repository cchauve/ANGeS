# ANGES 1.01, reconstructing ANcestral GEnomeS maps
# July 2012.
# Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca

import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))) + '/COMMON')

import bm	# for precision
import datetime

from decimal import *

class Parameters:
	def __init__(self):
		pass
	#enddef
	
	def from_file(self, file_name):
		# TREES
		self.species_tree_provided = False	# True if species_tree parameter set
		self.species_tree = ""				# species tree


		# MARKERS
		self.markers_input = ""				 # file of markers
		self.markers_provided = False        # True if markers parameter set
		self.markers_doubled = False         # True if markers are doubled
		self.markers_unique = 2              # Unique markers are considered by default
		self.markers_universal = 2           # Universal markers are considered by default

		# ACS
		self.acs_sa = False		# True if using supported adjacencies
		self.acs_mci = False		# True if using maximal common intervals
		self.acs_sci = False		# True if using strong common intervals
		self.acs_aci = False		# True if using all common intervals
		self.acs_ra = False		# True if using reliable adjacencies
		self.acs_pairs_provided = False       # True if file provided
		self.acs_file_provided = False        # True if file provided
		self.acs_pairs = ""	# file of species pairs
		self.acs_file = ""	# provided acs
		self.acs_correction = 0     # 0: no correction, 1 missing markers spanned in intervals are added, 2 X markers
		self.acs_weight = 1         # Linear interpolation by default (0 means none)

		# C1P: type
		self.c1p_linear = False	# True if doing linear C1P
		self.c1p_circular = False	# True if doing circular C1P
#		self.c1p_sandwich = False    # True if sandwich C1P
		# C1P telomeres
		self.c1p_telomeres = 0    # 0: no telomere, 1: added after optimization (greedy heuristic), 2: added after optimization (bab), 3: added during optimization (bab)

		# C1P: optimization
		self.c1p_bab = False		# True if using branch and bound
		self.c1p_heuristic = False   # True if using the greedy heuristic
		self.c1p_spectral = False    # True if using the spectral seriation method
		
		# Spectral seriation options
		self.c1p_spectral_alpha = Decimal("1.0")	# Spectral seriation alpha value

		# OUTPUT
		self.output_dir = ""       # Output directory
		self.output_ancestor = "" # Name of the ancestor to infer

		# CONTROL
		self.quiet = False		# True if in quiet mode
		self.redo = False		# True to recompute everything and False to only compute what is needed/changed

		## READING PARAMETERS -----------------------------------------------------------------------

		parameters_lines=file(file_name,'r').readlines()
		for line in parameters_lines:

			if len(line)>0 and line[0]!='#':
				split = line.split()
			
				parameter = split[0]
				value = split[1].rstrip('\n')

				if parameter=="markers_file":
					self.markers_input=value
					self.markers_provided = True
				elif parameter=="tree_file":
					self.species_tree=value
					self.species_tree_provided = True
		
				elif parameter=="output_directory":
					self.output_dir=value
				elif parameter=="output_ancestor":
					self.output_ancestor=value
	
				elif parameter=="markers_doubled":
					self.markers_doubled=bool(int(value))
				elif parameter=="markers_unique":
					self.markers_unique=int(value)
				elif parameter=="markers_universal":
					self.markers_universal=int(value)
		
				elif parameter=="acs_pairs":
					self.acs_pairs=value
					self.acs_pairs_provided=True
				elif parameter=="acs_sa":
					self.acs_sa=bool(int(value))
				elif parameter=="acs_ra":
					self.acs_ra=bool(int(value))
				elif parameter=="acs_sci":
					self.acs_sci=bool(int(value))	
				elif parameter=="acs_mci":
					self.acs_mci=bool(int(value))
				elif parameter=="acs_aci":
					self.acs_aci=bool(int(value))
				elif parameter=="acs_file":
					self.acs_file=value
					self.acs_file_provided=True
				elif parameter=="acs_correction":
					self.acs_correction=int(value)
				elif parameter=="acs_weight":
					self.acs_weight=int(value)
				
				elif parameter=="c1p_linear":
					self.c1p_linear=bool(int(value))
				elif parameter=="c1p_circular":
					self.c1p_circular=bool(int(value))
#				elif parameter=="c1p_sandwich":
#					self.c1p_sandwich=bool(int(value))
				elif parameter=="c1p_telomeres":
					self.c1p_telomeres=int(value)
					
				elif parameter=="c1p_heuristic":
					self.c1p_heuristic=bool(int(value))
				elif parameter=="c1p_bab":
					self.c1p_bab=bool(int(value))	
				elif parameter=="c1p_spectral":
					self.c1p_spectral=bool(int(value))
					
				elif parameter=="c1p_spectral_alpha":
					self.c1p_spectral_alpha = Decimal(value)
					
				elif parameter=="redo":
					self.redo=bool(int(value))
					
				elif parameter=="quiet":
					self.quiet=bool(int(value))
				
				else:
					print("Unrecognized parameter: " + parameter + " " + value)
				#endif
			#endif
		#endfor
	#enddef
	
	def to_file(self, file_name):
		f = file(file_name, 'w')
		
		f.write("# Parameters File\n")
		f.write("# Project: " + self.output_ancestor + "\n") 
		f.write("# Date: " + str(datetime.date.today()) + "\n")
		f.write("#\n")
		f.write("## Input files\n")
		f.write("#\n")
		f.write("markers_file " + self.markers_input + " # file name for markers\n")
		f.write("tree_file " + self.species_tree + " # file name for species tree\n")
		f.write("#\n")
		f.write("## Output options\n")
		f.write("#\n")
		f.write("output_directory " + self.output_dir + " # directory name\n")
		f.write("output_ancestor " + self.output_ancestor + " # ancestor name\n")
		f.write("#\n")
		f.write("## Markers options\n")
		f.write("#\n")
		f.write("markers_doubled " + str(int(self.markers_doubled)) + " # 0 = original markers, 1 = doubled markers\n")
		f.write("markers_unique " + str(self.markers_unique) + " # 0 = repeated markers allowed, 1 = ingroup unique, 2 = unique\n")
		f.write("markers_universal " + str(self.markers_universal) + " # 0 = missing markers allowed, 1 = ingroup universal, 2 = universal\n")
		f.write("#\n")
		f.write("## ACS options\n")
		f.write("#\n")
		if self.acs_pairs_provided:
			f.write("acs_pairs " + self.acs_pairs + " # name of file containing the species pairs to compare\n")
		else:
			f.write("#acs_pairs  # name of file containing the species pairs to compare\n")
		#endif
		f.write("acs_sa " + str(int(self.acs_sa)) + " # supported adjacencies: 0 = not computed, 1 = computed\n")
		f.write("acs_ra " + str(int(self.acs_ra)) + " # reliable adjacencies: 0 = not computed, 1 = computed\n")
		f.write("acs_sci " + str(int(self.acs_sci)) + " # strong common intervals: 0 = not computed, 1 = computed\n")
		f.write("acs_mci " + str(int(self.acs_mci)) + " # maximal common intervals: 0 = not computed, 1 = computed\n")
		f.write("acs_aci " + str(int(self.acs_aci)) + " # all common intervals: 0 = not computed, 1 = computed\n")
		if self.acs_file_provided:
			f.write("acs_file " + self.acs_file +  " # ACS provided by user\n")
		else:
			f.write("#acs_file  # ACS provided by user\n")
		#endif
		f.write("acs_weight " + str(self.acs_weight) + " # weighting ACS: 1 = linear interpolation\n")
		f.write("acs_correction " + str(self.acs_correction) + " # Correcting for missing markers: 0 = none, 1 = adding markers spanned by intervals, 2 = X, 3 = both\n")
		f.write("#\n")
		f.write("## C1P model\n")
		f.write("#\n")
		f.write("c1p_linear " + str(int(self.c1p_linear)) + " # 1 f or working with linear chromosomes\n")
		f.write("c1p_circular " + str(int(self.c1p_circular)) + " # 1 for working with a unique circular chromosomes\n")
#		f.write("c1p_sandwich " + str(int(self.c1p_sandwich)) + " # 1 for working with ternary matrices\n")
		f.write("#\n")
		f.write("## Telomeres (model+optimization)\n")
		f.write("#\n")
		f.write("c1p_telomeres " + str(self.c1p_telomeres) + " # 0: no telomere, 1: added after optimization (greedy heuristic), 2: added after optimization (bab), 3: added during optimization (bab)\n")
		f.write("#\n")
		f.write("## C1P optimization options\n")
		f.write("#\n")
		f.write("c1p_heuristic " + str(int(self.c1p_heuristic)) + " # Using a greedy heuristic\n")
		f.write("c1p_bab " + str(int(self.c1p_bab)) + " # Using a branch-and-bound\n")
		f.write("c1p_spectral " + str(int(self.c1p_spectral)) + " # Using spectral seriation\n")
		f.write("#\n")
		f.write("## Spectral seriation options\n")
		f.write("#\n")
		f.write("c1p_spectral_alpha " + str(self.c1p_spectral_alpha) + " # Spectral seriation alpha value\n")
		f.write("#\n")
#		f.wrire("## Control Options\n")
#		f.write("redo " + str(int(self.redo)) + " # 1 to recompute everything, 0 to only compute what is needed/changed\n")
#		f.write("#\n")
		f.write("# END\n")
		
		f.close()
	#enddef
#endclass
