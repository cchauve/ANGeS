# ANGES 1.01, reconstructing ANcestral GEnomeS maps
# July 2012.
# Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca

import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))) + '/COMMON')

import bm	# for precision
import params

from Tkinter import *
import tkFileDialog
import tkMessageBox
import ScrolledText
import subprocess
import threading
import thread
import StringIO
import time
from decimal import *

global code_dir
global app

class Tab(LabelFrame):
	def __init__(self, master, parent, text):
		LabelFrame.__init__(self, master)
		
		self._parent = parent
		
		self["text"] = text
		
		self.createWidgets()
	#enddef
	
	def createWidgets(self):
		pass
	#endddef
		
	def show(self):
		self.pack({"fill": "both", "expand": "yes"})
	#enddef
	
	def hide(self):
		self.pack_forget()
	#enddef
#endclass

class Browsebutton(Button):
	def __init__(self, master, mode, filevar):
		Button.__init__(self, master)
	
		self._mode = mode
		self._file = filevar
		self["text"] = "Browse..."
		self["command"] = self.click
	#enddef
	
	def click(self):
		if self._mode == 0:
			ret = tkFileDialog.askopenfilename()
		elif self._mode == 1:
			ret = tkFileDialog.asksaveasfilename()
		else:
			ret = tkFileDialog.askdirectory()
		#endif
		
		if ret != "" and ret != None:
			self._file.set(ret)
		#endif
	#enddef
#enddef

class Viewbutton(Button):
	def __init__(self, master, filevar, text="View"):
		Button.__init__(self, master)
		
		self._file = filevar
		self["text"] = text
		self["command"] = self.click
	#enddef
	
	def click(self):
		f = self._file.get()
	
		if os.path.exists(f):
			if sys.platform.startswith('linux'):
				subprocess.Popen(["xdg-open", f])
			elif sys.platform == 'os2' or sys.platform == 'darwin':
				if os.path.isdir(f):
					subprocess.Popen(["open", f])
				else:
					subprocess.Popen(["edit", f])
				#endif
			#endif
		else:
			tkMessageBox.showerror(title="File not found", message="Could not find the file: \""+f+"\"")
		#endif
	#enddef
#enddef

class FileFrame(Frame):
	def __init__(self, master, text, mode, view = True, default = ""):
		Frame.__init__(self, master)
		
		self._file = StringVar()
		
		self._file.set(default)
		
		self.createWidgets(text, mode, view)
	#enddef
	
	def createWidgets(self, text, mode, view):
		if text != None:
			self._lb_text = Label(self)
			self._lb_text["text"] = text + ":"
			
			self._lb_text.pack({"side": "left"})
		#endif
		
		self._en_file = Entry(self)
		self._en_file["textvariable"] = self._file
		
		self._en_file.pack({"side": "left", "fill": "x", "expand": True})
		
		self._bb_file = Browsebutton(self, mode, self._file)
		
		self._bb_file.pack({"side": "left"})
		
		if view:
			self._vb_file = Viewbutton(self, self._file)
			
			self._vb_file.pack({"side": "left"})
		#endif
	#enddef
	
	def getFile(self):
		return self._file.get()
	#enddef
	
	def setFile(self, file_name):
		self._file.set(file_name)
	#enddef
#endclass

class EditFrame(Frame):
	def __init__(self, master, text, default = ""):
		Frame.__init__(self, master)
		
		self._text = StringVar()
		
		self._text.set(default)
		
		self.createWidgets(text)
	#enddef
	
	def createWidgets(self, text):
		if text != None:
			self._lb_text = Label(self)
			self._lb_text["text"] = text + ":"
			
			self._lb_text.pack({"side": "left"})
		#endif
		
		self._en_text = Entry(self)
		self._en_text["textvariable"] = self._text
		
		self._en_text.pack({"side": "left", "fill": "x", "expand": True})
	#enddef
	
	def getValue(self):
		return self._text.get()
	#enddef
	
	def setValue(self, val):
		self._text.set(val)
	#enddef
#endclass

class InfoTab(Tab):
	def __init__(self, master, parent):
		Tab.__init__(self, master, parent, "Info")
	#enddef
	
	def createWidgets(self):
		self._lb_app_name = Label(self)
		self._lb_app_name["text"] = "ANGES 1.01: reconstructing ANcestral GEnomeS maps"
		
		self._lb_app_name.pack({"side": "top", "anchor": "nw"})
		
		self._lb_date = Label(self)
		self._lb_date["text"] = "July 2012"
		
		self._lb_date.pack({"side": "top", "anchor": "nw"})
		
		self._lb_desc = Label(self)
		self._lb_desc["text"] = "Author: Bradley Jones, Ashok Rajaraman, Eric Tannier, Cedric Chauve"	
		self._lb_desc.pack({"side": "top", "anchor": "nw"})

		self._lb_contact = Label(self)
		self._lb_contact["text"] = "Contact: Cedric Chauve, Department of Mathematics, Simon Fraser University (Burnaby, BC, Canada), cedric.chauve@sfu.ca"	
		self._lb_contact.pack({"side": "top", "anchor": "nw"})
	#enddef
#enddef

class GenTab(Tab):
	def __init__(self, master, parent):
		Tab.__init__(self, master, parent, "I/O Files")
	#enddef
	
	def createWidgets(self):
		##########
		
		self._lf_input = LabelFrame(self)
		self._lf_input["text"] = "Input"
		
		self._lf_input.pack({"side": "top", "fill": "x"})
		
		self._ff_markers_file = FileFrame(self._lf_input, "Markers file", 0, True)
		
		self._ff_markers_file.pack({"side": "top", "fill": "x", "anchor": "nw"})
		
		self._ff_tree_file = FileFrame(self._lf_input, "Species tree file", 0, True)
		
		self._ff_tree_file.pack({"side": "top", "anchor": "nw", "fill": "x"})
		
		##########
		
		##########
		
		self._lf_output = LabelFrame(self)
		self._lf_output["text"] = "Output"
		
		self._lf_output.pack({"side": "top", "fill": "x"})
		
		self._ff_output_directory = FileFrame(self._lf_output, "Output directory", 2, False)
		
		self._ff_output_directory.pack({"side": "top", "fill": "x", "anchor": "nw"})
		
		self._ef_ancestor_name = EditFrame(self._lf_output, "Ancestor Name")
		
		self._ef_ancestor_name.pack({"side": "top", "fill": "x", "anchor": "nw"})
		
		##########
		
		##########
		
		self._var_quiet = IntVar()
		self._var_suppress = IntVar()
		self._var_needed = IntVar()
		
		self._fr_logging = Frame(self)

#		self._fr_logging.pack({"side": "top", "anchor": "w"})
		
		self._ch_quiet = Checkbutton(self._fr_logging)
		self._ch_quiet["text"] = "Quiet mode"
		self._ch_quiet["variable"] = self._var_quiet
		
		self._ch_quiet.pack({"side": "left"})
		
		self._ch_suppress = Checkbutton(self._fr_logging)
		self._ch_suppress["text"] = "Suppress subprocess output"
		self._ch_suppress["variable"] = self._var_suppress
		
		self._ch_suppress.pack({"side": "left"})
		
		##########
		
		self._ch_needed = Checkbutton(self)
		self._ch_needed["text"] = "Compute only what is needed"
		self._ch_needed["variable"] = self._var_needed
		
#		self._ch_needed.pack({"side": "top", "anchor": "w"})
	#enddef
#endclass

# unused
class TreeTab(Tab):
	def __init__(self, master, parent):
		Tab.__init__(self, master, parent, "Phylogenetic tree")
	#enddef
	
	def createWidgets(self):
		self._ff_tree = FileFrame(self, "Phylogenetic tree file", 0, True)
		
		self._ff_tree.pack({"side": "top", "anchor": "nw", "fill": "x"})
		
		##########
		
		self._var_species = IntVar()
		
		self._lf_species = LabelFrame(self)
		self._lf_species["text"] = "Compute species pairs"
		
		self._lf_species.pack({"side": "top", "anchor": "w", "fill": "x"})		
		
		self._rb_anc = Radiobutton(self._lf_species)
		self._rb_anc["text"] = "For 1 ancestor"
		self._rb_anc["variable"] = self._var_species
		self._rb_anc["value"] = 0
		
		self._rb_anc.pack({"side": "top", "anchor": "w"})
		
		self._rb_all = Radiobutton(self._lf_species)
		self._rb_all["text"] = "For all ancestors"
		self._rb_all["variable"] = self._var_species
		self._rb_all["value"] = 1
		
		self._rb_all.pack({"side": "top", "anchor": "w"})
		
		 ##########
		
		self._fr_custom = Frame(self._lf_species)
		
		self._fr_custom.pack({"side": "top", "anchor": "w", "fill": "x"})
				
		self._rb_custom = Radiobutton(self._fr_custom)
		self._rb_custom["text"] = "Provide species pairs"
		self._rb_custom["variable"] = self._var_species
		self._rb_custom["value"] = 2
		
		self._rb_custom.pack({"side": "left"})
		
		self._ff_custom = FileFrame(self._fr_custom, None, 0, True)
		
		self._ff_custom.pack({"side": "left", "fill": "x", "expand": True})
		
		 ##########
		
		##########
		
		##########
		
		self._fr_anc_name = Frame(self)
		
		self._fr_anc_name.pack({"side": "top", "anchor": "w", "fill": "x"})
	
		self._lb_anc_name = Label(self._fr_anc_name)
		self._lb_anc_name["text"] = "Ancestor name:"
		
		self._lb_anc_name.pack({"side": "left", "anchor": "w"})
		
		self._en_anc_name = Entry(self._fr_anc_name)
		
		self._en_anc_name.pack({"side": "left", "anchor": "w", "fill": "x"})
		
		self._ch_name_derive = Checkbutton(self._fr_anc_name)
		self._ch_name_derive["text"] = "Read from phylogenetic tree"
		
		self._ch_name_derive.pack({"side": "left", "anchor": "w"})
		
		##########
	#enddef
#endclass

class MarkersTab(Tab):
	def __init__(self, master, parent):
		Tab.__init__(self, master, parent, "Markers")
	#endef
	
	def createWidgets(self):
		##########
	
		self._var_markers_unique = IntVar()
		self._var_markers_unique.set(2)
		
		self._lf_markers_unique = LabelFrame(self)
		self._lf_markers_unique["text"] = "Filtering Markers: Uniqueness"
		
		self._lf_markers_unique.pack({"side": "top", "anchor": "w", "fill": "x"})		
		
		self._rb_unique = Radiobutton(self._lf_markers_unique)
		self._rb_unique["text"] = "Unique"
		self._rb_unique["variable"] = self._var_markers_unique
		self._rb_unique["value"] = 2
		self._rb_unique["command"] = self._parent.modify_universal
		
		self._rb_unique.pack({"side": "left"})
		
		self._rb_inunique = Radiobutton(self._lf_markers_unique)
		self._rb_inunique["text"] = "Ingroup unique"
		self._rb_inunique["variable"] = self._var_markers_unique
		self._rb_inunique["value"] = 1
		self._rb_inunique["command"] = self._parent.modify_universal
		
		self._rb_inunique.pack({"side": "left"})
		
		self._rb_nonunique = Radiobutton(self._lf_markers_unique)
		self._rb_nonunique["text"] = "No filtering"
		self._rb_nonunique["variable"] = self._var_markers_unique
		self._rb_nonunique["value"] = 0
		self._rb_nonunique["command"] = self._parent.modify_universal
		
		self._rb_nonunique.pack({"side": "left"})
		
		##########
	
		##########
	
		self._var_markers_universal = IntVar()
		self._var_markers_universal.set(2)
		
		self._lf_markers_universal = LabelFrame(self)
		self._lf_markers_universal["text"] = "Filtering Markers: Universality"
		
		self._lf_markers_universal.pack({"side": "top", "anchor": "w", "fill": "x"})		
		
		self._rb_universal = Radiobutton(self._lf_markers_universal)
		self._rb_universal["text"] = "Universal"
		self._rb_universal["variable"] = self._var_markers_universal
		self._rb_universal["value"] = 2
		self._rb_universal["command"] = self._parent.modify_universal
		
		self._rb_universal.pack({"side": "left"})
		
		self._rb_inuniversal = Radiobutton(self._lf_markers_universal)
		self._rb_inuniversal["text"] = "Ingroup universal"
		self._rb_inuniversal["variable"] = self._var_markers_universal
		self._rb_inuniversal["value"] = 1
		self._rb_inuniversal["command"] = self._parent.modify_universal
		
		self._rb_inuniversal.pack({"side": "left"})
		
		self._rb_nonuniversal = Radiobutton(self._lf_markers_universal)
		self._rb_nonuniversal["text"] = "No filtering"
		self._rb_nonuniversal["variable"] = self._var_markers_universal
		self._rb_nonuniversal["value"] = 0
		self._rb_nonuniversal["command"] = self._parent.modify_universal
		
		self._rb_nonuniversal.pack({"side": "left"})
		
		##########
	
		self._var_markers_doubled = IntVar()
	
		self._ch_double = Checkbutton(self)
		self._ch_double["text"] = "Double markers"
		self._ch_double["variable"] = self._var_markers_doubled
		self._ch_double["command"] = self._parent.modify_universal
		
		self._ch_double.pack({"side": "top", "anchor": "nw"})
	#enddef
#endclass

class ACSTab(Tab):
	def __init__(self, master, parent):
		Tab.__init__(self, master, parent, "Ancestral Contiguous Sets (ACS)")
	#enddef
	
	def createWidgets(self):
		##########
		
		self._var_use_acs_pairs = IntVar()
		
		self._fr_acs_pairs = Frame(self)
		
		self._fr_acs_pairs.pack({"side": "top", "anchor": "w", "fill": "x"})
				
		self._ch_acs_pairs = Checkbutton(self._fr_acs_pairs)
		self._ch_acs_pairs["text"] = "User-provided species pairs file"
		self._ch_acs_pairs["variable"] = self._var_use_acs_pairs
		
		self._ch_acs_pairs.pack({"side": "left"})
		
		self._ff_acs_pairs = FileFrame(self._fr_acs_pairs, None, 0, True)
		
		self._ff_acs_pairs.pack({"side": "left", "fill": "x", "expand": True})
		
		##########
		
		##########
		
		self._var_aci = IntVar()
		self._var_mci = IntVar()
		self._var_sci = IntVar()
		self._var_sa = IntVar()
		self._var_ra = IntVar()
		
		self._lf_set = LabelFrame(self)
		self._lf_set["text"] = "ACS models"
		
		self._lf_set.pack({"side": "top", "anchor": "w", "fill": "x"})		
				
		self._ch_aci = Checkbutton(self._lf_set)
		self._ch_aci["text"] = "All common intervals"
		self._ch_aci["variable"] = self._var_aci
		
		self._ch_aci.pack({"side": "top", "anchor": "w"})

		self._ch_mci = Checkbutton(self._lf_set)
		self._ch_mci["text"] = "Maximal common intervals"
		self._ch_mci["variable"] = self._var_mci
		
		self._ch_mci.pack({"side": "top", "anchor": "w"})
		
		self._ch_sci = Checkbutton(self._lf_set)
		self._ch_sci["text"] = "Strong common intervals"
		self._ch_sci["variable"] = self._var_sci
		
		self._ch_sci.pack({"side": "top", "anchor": "w"})
		
		self._ch_sa = Checkbutton(self._lf_set)
		self._ch_sa["text"] = "Supported adjacencies"
		self._ch_sa["variable"] = self._var_sa
		
		self._ch_sa.pack({"side": "top", "anchor": "w"})
		
		self._ch_ra = Checkbutton(self._lf_set)
		self._ch_ra["text"] = "Reliable adjacencies"
		self._ch_ra["variable"] = self._var_ra
		
		self._ch_ra.pack({"side": "top", "anchor": "w"})
		
		##########

		##########

		self._var_acs_correction = IntVar()
		
		self._lf_acs_correction = LabelFrame(self)
		self._lf_acs_correction["text"] = "Missing marker correction"
		
		self._lf_acs_correction.pack({"side": "top", "fill": "x"})
		
		self._fr_acs_correction_top = Frame(self._lf_acs_correction)
		
		self._fr_acs_correction_top.pack({"side": "top", "fill": "x"})
		
		self._lf_correct_none = Radiobutton(self._fr_acs_correction_top)
		self._lf_correct_none["text"] = "None"
		self._lf_correct_none["variable"] = self._var_acs_correction
		self._lf_correct_none["value"] = 0
		self._lf_correct_none["command"] = self._parent.modify_x
		
		self._lf_correct_none.pack({"side": "left"})
		
		self._lf_correct_span = Radiobutton(self._fr_acs_correction_top)
		self._lf_correct_span["text"] = "Add markers spanned by intervals"
		self._lf_correct_span["variable"] = self._var_acs_correction
		self._lf_correct_span["value"] = 1
		self._lf_correct_span["command"] = self._parent.modify_x
		
		self._lf_correct_span.pack({"side": "left"})
		
		self._lf_correct_x = Radiobutton(self._fr_acs_correction_top)
		self._lf_correct_x["text"] = "Add X's"
		self._lf_correct_x["variable"] = self._var_acs_correction
		self._lf_correct_x["value"] = 2
		self._lf_correct_x["command"] = self._parent.modify_x
		
		self._lf_correct_x.pack({"side": "left"})
		
		self._lf_correct_both = Radiobutton(self._lf_acs_correction)
		self._lf_correct_both["text"] = "Add markers spanned by intervals and X's"
		self._lf_correct_both["variable"] = self._var_acs_correction
		self._lf_correct_both["value"] = 3
		self._lf_correct_both["command"] = self._parent.modify_x
		
		self._lf_correct_both.pack({"side": "top", "anchor": "w"})
		
		##########

		##########
		
		self._var_use_acs_file = IntVar()
				
		self._fr_acs_file = Frame(self)
		
		self._fr_acs_file.pack({"side": "top", "anchor": "w", "fill": "x"})
		
		self._ch_acs_file = Checkbutton(self._fr_acs_file)
		self._ch_acs_file["text"] = "User-provided ACS file"
		self._ch_acs_file["variable"] = self._var_use_acs_file
		
		self._ch_acs_file.pack({"side": "left"})
		
		self._ff_acs_file = FileFrame(self._fr_acs_file, None, 0, True)
		
		self._ff_acs_file.pack({"side": "left", "fill": "x", "expand": True})
		
		##########
		
		##########
		
		self._var_acs_weight = IntVar()
		self._var_acs_weight.set(1)
		
		self._lf_acs_weight = LabelFrame(self)
		self._lf_acs_weight["text"] = "Weighting"
		
		self._lf_acs_weight.pack({"side": "top", "fill": "x"})
		
		self._rb_linear = Radiobutton(self._lf_acs_weight)
		self._rb_linear["text"] = "Linear interpolation"
		self._rb_linear["value"] = 1
		self._rb_linear["variable"] = self._var_acs_weight
		
		self._rb_linear.pack({"side": "left", "anchor": "w"})
		
		self._rb_provided = Radiobutton(self._lf_acs_weight)
		self._rb_provided["text"] = "Provided in ACS file"
		self._rb_provided["value"] = 0
		self._rb_provided["variable"] = self._var_acs_weight
		
		self._rb_provided.pack({"side": "left", "anchor": "w"})
		
		##########		
	#enddef
#endclass

class WeightTab(Tab):
	def __init__(self, master, parent):
		Tab.__init__(self, master, parent, "Weighting")
	#enddef
	
	def createWidgets(self):
		self._var_weight = IntVar()
		
		self._rb_linear = Radiobutton(self)
		self._rb_linear["text"] = "Use linear interpolation"
		self._rb_linear["value"] = 0
		self._rb_linear["variable"] = self._var_weight
		
		self._rb_linear.pack({"side": "top", "anchor": "w"})
		
		self._rb_provided = Radiobutton(self)
		self._rb_provided["text"] = "Provided in ACS"
		self._rb_provided["value"] = 1
		self._rb_provided["variable"] = self._var_weight
		
		self._rb_provided.pack({"side": "top", "anchor": "w"})
	#enddef
#endclass

class C1PTab(Tab):
	def __init__(self, master, parent):
		Tab.__init__(self, master, parent, "C1P")
	#enddef
	
	def modify(self):
		if self._var_c1p_type.get() == 3:
			self._rb_after_heur["state"] = NORMAL
			self._rb_after_bab["state"] = NORMAL
			self._rb_during_bab["state"] = NORMAL
			self._parent._acs._ch_aci["state"] = NORMAL
			self._parent._acs._ch_sci["state"] = NORMAL
			self._parent._acs._ch_mci["state"] = NORMAL
			self._parent._acs._var_aci.set(0)
			self._parent._acs._var_mci.set(0)
			self._parent._acs._var_sci.set(0)
			self._ch_bab["state"] = DISABLED
			self._ch_heur["state"] = DISABLED
			self._ch_ser["state"] = DISABLED
		else:
			self._rb_after_heur["state"] = DISABLED
			self._rb_after_bab["state"] = DISABLED
			self._rb_during_bab["state"] = DISABLED
			
			if self._var_c1p_type.get() == 2:
				self._ch_bab["state"] = DISABLED
				self._ch_heur["state"] = DISABLED
				self._var_bab.set(0)
				self._var_heur.set(0)
				
			else:
				self._ch_bab["state"] = NORMAL
				self._ch_heur["state"] = NORMAL
				if self._var_c1p_type.get() == 1:
					self._parent._acs._ch_aci["state"] = DISABLED
					self._parent._acs._ch_sci["state"] = DISABLED
					self._parent._acs._ch_mci["state"] = DISABLED
					self._parent._acs._var_aci.set(0)
					self._parent._acs._var_mci.set(0)
					self._parent._acs._var_sci.set(0)
				else:
					self._parent._acs._ch_aci["state"] = NORMAL
					self._parent._acs._ch_sci["state"] = NORMAL
					self._parent._acs._ch_mci["state"] = NORMAL
					self._parent._acs._var_aci.set(0)
					self._parent._acs._var_mci.set(0)
					self._parent._acs._var_sci.set(0)
			#endif
			if not(self._parent._markers._var_markers_doubled.get()) and self._var_c1p_type.get() != 1:
				self._ch_ser["state"] = NORMAL
			else:
				self._ch_ser["state"] = DISABLED
				self._var_ser.set(0)
		
		#endif
		
		if self._var_ser.get():
			if self._parent._acs._var_acs_correction.get()>=2:
			    self._ef_ser_alpha._en_text["state"] = NORMAL
			else:
			    self._ef_ser_alpha._en_text["state"] = DISABLED
		else:
			self._ef_ser_alpha._en_text["state"] = DISABLED
		#endif
		
		if self._parent._markers._var_markers_unique.get() == 0:
			self._rb_telomere["state"] = DISABLED
			if self._var_c1p_type.get() == 3:
				self._var_c1p_type.set(0)
			#endif
		#endif
	#enddef
	
	def createWidgets(self):
		##########
		
		self._var_c1p_type = IntVar()
		
		self._lf_c1p_type = LabelFrame(self)
		self._lf_c1p_type["text"] = "Ancestral Map Model"
		
		self._lf_c1p_type.pack({"side": "top", "anchor": "w", "fill": "x"})		
		
		self._rb_c1p = Radiobutton(self._lf_c1p_type)
		self._rb_c1p["text"] = "Linear Chromosome(s)"
		self._rb_c1p["variable"] = self._var_c1p_type
		self._rb_c1p["value"] = 0
		self._rb_c1p["command"] = self.modify
		
		self._rb_c1p.pack({"side": "left"})
		
		self._rb_circular = Radiobutton(self._lf_c1p_type)
		self._rb_circular["text"] = "Circular Chromosome"
		self._rb_circular["variable"] = self._var_c1p_type
		self._rb_circular["value"] = 1
		self._rb_circular["command"] = self.modify
		
		self._rb_circular.pack({"side": "left"})
		
#		self._rb_sandwich = Radiobutton(self._lf_c1p_type)
#		self._rb_sandwich["text"] = "Sandwich"
#		self._rb_sandwich["variable"] = self._var_c1p_type
#		self._rb_sandwich["value"] = 2
#		self._rb_sandwich["command"] = self.modify
#		
#		self._rb_sandwich.pack({"side": "left"})
		
		self._rb_telomere = Radiobutton(self._lf_c1p_type)
		self._rb_telomere["text"] = "Linear Chromosome(s) with Telomeres"
		self._rb_telomere["variable"] = self._var_c1p_type
		self._rb_telomere["value"] = 3
		self._rb_telomere["command"] = self.modify
		
		self._rb_telomere.pack({"side": "left"})
		
		##########
		
		##########
		
		self._var_bab = IntVar()
		self._var_heur = IntVar()
		self._var_ser = IntVar()
		
		self._lf_method = LabelFrame(self)
		self._lf_method["text"] = "Ancestral Map Computation Method"
		
		self._lf_method.pack({"side": "top", "fill": "x"})
		
		self._ch_bab = Checkbutton(self._lf_method)
		self._ch_bab["text"] = "Branch and bound"
		self._ch_bab["variable"] = self._var_bab
		self._ch_bab["command"] = self.modify
		
		self._ch_bab.pack({"side": "left"})
		
		self._ch_heur = Checkbutton(self._lf_method)
		self._ch_heur["text"] = "Greedy Heuristic"
		self._ch_heur["variable"] = self._var_heur
		self._ch_heur["command"] = self.modify
		
		self._ch_heur.pack({"side": "left"})
		
		self._ch_ser = Checkbutton(self._lf_method)
		self._ch_ser["text"] = "Spectral Seriation"
		self._ch_ser["variable"] = self._var_ser
		self._ch_ser["command"] = self.modify
		
		self._ch_ser.pack({"side": "left"})
		
		##########

		self._lf_ser_alpha = LabelFrame(self)
		self._lf_ser_alpha["text"] = "Spectral seriation parameters"
		self._lf_ser_alpha.pack({"side": "top", "fill": "x"})
		self._ef_ser_alpha = EditFrame(self._lf_ser_alpha, u"\u03b1 value for X-markers", "1.0")
		
		self._ef_ser_alpha.pack({"side": "top", "anchor": "w"})
		
		##########
		
		self._var_telomere_method = IntVar()
		
		self._lf_telomere_method = LabelFrame(self)
		self._lf_telomere_method["text"] = "Linear Chromosome(s) with Telomeres: Computation Method"
		
		self._lf_telomere_method.pack({"side": "top", "fill": "x"})
		
		self._rb_after_heur = Radiobutton(self._lf_telomere_method)
		self._rb_after_heur["text"] = "Added after C1P optimization (Heuristic)"
		self._rb_after_heur["variable"] = self._var_telomere_method
		self._rb_after_heur["value"] = 1
		self._rb_after_heur["command"] = self.modify
		
		self._rb_after_heur.pack({"side": "top", "anchor": "w"})
		
		self._rb_after_bab = Radiobutton(self._lf_telomere_method)
		self._rb_after_bab["text"] = "Added after C1P optimization (Branch and bound)"
		self._rb_after_bab["variable"] = self._var_telomere_method
		self._rb_after_bab["value"] = 2
		self._rb_after_bab["command"] = self.modify
		
		self._rb_after_bab.pack({"side": "top", "anchor": "w"})
		
		self._rb_during_bab = Radiobutton(self._lf_telomere_method)
		self._rb_during_bab["text"] = "Added during C1P optimization (Branch and bound)"
		self._rb_during_bab["variable"] = self._var_telomere_method
		self._rb_during_bab["value"] = 3
		self._rb_during_bab["command"] = self.modify
		
		self._rb_during_bab.pack({"side": "top", "anchor": "w"})
		
		##########
		
		##########
		
		self._fr_PQtree = Frame(self)
		
#		self._fr_PQtree.pack({"side": "top", "anchor": "w"})
		
		self._ch_PQRtree = Checkbutton(self._fr_PQtree)
		self._ch_PQRtree["text"] = "Compute PQR-tree"
		
		self._ch_PQRtree.pack({"side": "left"})
		
		self._ch_PQtree = Checkbutton(self._fr_PQtree)
		self._ch_PQtree["text"] = "Compute PQ-tree"
		
		self._ch_PQtree.pack({"side": "left"})
		
		##########
		
		##########
		
		self._lf_tucker = LabelFrame(self)
		self._lf_tucker["text"] = "Tucker patterns"
				
#		self._lf_tucker.pack({"side": "top", "fill": "x"})
		
		self._ch_PQRtree = Checkbutton(self._lf_tucker)
		self._ch_PQRtree["text"] = "Compute Tucker patterns"
		
		self._ch_PQRtree.pack({"side": "left", "anchor": "w"})
		
		self._ch_PQtree = Checkbutton(self._lf_tucker)
		self._ch_PQtree["text"] = "Count Tucker patterns"
		
		self._ch_PQtree.pack({"side": "left", "anchor": "w"})
		
		##########
	#enddef
#endclass

#######################
# Logger class
# Displays output from experiments and writes the output to a log file.
#######################
class Logger(Toplevel):
	def __init__(self, event, thread_id, folder, log_file_name, title):
		Toplevel.__init__(self)
				
		self.title(title)
		
		self._event = event
		self._event.clear()
		
		self._thread = thread_id
		
		# create output directory
		try:
			os.mkdir(folder)
		except:
			pass
		#endtry
		
#		self._log_file_name = log_file_name
#		self._log_file = None
#		try:
#			self._log_file = open(log_file_name, "w")
#		except:
#			self._log_file = None
#		#endtry
		
		self._var_folder = StringVar()
		self._var_folder.set(folder)
		
		self._log_txt_pos = 0
					
		self.createWidgets()
	#enddef
	
	def createWidgets(self):		
		self._st_text = ScrolledText.ScrolledText(self)
		
		self._st_text.pack({"side": "top", "fill": "both"})
		
		self._fr_btns = Frame(self)
		
		self._fr_btns.pack({"side": "top", "fill": "both"})
		
		self._btn_cancel = Button(self._fr_btns)
		self._btn_cancel["text"] = "Stop"
		self._btn_cancel["command"] = self.stop
		
		self._btn_cancel.pack({"side": "right", "anchor": "e"})
		
		self._btn_view = Viewbutton(self._fr_btns, self._var_folder, "View project folder")
		
		self._btn_view.pack({"side": "right", "anchor": "e"})
	#enddef
	
	def destroy(self):
		self.stop()
		
#		if self._log_file:
#			self._log_file.close()
#		#endif
		
		Toplevel.destroy(self)
	#enddef

	def stop(self):	
#		if self._log_file:
#			self._log_file.close()
#		#endif
		
		try:
			if self._event.is_set():
				self._st_text.insert(END, "Stopped\n")
				self._st_text.see(END)
			#endif
		except:
			pass
		#endif
		
		self._event.set()
	#enddef
	
	def update_logger(self):
		sys.stdout.flush()
				
		whole_txt = sys.stdout.getvalue()
		
		txt = whole_txt[self._log_txt_pos:]
		
		self._log_txt_pos = len(whole_txt)
	
		if txt and txt != "":
			self._st_text.insert(END, txt)
			self._st_text.see(END)
		
#			if self._log_file:
#				self._log_file.write(txt)
#				self._log_file.flush()
#			#endif
		#endif
	
#		if not self._log_file:
#			try:
#				self._log_file = open(self._log_file_name, "w")
#			except:
#				self._log_file = None
#			#endtry
#		#endif
#		
#		if self._log_file:
#			try:
#				sys.stdout.flush()
#			
#				txt = self._log_file.read()
#								
#				if txt and txt != "":
#					if '&' in txt:
#						lines = txt.split('&')
#						pos = END
#											
#						if txt[0] == '&':
#							try:
#								while self._st_text.get(END, None) != '\n':
#									sys,__stdout__.write(self._st_text.get(END, None) + '\n')
#									self._st_text.delete(END)
#								#endwhile							
#							except:
#								pass
#							#endtry
#						#endif
#						
#						sys.__stdout__.write('-->deleted\n')
#						
#						for line in lines:
#							if '\n' in line:							
#								sys.__stdout__.write('<<' + line + '>>\n')
#							
#								self._st_text.insert(END, line)
#							#endif
#						#endif
#						
#						if '\n' not in lines[-1]:
#							sys.__stdout__.write('<<<' + lines[-1] + '>>>\n')
#						
#							self._st_text.insert(END, lines[-1])
#						#endif
#					else:
#						self._st_text.insert(END, txt)
#					#endif
#					
#					self._st_text.insert(END, txt)
#					self._st_text.see(END)
#				#endif
#			except:			
#				pass
#			#endtry
#		#endif
		
		if not self._thread.isAlive():	
#			try:
#				self._st_text.insert(END, "Exited on error\n")
#				self._st_text.see(END)
#			except:
#				pass
#			#endif
		
			self._event.set()
		#endif
	#enddef
#endclass

class Application(Frame):
	def __init__(self, master=None):
		Frame.__init__(self, master)
		
		self._master = master
		self._params = params.Parameters()
		self._parameters_file = ""
		
		self.createWidgets()
		self.pack({"fill": "both", "expand": "yes"})
	#enddef
	
	def switchTab(self, tab, btn):
		self._active_btn["state"] = NORMAL
		self._active_btn["relief"] = RAISED
		btn["state"] = DISABLED
		btn["relief"] = SUNKEN
	
		self._active.hide()
		tab.show()
		
		self._active = tab
		self._active_btn = btn
	#enddef
	
	def clickInfo(self):
		self.switchTab(self._info, self._info_btn)
	#enddef
	
	def clickGen(self):
		self.switchTab(self._gen, self._gen_btn)
	#enddef
	
	def clickTree(self):
		self.switchTab(self._tree, self._tree_btn)
	#enddef
	
	def clickMarkers(self):
		self.switchTab(self._markers, self._markers_btn)
	#enddef
	
	def clickACS(self):
		self.switchTab(self._acs, self._acs_btn)
	#enddef
		
	def clickWeight(self):
		self.switchTab(self._weight, self._weight_btn)
	#enddef
	
	def clickC1P(self):
		self.switchTab(self._c1p, self._c1p_btn)
	#enddef
	
	def modify_x(self):
		if self._acs._var_acs_correction.get() >= 2:
			self._c1p._rb_c1p["state"] = DISABLED
			self._c1p._rb_circular["state"] = DISABLED
#			self._c1p._rb_sandwich["state"] = NORMAL
			self._c1p._rb_telomere["state"] = DISABLED
		
			self._c1p._var_c1p_type.set(2)
		else:
			self._c1p._rb_c1p["state"] = NORMAL
			self._c1p._rb_circular["state"] = NORMAL
#			self._c1p._rb_sandwich["state"] = DISABLED
			self._c1p._rb_telomere["state"] = NORMAL
		#endif
		
		self._c1p.modify()
	#enddef
	
	def modify_universal(self):
		if self._markers._var_markers_universal.get() == 2:
			self._acs._lf_correct_none["state"] = DISABLED
			self._acs._lf_correct_span["state"] = DISABLED
			self._acs._lf_correct_x["state"] = DISABLED
			self._acs._lf_correct_both["state"] = DISABLED
		
			self._acs._var_acs_correction.set(0)
			
			if self._markers._var_markers_unique.get() !=2:
				self._acs._ch_ra["state"] = DISABLED
			else:
				self._acs._ch_ra["state"] = NORMAL
					
		else:
			self._acs._ch_ra["state"] = DISABLED
			if self._markers._var_markers_unique.get() == 2:
				self._acs._lf_correct_span["state"] = NORMAL
				self._acs._lf_correct_both["state"]= NORMAL
			else:
				self._acs._lf_correct_span["state"] = DISABLED
				self._acs._lf_correct_both["state"]= DISABLED
			#endif
			
			self._acs._lf_correct_none["state"] = NORMAL
			self._acs._lf_correct_x["state"] = NORMAL
		#endif
		
		self.modify_x()
	#enddef


		
		
	def commit(self):
		self._gen._ff_tree_file.setFile(self._params.species_tree)
		self._gen._ff_markers_file.setFile(self._params.markers_input)
		self._gen._ff_output_directory.setFile(self._params.output_dir)
		self._gen._ef_ancestor_name.setValue(self._params.output_ancestor)
#		self._gen._var_quiet.setValue(int(self._params.quiet))
#		self._gen._var_suppress.getValue(int(self._params.supress))
		
		self._markers._var_markers_unique.set(self._params.markers_unique)
		self._markers._var_markers_universal.set(self._params.markers_universal)
		self._markers._var_markers_doubled.set(int(self._params.markers_doubled))
		
		self._acs._var_use_acs_pairs.set(int(self._params.acs_pairs_provided))
		self._acs._ff_acs_pairs.setFile(self._params.acs_pairs)
		self._acs._var_aci.set(int(self._params.acs_aci))		
		self._acs._var_mci.set(int(self._params.acs_mci))
		self._acs._var_sci.set(int(self._params.acs_sci))
		self._acs._var_sa.set(int(self._params.acs_sa))
		self._acs._var_ra.set(int(self._params.acs_ra))
		self._acs._var_acs_correction.set(self._params.acs_correction)
		self._acs._var_use_acs_file.set(int(self._params.acs_file_provided))
		self._acs._ff_acs_file.setFile(self._params.acs_file)
		self._acs._var_acs_weight.set(self._params.acs_weight)
		if self._params.c1p_linear:
			self._c1p._var_c1p_type.set(0)
		elif self._params.c1p_circular:
			self._c1p._var_c1p_type.set(1)
#		elif self._params.c1p_sandwich:
#			self._c1p._var_c1p_type.set(2)
		elif self._params.c1p_telomeres > 0:
			self._c1p._var_c1p_type.set(3)
			self._c1p.modify()
			self._c1p._var_telomere_method.set(self._params.c1p_telomeres)
		#endif
		self._c1p._var_bab.set(int(self._params.c1p_bab))
		self._c1p._var_heur.set(int(self._params.c1p_heuristic))
		self._c1p._var_ser.set(int(self._params.c1p_spectral))
		self._c1p._ef_ser_alpha.setValue(str(self._params.c1p_spectral_alpha))
		self.modify_universal()
		
	#enddef
	
	def update_params(self):
		self._params.species_tree_provided = True
		self._params.species_tree = self._gen._ff_tree_file.getFile()
		self._params.markers_provided = True
		self._params.markers_input = self._gen._ff_markers_file.getFile()
		self._params.output_dir = self._gen._ff_output_directory.getFile()
		self._params.output_ancestor = self._gen._ef_ancestor_name.getValue()
#		self._params.quiet = self._gen._var_quiet.setValue() == 1
#		self._params.supress = self._gen._var_suppress.getValue() == 1
		
		self._params.markers_unique = self._markers._var_markers_unique.get()
		self._params.markers_universal = self._markers._var_markers_universal.get()
		self._params.markers_doubled = self._markers._var_markers_doubled.get() == 1
		
		self._params.acs_pairs_provided = self._acs._var_use_acs_pairs.get() == 1
		self._params.acs_pairs = self._acs._ff_acs_pairs.getFile()
		self._params.acs_aci = self._acs._var_aci.get() == 1
		self._params.acs_mci = self._acs._var_mci.get() == 1
		self._params.acs_sci = self._acs._var_sci.get() == 1
		self._params.acs_sa = self._acs._var_sa.get() == 1
		self._params.acs_ra = self._acs._var_ra.get() == 1
		self._params.acs_correction = self._acs._var_acs_correction.get()
		self._params.acs_file_provided = self._acs._var_use_acs_file.get() == 1
		self._params.acs_file = self._acs._ff_acs_file.getFile()
		self._params.acs_weight = self._acs._var_acs_weight.get()
		
		self._params.c1p_linear = self._c1p._var_c1p_type.get() == 0
		self._params.c1p_circular = self._c1p._var_c1p_type.get() == 1
#		self._params.c1p_sandwich = self._c1p._var_c1p_type.get() == 2
		if self._c1p._var_c1p_type.get() == 3:
			self._params.c1p_telomeres = self._c1p._var_telomere_method.get()
		else:
			self._params.c1p_telomeres = 0
		#endif
		self._params.c1p_bab = self._c1p._var_bab.get() == 1
		self._params.c1p_heuristic = self._c1p._var_heur.get() == 1
		self._params.c1p_spectral = self._c1p._var_ser.get() == 1
		self._params.c1p_spectral_alpha = Decimal(self._c1p._ef_ser_alpha.getValue())
		#self.modify_universal()
	#enddef
	
	def reset(self):
		self._gen._ff_tree_file.setFile('')
		self._gen._ff_markers_file.setFile('')
		self._gen._ff_output_directory.setFile('')
		self._gen._ef_ancestor_name.setValue('')
#		self._gen._var_quiet.setValue(int(self._params.quiet))
#		self._gen._var_suppress.getValue(int(self._params.supress))
		
		self._markers._var_markers_unique.set(3)
		self._markers._var_markers_universal.set(3)
		self._markers._var_markers_doubled.set(0)
		
		self._acs._var_use_acs_pairs.set(0)
		self._acs._ff_acs_pairs.setFile('')
		self._acs._var_aci.set(0)		
		self._acs._var_mci.set(0)
		self._acs._var_sci.set(0)
		self._acs._var_sa.set(0)
		self._acs._var_ra.set(0)
		self._acs._var_acs_correction.set(0)
		self._acs._var_use_acs_file.set(0)
		self._acs._ff_acs_file.setFile('')
		self._acs._var_acs_weight.set(1)
		self._c1p._var_c1p_type.set(0)

		self._c1p._var_bab.set(0)
		self._c1p._var_heur.set(0)
		self._c1p._var_ser.set(0)
		self._c1p._ef_ser_alpha.setValue(1.0)
		self.modify_universal()
	#enddef
	
	def load(self):
		file_name = tkFileDialog.askopenfilename()
		
		if file_name:
			self._params_file = file_name
			
			self._params.from_file(file_name)
			
			wrk_dir = os.path.dirname(self._params_file)
			
			os.chdir(wrk_dir)
			
			self.commit()
		#endif
	#enddef
	
	def save(self):
		file_name = tkFileDialog.asksaveasfilename()
		
		if file_name:
			self._params_file = file_name
			self.update_params()
			
			wrk_dir = os.path.dirname(self._params_file)
			
			os.chdir(wrk_dir)
			
			self._params.to_file(file_name)
		#endif
	#enddef
	
	def run(self):
#		if self.is_updated():
#			subprocess.call(["python", code_dir + "/MASTER/anges_CAR.py", self._params_file])
#		else:
#			tkMessageBox.showwarning(title="Warning", message="Parameters file not saved")
#		#endif
		self.update_params()

		try:
			os.mkdir(self._params.output_dir)
		except:
			pass
		#endtry
				
		log_file_name = self._params.output_dir + "/" + self._params.output_ancestor + "_LOG"
		
#		log_file = open(log_file_name, 'w')
		log_file = StringIO.StringIO()
		sys.stdout = log_file
		sys.stderr = log_file
				
		logger = Logger(threading.Event(), None, self._params.output_dir, log_file_name, "Log - " + self._params.output_ancestor)
#		thread_id = thread.start_new_thread(execute, ())
		thread_id = threading.Thread(target=execute, kwargs={"logger": logger})
		
		logger._thread = thread_id
		thread_id.daemon = True
		
		thread_id.start()
						
		update(logger)
	#enddef
	
	def createWidgets(self):
		self._btns = Frame(self)
		self._tab = Frame(self)
		self._bot_btns = Frame(self)
		
		self._info = InfoTab(self._tab, self)
		
		self._info_btn = Button(self._btns)
		self._info_btn["text"] = "Info"
		self._info_btn["command"] = self.clickInfo
		
		self._gen = GenTab(self._tab, self)
		
		self._gen_btn = Button(self._btns)
		self._gen_btn["text"] = "I/O Files"
		self._gen_btn["command"] = self.clickGen
		
#		self._tree = TreeTab(self._tab, self)
		
#		self._tree_btn = Button(self._btns)
#		self._tree_btn["text"] = "Phylogenetic tree"
#		self._tree_btn["command"] = self.clickTree
		
		self._markers = MarkersTab(self._tab, self)
		
		self._markers_btn = Button(self._btns)
		self._markers_btn["text"] = "Markers"
		self._markers_btn["command"] = self.clickMarkers
		
		self._acs = ACSTab(self._tab, self)
	
		self._acs_btn = Button(self._btns)
		self._acs_btn["text"] = "Ancestral Contiguous Sets"
		self._acs_btn["command"] = self.clickACS

#		self._weight = WeightTab(self._tab, self)
						
#		self._weight_btn = Button(self._btns, self)
#		self._weight_btn["text"] = "Weighting"
#		self._weight_btn["command"] = self.clickWeight
		
		self._c1p = C1PTab(self._tab, self)
				
		self._c1p_btn = Button(self._btns)
		self._c1p_btn["text"] = "Ancestral Maps"
		self._c1p_btn["command"] = self.clickC1P
		
#		self._as = ASTab(self._tab, self)
		
#		self._as_btn = Button(self._btns)
#		self._as_btn["text"] = "AS"
#		self._as_btn["command"] = self.clickAS
		
#		self._asg = AGSTab(self._tab, self)
		
#		self._asg_btn = Button(self._btns)
#		self._asg_btn["text"] = "ASG"
#		self._asg_btn["command"] = self.clickASG
		
		self._btn_reset_params = Button(self._bot_btns)
		self._btn_reset_params["text"] = "Reset Parameters"
		self._btn_reset_params["command"] = self.reset

		self._btn_load_params = Button(self._bot_btns)
		self._btn_load_params["text"] = "Load parameters file"
		self._btn_load_params["command"] = self.load

#		self._btn_save_params = Button(self._bot_btns)
#		self._btn_save_params["text"] = "Save parameters file"
#		self._btn_save_params["command"] = self.save
		
		self._btn_run_script = Button(self._bot_btns)
		self._btn_run_script["text"] = "Run script"
		self._btn_run_script["command"] = self.run
		
		##########
		
		self._btns.pack({"side": "top"})
		
		self._info_btn.pack({"side": "left"})
		self._gen_btn.pack({"side": "left"})
#		self._tree_btn.pack({"side": "left"})
		self._markers_btn.pack({"side": "left"})
		self._acs_btn.pack({"side": "left"})
#		self._weight_btn.pack({"side": "left"})
		self._c1p_btn.pack({"side": "left"})
#		self._as_btn.pack({"side": "left"})
#		self._asg_btn.pack({"side": "left"})
		
		self._tab.pack({"side": "top", "fill": "both", "expand": "yes"})
		
		self._bot_btns.pack({"side": "top", "fill": "both"})
		
		self._btn_run_script.pack({"side": "right"})
#		self._btn_save_params.pack({"side": "right"})
		self._btn_load_params.pack({"side": "right"})
		self._btn_reset_params.pack({"side": "right"})
		
		self._active = self._gen
		self._active_btn = self._gen_btn
		
		self.switchTab(self._info, self._info_btn)
		
		self.modify_universal()
	#enddef
#endclass

def execute(logger):
	execfile(code_dir + "/MASTER/anges_CAR.py", {"event": logger._event, "app": app, "code_dir": code_dir})
		
	time.sleep(0.5)
	
	logger._event.set()
#enddef

def update(logger):
	while not logger._event.isSet():
		logger.update_logger()
		try:
			logger.update()
		except:
			sys.__stdout__.write("break\n")
			sys.__stdout__.flush()
		
			break
		#endtry
		
		logger._event.wait(0.05)
	#endwhile
	
	time.sleep(0.1)
	
	try:
		logger.update_logger()
	except:
		pass	
	#endtry
	
	logger.stop()
	
	sys.stdout.close()
				
	sys.stdout = sys.__stdout__
	sys.stderr = sys.__stderr__
#enddef

#sys.setrecursionlimit(1000)

code_dir = os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))) # code directory

root = Tk()
app = Application(master=root)

root.title("ANGES 1.0")

app.mainloop()

