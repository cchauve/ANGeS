# ANGES 1.01, reconstructing ANcestral GEnomeS maps
# Program: July 2012.
# Documentation: August 2014 (v2)
# Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca


1. WHAT IS ANGES?

ANGES is a suite of Python scripts aimed at reconstructing ancestral
genome maps from conserved genomic segments in extant genomes,
including ingroups and outgroups. 

ANGES takes for minimal input (1) a species tree describing the
relationship between a group of extant genomes and with a marked
ancestral node (the ancestor) and (2) a set of homologous markers (the
markers) in the considered extant genomes. It computes conserved
genomic segments between pair of species (the Ancestral
Contiguous/Consecutive Sets, ACS), weights them according to their
pattern of conservation in the extant genomes and combines these ACS
into a set of ancestral chromosomal segments using a combinatorial
framework that was widely used for computing physical maps, the
Consecutive-Ones Property (C1P).

ANGES can reconstruct both multichromosomal linear ancestral genome
maps (for eukaryotic genomes) and unichromosomal circular ancestral
genome maps (for prokaryotic genomes). 

ANGES can handle a wide variety of sets markers, obtained from gene
families trees or multiple whole-genome alignments. It uses several
combinatorial models of ACS, together with efficient algorithms for
detecting them in pairs of genomes. It also uses several algorithms
related to the C1P, from fast heuristics to exact branch-and-bound.

The methodological principles of ANGES have been described in the
following papers.

C. Chauve, E. Tannier. A methodological framework for the
reconstruction of contiguous regions of ances- tral genomes and its
application to mammalian genomes. PLoS Computational Biology
4(11):e1000234, 2008

C. Chauve, H. Gavranovic, A. Ouangraoua, E. Tannier. Yeast ancestral
genome reconstructions: the possibilities of computational methods
II. Journal of Computational Biology 17:1097–1112, 2010.

H. Gavranovic, C. Chauve, J. Salse, E. Tannier. Mapping ancestral
genomes with massive gene loss: A matrix sandwich problem
Bioinformatics 27:i257–i265, 2011.

C. Chauve, J. Manuch, M. Patterson, R. Wittler. Tractability results
for the consecutive-ones property with multiplicity. In CPM 2011,
Lecture Notes in Comput. Sci. 6661:90–103, 2011.


2. INSTALLING ANGES.


ANGES is composed of a set of Python scripts, located in the src
directory, and nothing needs to be done for its installation.

ANGES has been developed and tested on Unix (including Linux and
MacOS) systems, using Python (http://www.python.org/) version 2.7.1+,
with the numpy library (http://numpy.scipy.org/) and the Tkinter
library for the graphical interface
(http://wiki.python.org/moin/TkInter).


3. USING ANGES.


With the graphical interface: 
python src_path/MASTER/anges_CAR_UI.py

Without the graphical interface:
python src_path/MASTER/anges_CAR.py parameter_file

In both cases, src_path is the path to access the ANGES source code.

Detailed instructions are provided in the manual anges_1.0.pdf and in
the README files in the examples directories.


4. DIRECTORY STRUCTURE.


src:      directory containing the python code

examples: a few examples using real data
