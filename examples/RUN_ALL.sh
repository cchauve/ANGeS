# ANGES 1.01, reconstructing ANcestral GEnomeS maps
# July 2012.
# Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca

# Running all examples

cd fungal_genomes/
python ../../src/MASTER/anges_CAR.py PARAMETERS_TEL_BAB
python ../../src/MASTER/anges_CAR.py PARAMETERS_TEL_HEUR
python ../../src/MASTER/anges_CAR.py PARAMETERS_SERIATION
cd ../mammalian_genomes/
python ../../src/MASTER/anges_CAR.py PARAMETERS_OUANGRAOUA_ET_AL_2011_AMNIOTE_TEL_BAB
python ../../src/MASTER/anges_CAR.py PARAMETERS_OUANGRAOUA_ET_AL_2011_BOREOEUTHERIAN_SERIATION
python ../../src/MASTER/anges_CAR.py PARAMETERS_OUANGRAOUA_ET_AL_2011_BOREOEUTHERIAN_TEL_HEUR
python ../../src/MASTER/anges_CAR.py PARAMETERS_MA_ET_AL_2006_BOREOEUTHERIAN_HEUR
python ../../src/MASTER/anges_CAR.py PARAMETERS_MA_ET_AL_2006_BOREOEUTHERIAN_SERIATION
python ../../src/MASTER/anges_CAR.py PARAMETERS_GAVRANOVIC_ET_AL_2011_SERIATION
cd ../bacterial_genomes/
python ../../src/MASTER/anges_CAR.py PARAMETERS_BURKHOLDERIA_1
python ../../src/MASTER/anges_CAR.py PARAMETERS_BURKHOLDERIA_2
cd ../plant_genomes/
python ../../src/MASTER/anges_CAR.py PARAMETERS_MONOCOTS_BAB
python ../../src/MASTER/anges_CAR.py PARAMETERS_MONOCOTS_TEL_BAB1
python ../../src/MASTER/anges_CAR.py PARAMETERS_MONOCOTS_TEL_BAB2
cd ../