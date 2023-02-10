#!/bin/python
#requires the file CHOMM.py, which is a simple wrapper to run MD using OpenMM using CHARMM parameters
#=====================================================================================
#platformName='CPU' ; # default is CUDA; CPU is usually too slow
#
# aux parameters (e.g. they help define the required ones, but are not themselves used by CHOMM)
firstrun=@nrun   ;# initial run index
numrun=1    ;# numbef of runs
name='deca' ; #prefix for output files

#==============================
# parameters required by CHOMM (some have default values)
restart=1 ; # 0 -- start from PDB coordinates; 1 -- restart from native xml file
restartfile=None ;
restartfile='.xml';

psffile='deca22c.psf' ;
pdbfile='deca22c.pdb' ;
paramfile='par_all22_prot.prm' ;
topfile='top_all22_prot.rtf' ;

xmlfile=restartfile ;        # to obtain cell vectors from xml file produced with OMM (default option if restart file is provided)
xscfile='.xsc' ;         # to obtain cell vectors from last line of xsc file

dx=1; # specify box size manually
dy=1;
dz=1;

hmass=1;         # amu, can use heavy hydrogens
friction=1     # 1/ps, thermostat coupling
dt=1;            # timestep in fs
pmefreq=1;       # >1 requires multiple timestepping, which _dramatically_ slows down the code
cutoff=12;       # nonbonded cutoff
switchdist=10;   # optional : default is cutoff - 1.5

constraints=0;   # harmonic positional restraints for equilibration
constraintscaling=1; # to scale hatmonic restraints uniformly
consfile=pdbfile; # as in NAMD/ACEMD, this file must have identical atom ordering to that in the system topology

shake=0; # whether to constrain bonds involving hydrogens

thermostat=1;  # whether to use a thermostat
temperature=300; # kelvin
andersen=0;    # to use Andersen instead of Langevin ; (Note that I see energy up drifts quite often with andersen)
barostat=0;
pressure=1;    # units of atm
membrane_on=0; # whether to use a barostat for membrane simulations (z-axis is the membrane normal)
pme=0; # whether to use PME
pbc=0; # whether periodic boundary conditions are on

dynamo=1; # whether string plugin is active
dynamoConfig=@smconfig; # config name for string plugin
dynamoLog=@smlog; # lof name for string plugin

mini=0;          # whether to minimize before dynamics
ministeps=100;   # number of minimization iterations

nsteps=1000000;    # number of simulation steps
outputfreq=1000; # frequency of generating output
dcdfreq=10000;   # frequency of dcd output

flag=''
outputName=name+str(firstrun)+flag ;
#
# run MD simulations with different parameters one after the other
#
for i in range(numrun):
 irun=i+firstrun;
 print(" =============================");
 print(" Run ", irun, "(will quit after", firstrun+numrun-1,")");
#
 outputName=name+str(irun)+flag ;
#
 restartfile=name+str(irun-1)+flag+'.xml';
 exec(open('CHOMM.py').read())
 restart=1
# to do : check if run was successful, rerun if not
