#!/bin/python
#requires the file CHOMM.py, which is simple wrapper function to run MD using OpenMM using CHARMM parameters
#=====================================================================================
#
# aux parameters (e.g. they help define the required ones, but are not themselves used by CHOMM)
name='alad' ; #prefix for output files
#platformName='CPU' ; #optional; default os 'CUDA'
#==============================
# parameters required by CHOMM (some have default values)
restart=1 ; # 0 -- start from PDB coordinates; 1 -- restart from native xml file
restartfile=None ;
restartfile=@{restartfile};

psffile='./alad_s.psf' ;
pdbfile='./alad_ms.pdb' ;
topfile='../param/alad36w.top';
paramfile='../param/alad36w.par' ;

implicitSolvent=0 ;# run OBC2 implicit solvent simulation

xmlfile=restartfile ;        # to obtain cell vectors from xml file produced with OMM (default option if restart file is provided)
#xscfile='.xsc' ;         # to obtain cell vectors from last line of xsc file
#boxfile='../alad.str'; # to obtain cell vectors from str file used in structure solvation

# specify larger box manually
#dx=112 ;
#dy=112 ;
#dz=112 ;

hmass=1;       # amu, can use heavy hydrogens
friction=1     # 1/ps, thermostat coupling
dt=2;        # timestep in fs
pmefreq=1;     # >1 requires multiple timestepping, which _dramatically_ slows down the code
cutoff=12;     # nonbonded cutoff

constraints=0;   # harmonic positional restraints for equilibration
constraintscaling=1; # to scale hatmonic restraints uniformly
consfile=pdbfile; # as in NAMD/ACEMD, this file must have identical atom ordering to that in the system topology
conscol=1 ; # 1 beta ; 2 occupancy

shake=1; # whether to constrain bonds involving hydrogens

thermostat=1;  # whether to use a thermostat
temperature=300; # kelvin
andersen=0;    # to use Andersen instead of Langevin ; (Note that I see energy up drifts quite often with andersen)
barostat=0;
pressure=1;    # units of atm
membrane_on=0; # whether to use a barostat for membrane simulations (z-axis is the membrane normal)
pme=1; # whether to use PME
pbc=1; # whether periodic boundary conditions are on

mini=0;          # whether to minimize before dynamics
ministeps=100;   # number of minimization iterations

dynamo=1
dynamoConfig='@{dynamoConfig}'
dynamoLog='@{dynamoLog}'

nsteps=@{nsteps};    # number of simulation steps
outputfreq=10000; # frequency of generating output
dcdfreq=10000;   # frequency of dcd output

outputName='scratch/'+name+'@{flag}' ;
#
from os.path import expanduser
exec(open(expanduser('~/scripts/openmm/CHOMM.py')).read())
