#!/bin/python
from __future__ import print_function
#requires the file CHOMM.py, which is simple wrapper function to run MD using OpenMM using CHARMM parameters
#=====================================================================================
#
#platformName='CPU' ; #optional; default os 'CUDA'
#==============================
# parameters required by CHOMM (some have default values)
psffile='./diala22-x.psf' ;
pdbfile='@{coordinates}' ;
paramfile='par_all22_prot.inp';
topfile='top_all22_prot.inp' ;

friction=1     # 1/ps, thermostat coupling
dt=1;        # timestep in fs
cutoff=12;     # nonbonded cutoff

thermostat=1;  # whether to use a thermostat
temperature=300; # kelvin

dynamo=1
dynamoConfig='@{input}'
dynamoLog='@{log}'

nsteps=500000;    # number of simulation steps
outputfreq=10000; # frequency of generating output
dcdfreq=10000;   # frequency of dcd output

outputName='diala.fts@{irep}';
#
from os.path import expanduser
exec(open(expanduser('CHOMM.py')).read())
