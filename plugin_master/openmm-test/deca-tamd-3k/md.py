#!/bin/python
from __future__ import print_function
import fileinput
import sys
from os import mkdir, path
#requires the file CHOMM.py, which is simple wrapper function to run MD using OpenMM using CHARMM parameters
#=====================================================================================
#
# aux parameters (e.g. they help define the required ones, but are not themselves used by CHOMM)
firstrun=16   ;# initial run index
numrun=10
name='deca'
#platformName='CPU' ; #optional; default is 'CUDA'
#==============================
# parameters required by CHOMM (some have default values)
if not (path.exists('scratch')):
 mkdir('scratch');
psffile='./'+name+'22c.psf' ;
pdbfile='./'+name+'22f.pdb' ;
paramfile='par_all22_prot.prm' ;
topfile='top_all22_prot.rtf' ;
3
implicitSolvent=0 ;# run OBC2 implicit solvent simulation

xscfile='.xsc' ;         # to obtain cell vectors from last line of xsc file
boxfile='./struc/'+name+'.str'; # to obtain cell vectors from str file used in structure solvation

# specify larger box manually
#dx=72 ;
#dy=72 ;
#dz=72 ;

hmass=1;       # amu, can use heavy hydrogens
friction=1   # 1/ps, thermostat coupling
dt=1;          # timestep in fs
pmefreq=1;     # >1 requires multiple timestepping, which _dramatically_ slows down the code
cutoff=12;      # nonbonded cutoff
switchdist=10 ; # (optional) switching distance

constraints=0;   # harmonic positional restraints for equilibration
constraintscaling=1; # to scale hatmonic restraints uniformly
consfile=pdbfile; # as in NAMD/ACEMD, this file must have identical atom ordering to that in the system topology
conscol=0 ; # 1 beta ; 2 occupancy

shake=0; # whether to constrain bonds involving hydrogens

thermostat=1;  # whether to use a thermostat
temperature=300; # kelvin
andersen=0;    # to use Andersen instead of Langevin ; (Note that I see energy up drifts quite often with andersen)
barostat=0;
pressure=1;    # units of atm
membrane_on=0; # whether to use a barostat for membrane simulations (z-axis is the membrane normal)
pme=0; # whether to use PME
pbc=0; # whether periodic boundary conditions are on

dynamo=1
dynamoTemplate='deca.smcv'
watershell_restart='NONE'

mini=0;          # whether to minimize before dynamics
ministeps=0;   # number of minimization iterations

numeq=1             # number of equilibration runs
numeqsteps=10000000; # number of equilibration steps
nummdsteps=10000000; # number of production steps
#nummdsteps=20000
outputfreq=10000;  # frequency of generating output
dcdfreq=10000;     # frequency of dcd output

flag=''
nsteps=numeqsteps;
#
if (firstrun==0):
 restart=0 ; # 0 -- start from PDB coordinates; 1 -- restart from native xml file
 restartfile=None ;
else:
 restart=1 ;
 if (firstrun>1):
  flag=''
 restartfile='./'+name+str(firstrun-1)+flag+'.xml';
#restartfile= ;# to override
 xmlfile=restartfile ;        # to obtain cell vectors from xml file produced with OMM (default option if restart file is provided)
#
# run MD simulations with different parameters one after the other
#
irun=firstrun
while irun < firstrun + numrun :


 print(" =============================");
 print(" Run ", irun, "(will quit after", firstrun+numrun-1,")");
# set some run-specific options
# constraintscaling = (90-10*irun) ;# turn off gradually by run 10
 if irun >= numeq:
   flag=''
   constraints=0 ;# to remove equilibration restraints
   nsteps=nummdsteps ;# increase number of steps
   hmass=1.0
   dt=1.0
   mini=0
   shake=0
   barostat=0 ;# turn off barostat
   pmefreq=1 ; # using MTS does not improve speed in my experience
#   cutoff=9 ;# to change cutoff
#
 if dynamo == 1:
# dynamo section :
  dynamoConfig='deca'+str(irun)+'.in'
# modify config template :
  df=open(dynamoTemplate,'r');
  dd=df.read()
  dd=dd.replace('@{irun}',str(irun))
  dd=dd.replace('@{iprev}',str(irun-1))
  df=open(dynamoConfig,'w');
  df.write(dd);
  df.close();
  dynamoLog=dynamoConfig+'.log';
#
 outputName='./'+name+str(irun)+flag ;
 from os.path import expanduser
# sys.exit()
 exec(open(expanduser('./CHOMM.py')).read())
# to do : check if run was successful, rerun if not
 irun=irun+1
 restartfile=outputName+'.xml'
 xmlfile=restartfile
 restart=1
