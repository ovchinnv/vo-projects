#!/bin/python
#
import simtk.openmm.app as app
import simtk.openmm as mm
import simtk.unit as u
from sys import stdout, stderr, exit
from shutil import move

#=====================================================================#
# define parameters that may not have been defined by user
#
if 1:
 try :
  corfile
 except NameError:
  corfile=None
#
 try :
  outputName
 except NameError:
  outputName='output'
#
 try :
  mini
 except NameError:
  mini=0
#
 try :
  constraints
 except NameError:
  constraints=0
 try :
  conscol
 except NameError:
  conscol=1
#
 try :
  switchdist
 except NameError:
  switchdist=cutoff-1.5;
#
 try :
  pbc
 except NameError:
  pbc=0
#
 try :
  pme
 except NameError:
  pme=0
#
 try :
  thermostat
 except NameError:
  thermostat=0
 try :
  andersen
 except NameError:
  andersen=0
#
 try :
  barostat
 except NameError:
  barostat=0
 try :
  membrane_on
 except NameError:
  membrane_on=0
#
 try :
  implicitSolvent
 except NameError:
  implicitSolvent=0
# string method
 try :
  struna
 except NameError:
  struna=0
# tempering plugin
 try :
  adTemp
 except NameError:
  adTemp=0
#
 try :
  platformName
 except NameError:
# use CUDA unless variable 'platformName' defined
  platformName="CUDA"
#
#========================== Subroutines
#==========================
 def dprint(*args):
  print(" ===> CHOMMPy : ",end="");
  for arg in args:
   print(arg,end="")
  print(""); # flush
#==========================
 def derror(*args):
  print(" ===> CHOMMPy ERROR : ",end="");
  for arg in args:
   print(arg,end="")
  print(""); # flush
#==========================
 def printe(simulation):
  forceGroups={'Bond':0, 'Angle':1, 'Dihed':2, 'UB':3, 'IMPR':4, 'CMAP':5, 'NBOND':6};
  ener={};
# evaluate
  for key in forceGroups:
   state=simulation.context.getState(getEnergy=True, groups=1<<forceGroups[key]) ;
   ener[key]=state.getPotentialEnergy().value_in_unit(u.kilocalories_per_mole);
# total potential energy
  ener['PE']=simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(u.kilocalories_per_mole);
# print
  for key in ['Bond', 'Angle', 'Dihed', 'UB', 'IMPR', 'CMAP', 'NBOND', 'PE']:
   print(key, end="\t\t\t");
  print(); # newline
  for key in ['Bond', 'Angle', 'Dihed', 'UB', 'IMPR', 'CMAP', 'NBOND', 'PE']:
   print(ener[key], end="\t");
  print();
#==========================
 def pote(simulation):
  "Return potential energy"
  pener=simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(u.kilocalories_per_mole);
  return pener;
#========================== box dimensions from .str file produced during system preparation
 def get_box_size_str(boxfile):
  with open(boxfile) as f:
     lines=f.readlines()
     i=0;
     for line in lines:
      words = line.split()
      if (words[0].upper() == 'SET') :
       break;
      i=i+1;
    # get box size
     dx = float(lines[i].split()[2]);
     dy = float(lines[i+1].split()[2]);
     dz = float(lines[i+2].split()[2]);
  return(dx, dy, dz);
#========================== box dimensions from .xml restart file produced with openmm
 def get_box_size_xml(xmlfile):
  with open(xmlfile) as f:
     lines=f.readlines()
     i=1;
     for line in lines:
      words = line.split()
      if (words[0].upper() == '<PERIODICBOXVECTORS>') :
       break;
      i=i+1;
    # get  box size
     dx = (float(lines[i].split('"')[1]))*u.nanometers.conversion_factor_to(u.angstrom)  ; # remember that these are in nanoneters
     dy = (float(lines[i+1].split('"')[3]))*u.nanometers.conversion_factor_to(u.angstrom) ;
     dz = (float(lines[i+2].split('"')[5]))*u.nanometers.conversion_factor_to(u.angstrom) ;
  return(dx, dy, dz);
#========================== box dimensions from .xsc file produced with ACEMD/NAMD programs
 def get_box_size_xsc(xscfile):
  with open(xscfile) as f:
     lines=f.readlines()
     for line in lines:
      pass
     cell=line.split(' ')
     dx = float(cell[1]);
     dy = float(cell[5]);
     dz = float(cell[9]);
  return(dx, dy, dz);
#
#========================== Initialize simulation system
 dprint("Reading PSF from file '", psffile, "'");
 psf=app.CharmmPsfFile(psffile);
 dprint("Reading topology from file '",topfile,"' and parameters from file '", paramfile,"'");
 params=app.CharmmParameterSet(topfile, paramfile, permissive=False); # running without atom typing (via mass entries in the topology) often leads to ERRORS !
#========================================================
 if (pbc):
  dprint("Periodic boundary conditions will be used")
  if (not restart):
   try :
    dx; dy; dz; # check if dimensions are specified manually
   except NameError:
    try :
     boxfile;
     dprint("Setting orthorhombic cell lengths from file '",boxfile,"'")
     dx, dy, dz=get_box_size_str(boxfile);
    except NameError:
     try:
      xscfile;
      dprint("Setting orthorhombic cell lengths from file '",xscfile,"'")
      dx, dy, dz=get_box_size_xsc(xscfile);
     except NameError:
      try:
       xmlfile;
       dprint("Setting orthorhombic cell lengths from file '",xmlfile,"'")
       dx, dy, dz=get_box_size_xml(xmlfile);
      except NameError:
       derror("Could not set periodic cell size.")
  else:
   try:
    xmlfile;
    dprint("Setting orthorhombic cell lengths from file '",xmlfile,"'")
    dx, dy, dz=get_box_size_xml(xmlfile);
   except Die:
    derror("Could not set periodic cell size.")
#
  try:
   dprint("Periodic cell dimensions are (", dx*u.angstrom, ")x(", dy*u.angstrom, ")x(", dz*u.angstrom,")")
   psf.setBox(dx*u.angstrom, dy*u.angstrom, dz*u.angstrom);
  except Exception:
   sys.exit(-1)
#===========================================================
  if (pme):
   nbondMethod=app.PME
   dprint("PME is on");
  else:
   nbondMethod=app.CutoffPeriodic
   dprint("PME is off");
 else:
  psf.setBox(1000*u.angstrom, 1000*u.angstrom, 1000*u.angstrom) # set to a very large box to eliminate wrapping
  nbondMethod=app.CutoffNonPeriodic
#===================================================== SHAKE
 if (shake==1):
  cons=app.HBonds
  rigidWater=True
  dprint("Will constrain all bonds involving hydrogens");
 elif (shake>1):
  cons=app.AllBonds
  rigidWater=True
  dprint("Will constrain all bond lengths");
 else:
  cons=None
  rigidWater=False
 dprint("Initializing simulation system");
 dprint("Nonbonded cutoff is ",cutoff*u.angstrom,". Switching is active at ",switchdist*u.angstrom)
 if (hmass>1):
  dprint("Hydrogen mass is ",hmass*u.amu)

 if (implicitSolvent==1):
  system=psf.createSystem(params,
                         nonbondedMethod=nbondMethod, nonbondedCutoff=cutoff*u.angstrom, switchDistance=switchdist*u.angstrom,
                         constraints=cons, rigidWater=rigidWater, removeCMMotion=False, hydrogenMass=hmass*u.amu,
                         implicitSolvent=app.OBC2,
                         verbose=False);
 else:
  system=psf.createSystem(params,
                         nonbondedMethod=nbondMethod, nonbondedCutoff=cutoff*u.angstrom, switchDistance=switchdist*u.angstrom,
                         constraints=cons, removeCMMotion=False, hydrogenMass=hmass*u.amu, rigidWater=rigidWater,
                         verbose=False);

#================= harmonic restraints from file, a la NAMD/ACEMD
 if (constraints) :
  if (conscol==1): # beta
   dprint("Adding absolute positional harmonic restraints to atoms marked in the beta column of PDB file '"+consfile+"'");
  elif (conscol==2): #occupancy
   dprint("Adding absolute positional harmonic restraints to atoms marked in the occupancy column of PDB file '"+consfile+"'");

  force=mm.CustomExternalForce("s*0.5*k*periodicdistance(x,y,z,x0,y0,z0)^2");
#  force=mm.CustomExternalForce("s*0.5*k*( (x-x0)^2 + (y-y0)^2 + (z-z0)^2 )");
  force.addPerParticleParameter("k");
  force.addPerParticleParameter("x0");
  force.addPerParticleParameter("y0");
  force.addPerParticleParameter("z0");
  force.addGlobalParameter("s", constraintscaling);
# read per atom restraints :
  res=app.PDBFile(consfile);
  iatom=0; icons=0;
  for r, o, b  in zip(res.positions, res.occupancy, res.temperature_factor) :
   if (conscol==1): # beta
    bnodim=b/u.angstrom/u.angstrom; # have to deal with units, which are A^2 for B-factors
   elif (conscol==2): #occupancy
    bnodim=o

   if (bnodim > 0) :
    icons+=1;
    k=bnodim*u.kilocalorie/u.mole/u.angstrom/u.angstrom
#   dprint(" Adding restraint on atom ",iatom," with force constant ", k );
    x0=r[0].value_in_unit(u.nanometer)
    y0=r[1].value_in_unit(u.nanometer)
    z0=r[2].value_in_unit(u.nanometer)
    force.addParticle(iatom, [k,x0,y0,z0]);
   iatom+=1;
  dprint("Added restraints on ", icons, " atoms");
  dprint("Harmonic force constants will be scaled uniformly by x"+str(constraintscaling));
  system.addForce(force)
#
#================= string plugin
 if (struna==1) :
  from openmmstruna import *
  system.addForce(StrunaForce(strunaConfig, strunaLog))
#================= add integrator :
 dprint("Configuring integrator");
# first, add barostat if requested :
 if (thermostat and barostat):
  if (membrane_on):
   dprint("Initializing Monte-Carlo membrane barostat at pressure ",pressure*u.atmosphere, " with no surface tension in the plane of membrane (XY)");
   barostatForce=mm.MonteCarloMembraneBarostat(pressure*u.atmosphere, 0*u.atmosphere*u.angstrom, temperature*u.kelvin,
                                               mm.MonteCarloMembraneBarostat.XYIsotropic, mm.MonteCarloMembraneBarostat.ZFree);
  else:
   dprint("Initializing Monte-Carlo barostat at pressure ",pressure*u.atmosphere)
   barostatForce=mm.MonteCarloBarostat(pressure*u.atmosphere, temperature*u.kelvin)
  system.addForce(barostatForce) ;
#
#================= deal with multiple timestepping :
 if (pme and pmefreq > 1) :
#================= add thermostat force if needed
  if (thermostat):
   dprint("Initializing Andersen thermostat with coupling to bath with friction ",friction/u.picosecond," at temperature ",temperature*u.kelvin)
   thermostatForce=mm.AndersenThermostat(temperature*u.kelvin, friction/u.picosecond);
#   thermostatForce.setForceGroup(0) ; # make sure to assign a group for MTS
   system.addForce(thermostatForce);

  dprint("Multiple time stepping (MTS) will be used");
  dprint("All forces except nonbonded reciprocal forces will be assigned a ForceGroup of zero");
# split force objects into a multiple groups
# note that createSystem puts different psf sections into different force groups for ease of energy decomposition;
  for f in system.getForces() :
# put all forces into the same group (note that this will make energy decomposition impossible)
# we should be able to use many groups with the same substep in the RESPA init, but that might slow it down
   f.setForceGroup(0);
# reciprocal forces get a separate group for RESPA
   if isinstance(f,mm.NonbondedForce) :
      f.setReciprocalSpaceForceGroup(31);
  dprint("Initializing RESPA MTS integrator");
  integrator=mm.MTSIntegrator(dt*pmefreq*u.femtosecond, [(31,1), (0,pmefreq)]);
 else :
  if (thermostat):
   if (andersen):
    dprint("Initializing Andersen thermostat coupled to bath with friction ",friction/u.picosecond," at temperature ",temperature*u.kelvin);
    system.addForce(mm.AndersenThermostat(temperature*u.kelvin, friction/u.picosecond));
    dprint("Initializing Verlet integrator with timestep ",dt*u.femtosecond);
    integrator=mm.VerletIntegrator(dt*u.femtosecond);
   else:
    dprint("Initializing Langevin thermostatted integrator with timestep ",dt*u.femtosecond," coupled to bath with friction ",friction/u.picosecond," at temperature ",temperature*u.kelvin);
    integrator=mm.LangevinIntegrator(temperature*u.kelvin, friction/u.picosecond, dt*u.femtosecond);
  else:
   dprint("Initializing Verlet integrator with timestep ",dt*u.femtosecond);
   integrator=mm.VerletIntegrator(dt*u.femtosecond);
#====================================================
 if (adTemp):
  import tempering
  from math import ceil,  sqrt
  from numpy import array as narray
  dprint("Initializing Adaptive Tempering Plugin with input file '",adTempConfig,"' and output file '",adTempLog,"'");
  tempering.init(adTempConfig, adTempLog)
  adTemperature=temperature; # initial temperature, must match config file
#====================================================
#
 dprint("Initializing compute platform ",platformName);
 platform=mm.Platform.getPlatformByName(platformName);
 properties={'CudaPrecision': 'mixed'};
 dprint("Preparing simulation topology");
 if (platformName=="CUDA") :
  simulation=app.Simulation(psf.topology, system, integrator, platform, properties);
 else :
  simulation=app.Simulation(psf.topology, system, integrator, platform);
#
 if (restart == 0) :
  if (corfile!=None):
   cor=app.CharmmCrdFile(corfile);
   dprint("Setting simulation coordinates from file '",corfile,"'");
   simulation.context.setPositions(cor.positions);
  else:
   pdb=app.PDBFile(pdbfile);
   dprint("Setting simulation coordinates from file '",pdbfile,"'");
   simulation.context.setPositions(pdb.positions);
 else :
  dprint("Setting simulation restart data from file '",restartfile,"'");
  with open(restartfile, 'r') as f:
   xml=f.read();
   oldstate=mm.XmlSerializer.deserialize(xml)
   simulation.context.setPositions(oldstate.getPositions());
   simulation.context.setVelocities(oldstate.getVelocities());
   simulation.context.setTime(oldstate.getTime());
#
#================ Print initial energy compoments :
 dprint("Initial Potential energy" );
 printe(simulation);
#================ Energy minimization
# NOTE : I have been unable to use the minimizer when both maxIterations and tolerance are specified (OpenMM7)
 if (mini) :
  dprint("Minimizing energy for ",ministeps," steps");
  simulation.minimizeEnergy(maxIterations=ministeps); # optional iterations, tolerance
  dprint("Potential energy after minimization");
  printe(simulation);
#=============== MD simulation
 if (nsteps>0):
  simulation.reporters.append(app.DCDReporter('output.dcd',dcdfreq));
  simulation.reporters.append(app.StateDataReporter(stdout, outputfreq, step=True, potentialEnergy=True, kineticEnergy=True, speed=True, temperature=True, 
                                                    volume=pbc, separator=' \t '));
  if (adTemp):
   num_inner_steps=adTempFreq ;
   num_outer_steps=ceil(nsteps/num_inner_steps) ;
   numsteps=num_inner_steps*num_outer_steps ;
   dprint("Running MD simulation for ",numsteps," steps");
   for n_outer in range(0,num_outer_steps):
    simulation.step(num_inner_steps); # run system
# compute potential energy
    state=simulation.context.getState(getEnergy=True)
    pener=state.getPotentialEnergy().value_in_unit(u.kilocalories_per_mole);
    newAdTemperature=tempering.update(simulation.currentStep, pener, adTemperature);
    if (abs(newAdTemperature-adTemperature) > 1.0e-10) :
     state=simulation.context.getState(getVelocities=True)
#     dprint("Updating thermostat temperature ",adTemperature," --> ",newAdTemperature);
     velocity_scale = sqrt(newAdTemperature/adTemperature);
     adTemperature=newAdTemperature;
     simulation.integrator.setTemperature(adTemperature);
     velocities=state.getVelocities().value_in_unit(u.nanometers/u.picoseconds);
     new_velocities=[v*velocity_scale for v in velocities];
     simulation.context.setVelocities(new_velocities);
  else:
   dprint("Running MD simulation for ",nsteps," steps");
   simulation.step(nsteps);


  dprint("Writing simulation restart files");
  simulation.saveState(outputName+'.xml');
  simulation.saveCheckpoint(outputName+'.chk');
#==== write periodic box vectors
  state=simulation.context.getState();
  a,b,c=state.getPeriodicBoxVectors();
  fxsc=open(outputName+'.xsc','w');
  fxsc.write("#CHOMMPy.xsc stub\n");
  fxsc.write(str(nsteps)+" "+str(a[0].value_in_unit(u.angstrom))+" 0 0 0 "+str(b[1].value_in_unit(u.angstrom))+" 0 0 0 "+str(c[2].value_in_unit(u.angstrom))+" 0 0 0 0 0 0 0 0 0\n");
  fxsc.close();

#==== reset switching distance
 del switchdist;
#==== move dcd file to destination file
 move('output.dcd', outputName+'.dcd');
#==== finalize tempering plugin
 if (adTemp):
  dprint("Finalizing Adaptive Tempering Plugin")
  tempering.done()
