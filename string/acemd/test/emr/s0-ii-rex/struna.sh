#!/bin/bash
# shell utility functions for ftsm
#======================================================================================================================
 function dpar() { # returns dpar for ftsm ; takes arguments irep(1) brep(2) erep(3) : this, first, and last replica
  irep=$1
  brep=$2
  erep=$3
  if [ $irep -eq $brep ] ; then
   dpar=0
  elif [ $irep -eq $erep ] ; then
   dpar=1
  else
   dpar=0.5
  fi
  echo $dpar
 }
#======================================================================================================================
 function ileft() { # returns left replica index for ftsm ; takes arguments irep(1) brep(2) erep(3) : this, first, and last replica
  irep=$1
  brep=$2
  erep=$3
  if [ $irep -eq $brep ] ; then
   ileft=$irep
  else
   ileft=$((irep-1))
  fi
  echo $ileft
 }
#======================================================================================================================
 function iright() { # returns left replica index for ftsm ; takes arguments irep(1) brep(2) erep(3) : this, first, and last replica
  irep=$1
  brep=$2
  erep=$3
  if [ $erep -eq $irep ] ; then
   iright=$irep
  else
   iright=$((irep+1))
  fi
  echo $iright
 }
#======================================================================================================================
 function update() { # replaces a template parameter in an input file
  f=$1
  p=$2
  if [ -z "$3" ]; then
   v=`eval "echo \\$$p"`
  else
   v=$3
  fi
# print assignment
#  echo ": $p <= $v"
# build sed command
  scmd="s|@{$p}|$v|g"
#  echo $scmd
# execute substitution command
  eval 'sed -i "$scmd" $f'
 }
#======================================================================================================================
 function sm_rex_dir() { # generates randomly dir=0 or 1 for replica exchange direction (i.e. replica n+dir exchanges with replica n+dir+1 ; n=n+2)
  xdir=$((RANDOM%2)) ; #returns either 0 or 1 
  echo $xdir
 }
#======================================================================================================================
 function metropolis() { # Metropolis acceptance criterion $1 = dE (kcal/mol), $2 = temperature (K)
  kboltz=.001987191; # Boltzmann's constant in kcal/mol/K
  dE=$1
  temp=$2
# generate random number
#  random=`echo " scale=10 ; $RANDOM/32767" |bc -l`;
  random=`echo "import random; print random.random() " | python2`
  cmd="if ( $dE < 0 ) { print 1 } else { if ( $random <= e(-($dE)/($kboltz*$temp)) ) { print 1 } else { print 0 } }"
#  echo $cmd >> rex.dbg
  success=`echo $cmd | bc -l`
#  echo $success >> rex.dbg
  echo $success
 }
#======================================================================================================================
 function run_acemd() { #implement as a function to detect crashes
  qok=0;
  while [ 1 ] ; do
   acemd $config > $output
   if [ $? -ne 0 ];  then
    echo " ===> ACEMD execution interrupted due to errors during run. Abort."
    exit
   else
    qunstable=`tail $output | grep -i unstable | wc -m`
    if [ "$qunstable" -gt "0" ]; then
     echo Simulation became unstable. Repeating this run.
     mv $output ${output}.crash
    else
     return
    fi
   fi
  done
 }
 