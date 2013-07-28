/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

// written by David Hurwitz, March to May 1998.

#include <memory.h>
#include <iostream>
#include <iomanip>
using namespace std;
#include "charm++.h"
#include "FreeEnergyAssert.h"
#include "FreeEnergyGroup.h"


AGroup::AGroup() {
//------------------------------------------------------------------------
// make room for some ints.
//------------------------------------------------------------------------
  size = 0;
  maxSize = kGroupNumToStart;
  inds = new int[maxSize];
  weights = new double[maxSize];
  normWeights = new double[maxSize];
  r = new double[3*maxSize];
//
  std::fill_n(r, 3*maxSize, 0.); // initialize to zero
//
  isRset = 1;
  sumWeights = 0.0;
  normalizeWeights=1;
}


AGroup::~AGroup() {
//------------------------------------------------------------------------
// return borrowed memory to the free store.
//------------------------------------------------------------------------
  delete []inds;
  delete []weights;
  delete []normWeights;
  delete []r;
}


void AGroup::Clear() {
//------------------------------------------------------------------------
// leave memory allocation alone.
//------------------------------------------------------------------------
  size = 0;
  sumWeights=0.0;
  normalizeWeights=1;
  isRset=1;
}


void AGroup::Add(int i, double w, double x, double y, double z) {
//------------------------------------------------------------------------
// add an int to the list.  if there's not enough room, make room.
//------------------------------------------------------------------------
  int*  new_inds;
  double* new_weights;
  double* new_normWeights;
  double* new_r;

  // if there's no room for a new int
  if (size == maxSize) {
    // create an array with more space
    maxSize *= kGroupMultiplier;
    new_inds = new int[maxSize];
    new_weights = new double[maxSize];
    new_normWeights = new double[maxSize];
    new_r = new double[3*maxSize];
    // fast copy from the full array to the new one (memcpy(dest, src, bytes))
    memcpy(new_inds, inds, sizeof(int)*size);
    memcpy(new_weights, weights, sizeof(double)*size);
    memcpy(new_r, r, sizeof(double)*size*3);
    // no need to copy normalized weights, since they will be recomputed
    // return the space used for the full array
    delete []inds;
    delete []weights;
    delete []normWeights;
    delete []r;
    // point to the bigger arrays
    inds = new_inds;
    weights = new_weights;
    normWeights = new_normWeights;
    r = new_r;
  }

  // add the int to the int array, etc.
  inds[size] = i;
  weights[size] = w;
  int j = size * 3;
  r[j++]=x;  r[j++]=y; r[j]=z;
  size++;
  sumWeights += w;
  normalizeWeights=1;
}


int AGroup::operator[] (int Index) {
//------------------------------------------------------------------------
// return an int from this group of ints.
// note, by returning int, rather than int&, this function only allows
// retrieval of an int, not setting an int.  to set, use add.
//------------------------------------------------------------------------
  ASSERT((Index>=0) && (Index<size));
  return(inds[Index]);
}

AGroup& AGroup::operator= (AGroup& Group) {
//------------------------------------------------------------------------
// make this object identical to the passed one.
//------------------------------------------------------------------------
  // if there's not enough room here for Group's array
  if (size < Group.size) {
    // free this array and allocate space for the bigger one
    delete []inds;
    delete []r;
    delete []weights;
    delete []normWeights;
    maxSize = Group.maxSize;
    inds = new int[maxSize];
    weights = new double[maxSize];
    r = new double[maxSize*3];
    normWeights = new double[maxSize];
  }
  // fast copy from the passed array to this one (memcpy(dest, src, bytes))
  size = Group.size;
  sumWeights=Group.sumWeights;
  normalizeWeights=Group.normalizeWeights;
  memcpy(inds, Group.inds, sizeof(int)*size);
  memcpy(weights, Group.weights, sizeof(double)*size);
  memcpy(r, Group.r, sizeof(double)*size*3);
  memcpy(normWeights, Group.normWeights, sizeof(double)*size);
  return(*this);
}
//
void AGroup::List(int NumToList) {
//------------------------------------------------------------------------
// list NumToList integers in the group to standard out.
// if NumToList is negative, list them all.
//------------------------------------------------------------------------
  int  i;

  if (NumToList < 0) {
    NumToList = size;
  }
  for (i=0; i<NumToList; i++) {
//    cout << setw(10) << i << "   " << setw(10) << (*this)[i] << std::endl;
    CkPrintf("%10d   %10d\n", i, (*this)[i]);
  }
}

double const * AGroup::GetWeights() {
  if (normalizeWeights) { //normalize
//normalize weights
   double omass = ( sumWeights > 0.0 ? 1.0/sumWeights : 1.0 ) ;
   for (int i=0 ; i < size ; i++ ) { normWeights[i] = omass * weights[i]; }
   normalizeWeights=0;
  }
  return (double const *) normWeights;
}

double * AGroup::GetCoords() {
 return r;
}

int const * AGroup::GetInds() {
 return (int const *) inds;
}

void AGroup::SetCoords(int i, double* xyz) {
 ASSERT( i > -1 && i < size );
 int j = i*3;
 r[j++]=xyz[0];
 r[j++]=xyz[1];
 r[j++]=xyz[2];
}

void AGroup::SetCoords(int i, double x, double y, double z) {
 ASSERT( i > -1 && i < size );
 int j = i*3;
 r[j++]=x;
 r[j++]=y;
 r[j++]=z;
}
