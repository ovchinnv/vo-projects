/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

//--------------------------------------------------------------------
// AGroup contains an indexed list of atomic coordinates and weights
// written by David Hurwitz, March to May 1998.
// modifications by VO, 2013
//--------------------------------------------------------------------
#if !defined(GROUP_HPP)
  #define GROUP_HPP

const int kGroupNumToStart = 16;    // to start, there's room for this num ints.
const int kGroupMultiplier = 4;     // each time array size is exceeded,
                                    // its size is increased by this many times.

class AGroup {
private:
  int*  inds;      // list of integer indices
  int   size;      // the number of integers in the list
  int   maxSize;   // the maximum number of integers allowed in
                   // the list without allocating more memory
// VO 2013

  double* weights;     // list of weights corresponding to the indices
  double* normWeights ;// list of normalized weights
  double  sumWeights;  // sum of weights
//
  double* r;    // coordinate triplets
  bool isRset;  // flag that determines whether coordinates have been set
//

  bool normalizeWeights;
public:
  AGroup();
  ~AGroup();
  void    Clear();
  void    Add(int index, double weight, double x, double y, double z);
  void    Add(int index, double weight=1.0) { Add(index, weight, 0.0, 0.0, 0.0) ; isRset=0; }
  AGroup& operator= (AGroup& Group);
  int     operator[] (int Index);
  int     GetSize() { return(size); }
  void    List(int NumToList=-1);
  double const * GetWeights();
  double const * GetCoords();
  void SetCoords(int, double*);
  void SetCoords(int, double, double, double);
  int const * GetInds();
  bool HaveCoords() {return isRset;}
};

#endif
