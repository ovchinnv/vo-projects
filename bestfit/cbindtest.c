// test of bestfit interface
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "bestfit.h"

#define __TOL 1.e-8

int main()
{
// declare data pointers
 __CFLOAT *r=NULL, *r_rot=NULL, *u=NULL, *w=NULL, *eigval=NULL;
 __CFLOAT *ugrad=NULL, *ugrad_fd=NULL, *rtemp=NULL;
 __CINT natom;
 __CFLOAT uplus[9], uminus[9];
 __CBOOL qswapdim=1 ; // if true, coordinates are stored as consecutive x,y,z triplets
 int i, j, n;
 int ix, iy, iz;
 __CFLOAT *rmsd_=NULL;
//
 double rmsd_charmm[31]=
 {0.161762, 0.316171, 0.462415, 0.606288, 0.753354, 0.887121, 0.977275, 1.028353,
  1.058443, 1.080392, 1.099811, 1.118947, 1.138144, 1.158221, 1.179165, 1.200457,
  1.222006, 1.243396, 1.264580, 1.284780, 1.305111, 1.324648, 1.343786, 1.362091,
  1.379977, 1.396627, 1.412621, 1.428078, 1.445545, 1.468143, 1.493650 };
 double rmsd_vmd[31]=
{ 0.16176186501979828, 0.3161707818508148, 0.4624148905277252, 0.6062883138656616,
  0.7533543109893799, 0.8871204853057861, 0.9772745370864868, 1.0283526182174683,
  1.0584430694580078, 1.0803922414779663, 1.0998106002807617, 1.1189470291137695,
  1.1381436586380005, 1.158220648765564, 1.179165005683899, 1.2004566192626953, 
  1.2220057249069214, 1.2433956861495972, 1.2645800113677979, 1.2847803831100464, 
  1.3051105737686157, 1.3246477842330933, 1.3437864780426025, 1.362090826034546, 
  1.379976749420166, 1.396626591682434, 1.4126205444335938, 1.428078055381775, 
  1.445544719696045, 1.4681434631347656, 1.4936503171920776};
//
 puts(" --------------------------------------- ");
 puts(" Reading structures for best-fit test... "); 

 int nrep=32;

 char *basename=(char*)"testdata/diala22_zts_";
 char *ext=(char*)".cor";
 int const slen=100;
 char line[slen], fname[slen];
//
 int size; // global
 for (n=0; n<nrep ;n++) {
  sprintf(line,__CINTFMT, n);
  strcpy(fname,basename);
  strcat(fname,line);
  strcat(fname,ext);
  FILE *fid = fopen(fname, "r");
//
  do {
   fgets(line, slen, fid);
  } while (line[0]=='*');
//
  sscanf(line,__CINTFMT,&natom);
//
  if (n==0) {
    size = 3*natom*sizeof(__CFLOAT);
    r = (__CFLOAT *) malloc(size*nrep);
    u = (__CFLOAT *) malloc(9*nrep*sizeof(__CFLOAT));
    w = (__CFLOAT *) malloc(natom*sizeof(__CFLOAT));
    __CFLOAT d=1.0/natom ; for (i=0;i<natom;i++) { w[i]=d; }
  }
//
  int dr=n*natom*3;
  r_rot=r+dr; // point to current structure
//
  for (i=0;i<3*natom;){
   fscanf(fid,"%*20c");    //skip first 20 chars
   ix=i++; iy=i++; iz=i++;
   fscanf(fid,__CFLOATFMT,r_rot+ix); //read x
   fscanf(fid,__CFLOATFMT,r_rot+iy); //read y
   fscanf(fid,__CFLOATFMT,r_rot+iz); //read z
   fscanf(fid,"%*[^\n]");    //skip to end of line
//printf("R:%20.10f %20.10f %20.10f \n", r[ix], r[iy],r[iz]);
  } //natom
// remove COM
//
//  for (i=0;i<3*natom;){ ix=i++;iy=i++;iz=i++; printf("R:%20.10f %20.10f %20.10f \n", r_rot[ix], r_rot[iy], r_rot[iz]); }
  eigval=(__CFLOAT*)com(r_rot,w,natom,qswapdim);   // pointer returned is of type void, cast to double
//  printf("COM: %20.15f %20.15f %20.15f\n",eigval[0],eigval[1],eigval[2]);
//return;
// subtract COM
  for (i=0;i<3*natom;){r_rot[i++]-=eigval[0];r_rot[i++]-=eigval[1];r_rot[i++]-=eigval[2]; }
 } //nrep
//
//
 puts (" Superposing all structures onto first structure using RMSBestFit ...");
//
 rmsd_=(__CFLOAT*)malloc(nrep*sizeof(__CFLOAT)); rmsd_[0]=0.;
 for (n=1;n<nrep;n++) {
  int dr=n*3*natom ; int du=n*9;
  RMSBestFit ( r + dr , r, w, natom, u + du, qswapdim ) ;
//  RMSBestFitEval( r + dr , r, w, natom, u + du, eigval, qswapdim ) ;
// rotate matrix
  r_rot=(__CFLOAT*)matmul( u+du, r+dr, 3, 3, natom);
//  for (i=0;i<3*natom;){ ix=i++;iy=i++;iz=i++; printf("R:%20.10f %20.10f %20.10f \n", r_rot[ix], r_rot[iy], r_rot[iz]); }
//
  rmsd_[n]=rmsd(r, r_rot, w, natom, qswapdim);
//  printf("%5d Rotation matrix :\n",n);
  int i=du; int j;
//  for (j=0;j<3;j++) { ix=i++;iy=i++;iz=i++;printf("%20.15f %20.15f %20.15f\n", u[ix],u[iy],u[iz]) ;}

// test orthogonality
//  puts(" -------------------------------------------------------------------- ");
//  for
//!  write(0,'(A,I5,3G25.15)') ' Orthogonality error: ',n,rmsd(matmul(u(:,:,n),transpose(u(:,:,n))),Kd,(/one,one,one/))/nine
//!  write(0, '(A,10G25.15)')   ' RMSD from eigenvalues  : ',sqrt(rmsd(r(:,:,n),rave,w)**2+rmsd(r(:,:,1),rave,w)**2-two*sum(eigval)) ! , eigval ! , maxval(abs(rave))
  printf(" RMSD from coordinates  : %25.15f\n",rmsd_[n]);
//return;
//!  write(0, '(A,3G25.15)')   ' RMSD from coordinates  : ',rmsd_(n-1)
//!   write(0, *) rmsd_(n-1), rmsd_charmm(n-1), rmsd_vmd(n-1)
//
 }
// compute differences between this calculation and CHARMM, VMD
 double dcharmm=0., dvmd=0., d=0;
 for (i=1;i<nrep;i++){
  d=fabs(rmsd_[i]-rmsd_charmm[i-1]) ; if ( d  > dcharmm ) dcharmm=d;
  d=fabs(rmsd_[i]-rmsd_vmd[i-1])    ; if ( d  > dvmd    ) dvmd=d;
 }
//
 puts(" -------------------------------------------------------");
 puts("Maximum RMS differences between different calculations:");
 printf("Present vs. CHARMM : %20.15e\n",dcharmm);
 printf("Present vs. VMD    : %20.15e\n",dvmd);
//
// compute gradients of rotation analytically and compare to finite-differences
//
 puts(" -------------------------------------- ");
 puts(" Comparing analytical and finite-difference gradients from RMSBestFit ...");
 puts("  max|u - ( u(x+h) - u(x-h) )/2h|");
 puts("  (NOTE: 2nd order FD should exhibit a quadratic reduction in max. error [for h not too small!])");
 puts;
 puts("  FD interval, Max. grad. error (all structures): ");
 puts("  ----------------------------------------------- ");
//
 ugrad=   (__CFLOAT *) malloc(size*9); //analytical
 ugrad_fd=(__CFLOAT *) malloc(size*9); //FD
//
// initial FD step
 double h=10.0;
//
 double err, maxerr;
 rtemp=(__CFLOAT *) malloc(size);
 for (i=0;i<3*natom;i++){ rtemp[i]=r[i]; } // copy first structure to rtemp 
 do
 {
  maxerr=0.;
  for (n=0;n<nrep;n++) { // over all structures
   int dr=3*natom*n ; int du=9*n;
//   RMSBestFitGrad(r+dr,rtemp,w,natom,u+du,ugrad,1,natom, qswapdim); // compute matrix and gradients
   RMSBestFitGradEval(r+dr,rtemp,w,natom,u+du,ugrad,1,natom, eigval, qswapdim); // compute matrix and gradients; also, eigenvalues
//
   for (i=0;i<3*natom;i++) { // over all atoms
     rtemp[i]=rtemp[i]+h; // fd perturbation
     RMSBestFit(r+dr, rtemp, w, natom, uplus, qswapdim);
     rtemp[i]=rtemp[i]-h*2; // fd perturbation
     RMSBestFit(r+dr, rtemp, w, natom, uminus, qswapdim);
//
     int l;
     for (l=0, iy=i*9 ; l<9 ;l++, iy++) {
      ugrad_fd[iy] = (uplus[l]-uminus[l])*0.5/h ; // second - order
     }
//
     rtemp[i]=r[i]; // restore coordinate
   } // i
// compare analytical and

   err=0;
   for (i=0;i<27*natom;i++){
    d = fabs(ugrad[i]-ugrad_fd[i]) ; if ( d  > err ) err=d;
   }
   if ( err  > maxerr ) maxerr=err;
//
//  printf("  %d,%20.13e\n",n,err);
  } // n-structures
  printf("  %23.13e,%20.13e\n",h,maxerr);
  h=h*0.1;
 } while (h>__TOL);
//
 puts("  --------------------------------- ");
//
// deallocate memory
 free(r);
 free(r_rot);
 free(rtemp);
 free(u);
 free(w);
 free(ugrad);
 free(ugrad_fd);
 free(rmsd_);
 free(eigval);
} //main


