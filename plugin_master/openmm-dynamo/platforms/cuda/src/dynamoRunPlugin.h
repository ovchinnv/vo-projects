// can be multiply included from CudaDynamoKernel.cpp for different precision
    if (atomlist==NULL) { // atomlist is not defined; therefore, provide all coords
     for (i=0, rptr=r ; i < natoms ; i++) {
      *(rptr) = pos[i][0]*nm2A; rptr++ ; //units
      *(rptr) = pos[i][1]*nm2A; rptr++ ;
      *(rptr) = pos[i][2]*nm2A; rptr++ ;
// cerr << "position of atom: "<<j<<"="<<pos[j][0]<<pos[j][1]<<pos[j][2]<<endl;
     }
    // compute plugin forces and energy
// cerr << "CALLING DYNAMO"<<endl;
     ierr = (sizeof(_FLOAT)==sizeof(double)) ? \
      master_dyna_plugin(iteration, 0, (double*)r, NULL, (double*)fr, NULL, (double*)&master_energy, NULL, &atomlist, usesPeriodic, (double*)&box, NULL) : \
      master_dyna_plugin(iteration, 1, NULL, (float*)r,  NULL, (float*)fr, NULL,  (float*)&master_energy, &atomlist, usesPeriodic, NULL, (float*)&box) ; // might return valid atomlist
    // copy plugin forces
     if (atomlist!=NULL) { // atom indices provided; use them for adding forces
// if we are here, this means that the atom indices were not provided in the initialization call, so we need to perform that part of init here
#ifdef __DYNAMO_SUBSET
// allocate & populate device atomlist
       natoms_requested=atomlist[0];
       cudaError_t cuda_err=cudaMalloc((void**)&(atomlist_d),(natoms_requested+1)*sizeof(int));
       if (cuda_err!=cudaSuccess) {ierr=1;cerr<<__STRNG(_whoami)<<"Error allocating CUDA memory:"<<cudaGetErrorString(cuda_err);}
       cudaMemcpy(atomlist_d,atomlist,(natoms_requested+1)*sizeof(int),cudaMemcpyHostToDevice);
//       memset(frc,0,3*natoms_requested*elementSize); //should not be needed
// copy subset of forces to array
       for ( ii=0, i=0 ; i++ < natoms_requested; ) { // increment i right after the comparison, b/c skipping 0th entry which stores list size
        j=3*(atomlist[i]-1);
        frc[ii]= fr[j]*str2omm_f; ++ii; ++j;
        frc[ii]= fr[j]*str2omm_f; ++ii; ++j;
        frc[ii]= fr[j]*str2omm_f; ++ii;
       }
#else
       memset(frc,0,3*natoms_requested*elementSize); // make sure all forces are zero, since we will be using an all-atom force assignment kernel
// NOTE: even though we communicated the entire atom array, we still only use the indices in the atomlist, hopefully with slight time savings
//       for (aptr=atomlist+1 ; aptr<=atomlist + atomlist[0] ; aptr++) // iterate until atomlist points to the last index
//        i=*aptr - 1; // for zero offset (first coordinate lives in r[0])
//        j=3*i;
       for (i=0 ; i++ < atomlist[0] ; ){ // same as above, but clearer ; note: comparing i, then immediately incrementing
        j=3*(atomlist[i]-1);
        frc[j]= fr[j]*str2omm_f;++j;
        frc[j]= fr[j]*str2omm_f;++j;
        frc[j]= fr[j]*str2omm_f;
       }
#endif
     } else { // no atomlist provided; loop over all atoms (natoms_requested=natoms)
       for (j=0 ; j < 3*natoms_requested ; j++) {
        frc[j]= fr[j]*str2omm_f;
       }
     } // atomlist
//============================================================================================================================
    } else { // atomlist not null : loop over only the desired indices
     for (aptr=atomlist+1 ; aptr<atomlist + 1 + (*atomlist) ; aptr++) { // iterate until atomlist points to the last index
      j=*aptr - 1;
      rptr=r + 3*j ;
      *(rptr) = pos[j][0]*nm2A; rptr++ ;//units
      *(rptr) = pos[j][1]*nm2A; rptr++ ;
      *(rptr) = pos[j][2]*nm2A;
// cerr << "position of atom: "<<j<<"="<<pos[j][0]<<pos[j][1]<<pos[j][2]<<endl;
     }
#ifdef __TIMEDEBUG
     clock_t dstart_time=clock();
#endif
     ierr = (sizeof(_FLOAT)==sizeof(double)) ? \
      master_dyna_plugin(iteration, 0, (double*)r, NULL, (double*)fr, NULL, (double*)&master_energy, NULL, &atomlist, usesPeriodic, (double*)&box, NULL) : \
      master_dyna_plugin(iteration, 1, NULL, (float*)r,  NULL, (float*)fr, NULL, (float*)&master_energy, &atomlist, usesPeriodic, NULL, (float*)&box) ; // atomlist should not be modified in this call
#ifdef __TIMEDEBUG
     clock_t dstop_time=clock();
     cout<<__STRNG(_whoami)<<"master_dyna_plugin took "<<dstop_time-dstart_time<<" cycles"<<endl;
     cout<<__STRNG(_whoami)<<"OMP INFO :"<<endl;
#endif
//
#ifdef __DYNAMO_SUBSET
//     memset(frc,0,3*natoms_requested*elementSize); // should not be needed
// copy forces to subset array
     for ( ii=0, i=0 ; i++ < natoms_requested; ) {
        j=3*(atomlist[i]-1);
        frc[ii]= fr[j]*str2omm_f; ++ii;++j;
        frc[ii]= fr[j]*str2omm_f; ++ii;++j;
        frc[ii]= fr[j]*str2omm_f; ++ii;
     }
#else
// copy forces of select atoms to full array
//     for (aptr=atomlist+1 ; aptr<atomlist + 1 + (*atomlist) ; aptr++)  // iterate until atomlist points to the last index
//        i=*aptr - 1; // for zero offset (first coordinate lives in r[0])
//        j=3*i;
     memset(frc,0,3*natoms_requested*elementSize);
     for (i=0 ; i++ < atomlist[0] ; ) {
        j=3*(atomlist[i]-1);
        frc[j]= fr[j]*str2omm_f;++j;
        frc[j]= fr[j]*str2omm_f;++j;
        frc[j]= fr[j]*str2omm_f;
     }
#endif
    } // atomlist == NULL
