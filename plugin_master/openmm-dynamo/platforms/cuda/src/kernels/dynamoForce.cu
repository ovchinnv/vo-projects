extern "C" __global__
void addForces(const real* __restrict__ forces, long long* __restrict__ forceBuffers, int* __restrict__ atomIndex) {
    for (int atom = blockIdx.x*blockDim.x+threadIdx.x; atom < NUM_ATOMS; atom += blockDim.x*gridDim.x) {
        int index = atomIndex[atom]; // this is the "psf" index of the atom
        forceBuffers[atom] += (long long) (forces[3*index]*0x100000000);
        forceBuffers[atom+PADDED_NUM_ATOMS] += (long long) (forces[3*index+1]*0x100000000);
        forceBuffers[atom+2*PADDED_NUM_ATOMS] += (long long) (forces[3*index+2]*0x100000000);
    }
}
extern "C" __global__
void expandSubForces(const real* __restrict__ subforces, real* __restrict__ forces, const int* __restrict__ atomlist) {
    for (int iatom = blockIdx.x*blockDim.x+threadIdx.x; iatom < atomlist[0]; iatom += blockDim.x*gridDim.x) {
        int atom = atomlist[iatom+1]-1; // psf index of atom
        forces[3*atom  ]=subforces[3*iatom  ];
        forces[3*atom+1]=subforces[3*iatom+1];
        forces[3*atom+2]=subforces[3*iatom+2];
    }
}
