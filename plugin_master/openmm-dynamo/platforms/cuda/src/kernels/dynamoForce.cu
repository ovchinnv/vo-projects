extern "C" __global__
void addForces(const real* __restrict__ forces, long long* __restrict__ forceBuffers, int* __restrict__ atomIndex) {
    for (int atom = blockIdx.x*blockDim.x+threadIdx.x; atom < NUM_ATOMS; atom += blockDim.x*gridDim.x) {
        int index = atomIndex[atom];
        forceBuffers[atom] += (long long) (forces[3*index]*0x100000000);
        forceBuffers[atom+PADDED_NUM_ATOMS] += (long long) (forces[3*index+1]*0x100000000);
        forceBuffers[atom+2*PADDED_NUM_ATOMS] += (long long) (forces[3*index+2]*0x100000000);
    }
}
extern "C" __global__
void addForcesSub(const real* __restrict__ forces, long long* __restrict__ forceBuffers, int* __restrict__ atomIndex, int * __restrict__ atomlist) {
    for (int iatom = blockIdx.x*blockDim.x+threadIdx.x; iatom < atomlist[0]; iatom += blockDim.x*gridDim.x) {
        int atom = atomlist[iatom+1]-1; // atomlist stores its length in position 0, which is why we add 1
        int index = atomIndex[atom];
        forceBuffers[atom] += (long long) (forces[3*index]*0x100000000);
        forceBuffers[atom+PADDED_NUM_ATOMS] += (long long) (forces[3*index+1]*0x100000000);
        forceBuffers[atom+2*PADDED_NUM_ATOMS] += (long long) (forces[3*index+2]*0x100000000);
    }
}
