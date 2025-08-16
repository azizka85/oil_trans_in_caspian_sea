inline void slae_direct_tridiagonal_calc(
    int ny, int nz, 
    int i, int j,
    __global const float* l,
    __global const float* c,
    __global const float* r,
    __global float* d,
    __global float* u
) {
    int shift = nz * (j + i*ny);

    __global const float* ls = l + shift;
    __global const float* cs = c + shift;
    __global const float* rs = r + shift;
    __global float* ds = d + shift;
    __global float* us = u + shift;

    if (nz > 0) {
        us[0] = rs[0] / cs[0];
        ds[0] = ds[0] / cs[0];

        for (int k = 1; k < nz; k++) {
            float c1 = cs[k] - ls[k]*us[k-1];
    
            us[k] = rs[k] / c1;
            ds[k] = (ds[k] - ls[k]*ds[k-1]) / c1;
        }

        us[nz-1] = ds[nz-1];

        for (int k = nz-2; k >= 0; k--) {
            us[k] = ds[k] - us[k]*us[k+1];
        }
    }
}