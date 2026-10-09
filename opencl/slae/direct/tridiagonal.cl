inline void slae_direct_tridiagonal_calc(
    int ny, int nz, 
    int i, int j,
    __global const double* l,
    __global const double* c,
    __global const double* r,
    __global double* d,
    __global double* u
) {
    int shift = nz * (j + i*ny);

    __global const double* ls = l + shift;
    __global const double* cs = c + shift;
    __global const double* rs = r + shift;
    __global double* ds = d + shift;
    __global double* us = u + shift;

    if (nz > 0) {
        us[0] = rs[0] / cs[0];
        ds[0] = ds[0] / cs[0];

        for (int k = 1; k < nz; k++) {
            double c1 = cs[k] - ls[k]*us[k-1];
    
            us[k] = rs[k] / c1;
            ds[k] = (ds[k] - ls[k]*ds[k-1]) / c1;
        }

        us[nz-1] = ds[nz-1];

        for (int k = nz-2; k >= 0; k--) {
            us[k] = ds[k] - us[k]*us[k+1];
        }
    }
}