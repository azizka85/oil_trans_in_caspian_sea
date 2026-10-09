#include <opencl/slae/direct/tridiagonal.cl>

__kernel void wind_induced_currents_davies85_variable_parameters_calc_ua(
    double f,
    double kb,
    double g,
    double rho,
    double dx,    
    double dt,
    int nx, int ny, int nz,
    __global const double* h,
    __global const double* qx,
    __global const double* ua,
    __global const double* va,
    __global const double* z,
    __global const double* uf,
    __global double* u1a
) {
    int i = get_global_id(0);
    int j = get_global_id(1);

    int k = j + i*ny;    
    int idx = nz - 1 + nz*k;

    if (h[k] > 0) {
        u1a[k] = ua[k] + f*va[k]*dt - kb*dt*(uf[idx] + ua[k])/h[k] + qx[k]*dt/rho/h[k];

        if (i > 0 && i < nx-1) {
            int kr = j + (i + 1)*ny;
            int kl = j + (i - 1)*ny;


            u1a[k] -= g*dt*(z[kr] - z[kl])/2/dx;
        }
    }
}

__kernel void wind_induced_currents_davies85_variable_parameters_calc_va(
    double f,
    double kb,
    double g,
    double rho,
    double dy,    
    double dt,
    int nx, int ny, int nz,
    __global const double* h,
    __global const double* qy,
    __global const double* ua,
    __global const double* va,
    __global const double* z,
    __global const double* vf,
    __global double* v1a
) {
    int i = get_global_id(0);
    int j = get_global_id(1);

    int k = j + i*ny;    
    int idx = nz - 1 + nz*k;

    if (h[k] > 0) {
        v1a[k] = va[k] - f*ua[k]*dt - kb*dt*(vf[idx] + va[k])/h[k] + qy[k]*dt/rho/h[k];

        if (j > 0 && j < ny-1) {
            int ku = j + 1 + i*ny;
            int kd = j - 1 + i*ny;

            v1a[k] -= g*dt*(z[ku] - z[kd])/2/dy;
        }
    }
}

__kernel void wind_induced_currents_davies85_variable_parameters_calc_rhs(
    double f,
    double kb,
    double rho,
    double dx, double dy,    
    double dt,
    int ny, int nz,
    __global const double* nu,
    __global const double* dz,
    __global const double* h,
    __global const double* qx,
    __global const double* qy,
    __global const double* ua,
    __global const double* u1a,
    __global const double* va,
    __global const double* v1a,
    __global const double* uf,
    __global const double* vf,
    __global double* ud,
    __global double* vd
) {
    int i = get_global_id(0);
    int j = get_global_id(1);
    int k = get_global_id(2);

    int p = j + i*ny;

    int id = k + p*nz;
    
    if (h[p] > 0) {
        int idu = k + 1 + p*nz;
        int idd = k - 1 + p*nz;

        int idb = nz - 1 + p*nz;

        if (k == 0) {
            double um1 = uf[id] + h[p]*qx[p]*dz[k]/rho/nu[id];
            double vm1 = vf[id] + h[p]*qy[p]*dz[k]/rho/nu[id];

            ud[id] = uf[id] + f*vf[id]*dt + dt*(
                nu[idu]*(uf[idu] - uf[id])/(dz[k+1] + dz[k]) - 
                nu[id]*(uf[id] - um1)/2/dz[k]
            )/h[p]/h[p]/dz[k] + kb*(uf[idb] + ua[p])*dt/h[p] - qx[p]*dt/rho/h[p]
            + dt*qx[p]/rho/h[p]/2/dz[k];

            vd[id] = vf[id] - f*uf[id]*dt + dt*(
                nu[idu]*(vf[idu] - vf[id])/(dz[k+1] + dz[k]) - 
                nu[id]*(vf[id] - vm1)/2/dz[k]
            )/h[p]/h[p]/dz[k] + kb*(vf[idb] + va[p])*dt/h[p] - qy[p]*dt/rho/h[p]
            + dt*qy[p]/rho/h[p]/2/dz[k];
        } else if (k == nz-1) {
            double un = (
                uf[id] - kb*h[p]*dz[k]*uf[id]/2/nu[idu] 
                - kb*h[p]*dz[k]*ua[p]/nu[idu]
            )/(1 + kb*h[p]*dz[k]/2/nu[idu]);

            double vn = (
                vf[id] - kb*h[p]*dz[k]*vf[id]/2/nu[idu] 
                - kb*h[p]*dz[k]*va[p]/nu[idu]
            )/(1 + kb*h[p]*dz[k]/2/nu[idu]);

            ud[id] = uf[id] + f*vf[id]*dt + dt*(
                nu[idu]*(un - uf[id])/2/dz[k] - 
                nu[id]*(uf[id] - uf[idd])/(dz[k] + dz[k-1])
            )/h[p]/h[p]/dz[k] + kb*(uf[idb] + ua[p])*dt/h[p] - qx[p]*dt/rho/h[p]
            - dt*kb*ua[p]/h[p]/2/dz[k]/(1 + kb*h[p]*dz[k]/2/nu[idu]);

            vd[id] = vf[id] - f*uf[id]*dt + dt*(
                nu[idu]*(vn - vf[id])/2/dz[k] - 
                nu[id]*(vf[id] - vf[idd])/(dz[k] + dz[k-1])
            )/h[p]/h[p]/dz[k] + kb*(vf[idb] + va[p])*dt/h[p] - qy[p]*dt/rho/h[p]
            - dt*kb*va[p]/h[p]/2/dz[k]/(1 + kb*h[p]*dz[k]/2/nu[idu]);
        } else {
            double u = uf[idb] + ua[p];
            double v = vf[idb] + va[p];

            ud[id] = uf[id] + f*vf[id]*dt + dt*(
                nu[idu]*(uf[idu] - uf[id])/(dz[k+1] + dz[k]) - 
                nu[id]*(uf[id] - uf[idd])/(dz[k] + dz[k-1])
            )/h[p]/h[p]/dz[k] + kb*u*dt/h[p] - qx[p]*dt/rho/h[p];   

            vd[id] = vf[id] - f*uf[id]*dt + dt*(
                nu[idu]*(vf[idu] - vf[id])/(dz[k+1] + dz[k]) - 
                nu[id]*(vf[id] - vf[idd])/(dz[k] + dz[k-1])
            )/h[p]/h[p]/dz[k] + kb*v*dt/h[p] - qy[p]*dt/rho/h[p];
        }
    }
}

__kernel void wind_induced_currents_davies85_variable_parameters_create_tridiagonal_matrix(
    double kb, double dt,
    int ny, int nz,
    __global const double* nu,
    __global const double* dz,
    __global const double* h,
    __global double* al,
    __global double* ac,
    __global double* ar
) {
    int i = get_global_id(0);
    int j = get_global_id(1);
    int k = get_global_id(2);

    int p = j + i*ny;

    int id = k + p*nz;
    int idu = k + 1 + p*nz;

    if (h[p] > 0) {
        if (k == 0) {
            ac[id] = 1 + dt*nu[idu]/h[p]/h[p]/dz[k]/(dz[k+1] + dz[k]);
            ar[id] = -dt*nu[idu]/h[p]/h[p]/dz[k]/(dz[k+1] + dz[k]);
        } else if (k == nz-1) {
            al[id] = -dt*nu[id]/h[p]/h[p]/dz[k]/(dz[k] + dz[k-1]);
            ac[id] = 1 + dt*(
                kb*h[p]*dz[k]/(1 + kb*h[p]*dz[k]/2/nu[idu])/2/dz[k] +
                nu[id]/(dz[k] + dz[k-1])
            )/h[p]/h[p]/dz[k];
        } else {
            al[id] = -dt*nu[id]/h[p]/h[p]/dz[k]/(dz[k] + dz[k-1]);
            ac[id] = 1 + dt*(
                nu[idu]/(dz[k+1] + dz[k]) +
                nu[id]/(dz[k] + dz[k-1])
            )/h[p]/h[p]/dz[k];
            ar[id] = -dt*nu[idu]/h[p]/h[p]/dz[k]/(dz[k+1] + dz[k]);
        }        
    }
}

__kernel void wind_induced_currents_davies85_variable_parameters_calc_uvf(
    int ny, int nz,
    __global const double* h,
    __global const double* al,
    __global const double* ac,
    __global const double* ar,
    __global double* uf,
    __global double* vf,
    __global double* ud,
    __global double* vd
) {
    int i = get_global_id(0);
    int j = get_global_id(1);

    int k = j + i*ny;

    if (h[k] > 0) {
        slae_direct_tridiagonal_calc(
            ny, nz, 
            i, j,
            al, ac, ar, 
            ud, uf
        );

        slae_direct_tridiagonal_calc(
            ny, nz, 
            i, j,
            al, ac, ar, 
            vd, vf
        );    
    }
}

__kernel void wind_induced_currents_davies85_variable_parameters_calc_z(
    double dx, double dy,
    int nx, int ny, int nz,
    __global const double* h,
    __global const double* ua,
    __global const double* va,
    __global double* z
) {
    int i = get_global_id(0);
    int j = get_global_id(1);

    int k = j + i*ny;    

    if (h[k] > 0) {
        if (i > 0 && i < nx-1) {
            int kr = j + (i + 1)*ny;
            int kl = j + (i - 1)*ny;

            z[k] -= (h[kr]*ua[kr] - h[kl]*ua[kl])/2/dx;
        }

        if (j > 0 && j < ny-1) {
            int ku = j + 1 + i*ny;  
            int kd = j - 1 + i*ny;

            z[k] -= (h[ku]*va[ku] - h[kd]*va[kd])/2/dy;
        }
    }
}
