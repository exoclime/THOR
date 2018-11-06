
#include "../headers/phy/valkyrie_conservation.h"

__device__ double atomicAddFunc(double *address, double val) {
    unsigned long long int *address_as_ull = (unsigned long long int *)address;
    unsigned long long int  old            = *address_as_ull, assumed;

    do {
        assumed = old;
        old     = atomicCAS(address_as_ull, assumed, __double_as_longlong(val + __longlong_as_double(assumed)));
        //  Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
    } while (assumed != old);

    return __longlong_as_double(old);
}

__global__ void CalcTotEnergy(double *Etotal_d,
                              double *GlobalE_d,
                              double *Mh_d,
                              double *W_d,
                              double *Rho_d,
                              double *temperature_d,
                              double  Gravit,
                              double  Cp,
                              double  Rd,
                              double  A,
                              double *Altitude_d,
                              double *Altitudeh_d,
                              double *lonlat_d,
                              double *areasT,
                              double *func_r_d,
                              int     num,
                              bool    DeepModel) {

    int id  = blockIdx.x * blockDim.x + threadIdx.x;
    int nv  = gridDim.y;
    int lev = blockIdx.y;

    if (id < num) {
        double Ek, Eint, Eg;
        double wx, wy, wz;
        double Cv = Cp - Rd;

        //calculate control volume
        double zup, zlow, Vol;
        if (DeepModel) {
            zup  = Altitudeh_d[lev + 1] + A;
            zlow = Altitudeh_d[lev] + A;
            Vol  = 0.5 * areasT[id] / A * (zup * zup - zlow * zlow);
        }
        else {
            zup  = Altitudeh_d[lev + 1];
            zlow = Altitudeh_d[lev];
            Vol  = areasT[id] * (zup - zlow);
        }
        //calc cartesian values of vertical wind
        wx = W_d[id * nv + lev] * cos(lonlat_d[id * 2 + 1]) * cos(lonlat_d[id * 2]);
        wy = W_d[id * nv + lev] * cos(lonlat_d[id * 2 + 1]) * sin(lonlat_d[id * 2]);
        wz = W_d[id * nv + lev] * sin(lonlat_d[id * 2 + 1]);

        //kinetic energy density 0.5*rho*v^2
        Ek = 0.5 * ((Mh_d[id * 3 * nv + lev * 3 + 0] + wx) * (Mh_d[id * 3 * nv + lev * 3 + 0] + wx) + (Mh_d[id * 3 * nv + lev * 3 + 1] + wy) * (Mh_d[id * 3 * nv + lev * 3 + 1] + wy) + (Mh_d[id * 3 * nv + lev * 3 + 2] + wz) * (Mh_d[id * 3 * nv + lev * 3 + 2] + wz))
             / Rho_d[id * nv + lev];

        //internal energy rho*Cv*T
        Eint = Cv * temperature_d[id * nv + lev] * Rho_d[id * nv + lev];

        //gravitation potential energy rho*g*altitude (assuming g = constant)
        Eg = Rho_d[id * nv + lev] * Gravit * Altitude_d[lev];

        //total energy in the control volume
        Etotal_d[id * nv + lev] = (Ek + Eint + Eg) * Vol;

#ifdef GLOBAL_CONSERVATION_ATOMICADD
        atomicAddFunc(GlobalE_d, Etotal_d[id * nv + lev]);
#endif // GLOBAL_CONSERVATION_ATOMICADD
        
    // printf("E = %e\n",Etotal_d[id*nv+lev]);
    }
}

__global__ void CalcAngMom(double *AngMomx_d,
                           double *AngMomy_d,
                           double *AngMomz_d,
                           double *GlobalAMx_d,
                           double *GlobalAMy_d,
                           double *GlobalAMz_d,
                           double *Mh_d,
                           double *Rho_d,
                           double  A,
                           double  Omega,
                           double *Altitude_d,
                           double *Altitudeh_d,
                           double *lonlat_d,
                           double *areasT,
                           int     num,
                           bool    DeepModel) {

    int id  = blockIdx.x * blockDim.x + threadIdx.x;
    int nv  = gridDim.y;
    int lev = blockIdx.y;

    if (id < num) {
        double AMx, AMy, AMz;
        double rx, ry, rz, r;

        //calculate control volume
        double zup, zlow, Vol;
        if (DeepModel) {
            zup  = Altitudeh_d[lev + 1] + A;
            zlow = Altitudeh_d[lev] + A;
            Vol  = 0.5 * areasT[id] / A * (zup * zup - zlow * zlow);
            //radius vector
            r  = (A + Altitude_d[lev]);
            rx = r * cos(lonlat_d[id * 2 + 1]) * cos(lonlat_d[id * 2]);
            ry = r * cos(lonlat_d[id * 2 + 1]) * sin(lonlat_d[id * 2]);
            rz = r * sin(lonlat_d[id * 2 + 1]);
        }
        else {
            zup  = Altitudeh_d[lev + 1];
            zlow = Altitudeh_d[lev];
            Vol  = areasT[id] * (zup - zlow);
            //radius vector
            r  = A;
            rx = r * cos(lonlat_d[id * 2 + 1]) * cos(lonlat_d[id * 2]);
            ry = r * cos(lonlat_d[id * 2 + 1]) * sin(lonlat_d[id * 2]);
            rz = r * sin(lonlat_d[id * 2 + 1]);
        }

        //angular momentum r x p (total x and y over globe should ~ 0, z ~ const)
        AMx = ry * Mh_d[id * 3 * nv + lev * 3 + 2] - rz * Mh_d[id * 3 * nv + lev * 3 + 1]
              - Rho_d[id * nv + lev] * Omega * r * rz * cos(lonlat_d[id * 2 + 1]) * cos(lonlat_d[id * 2]);
        AMy = -rx * Mh_d[id * 3 * nv + lev * 3 + 2] + rz * Mh_d[id * 3 * nv + lev * 3 + 0]
              - Rho_d[id * nv + lev] * Omega * r * rz * cos(lonlat_d[id * 2 + 1]) * sin(lonlat_d[id * 2]);
        AMz = rx * Mh_d[id * 3 * nv + lev * 3 + 1] - ry * Mh_d[id * 3 * nv + lev * 3 + 0]
              + Rho_d[id * nv + lev] * Omega * r * r * cos(lonlat_d[id * 2 + 1]) * cos(lonlat_d[id * 2 + 1]);
        //AMx, AMy should go to zero when integrated over globe
        // (but in practice, are just much smaller than AMz)

        //total in control volume
        AngMomx_d[id * nv + lev] = AMx * Vol;
        AngMomy_d[id * nv + lev] = AMy * Vol;
        AngMomz_d[id * nv + lev] = AMz * Vol;

#ifdef GLOBAL_CONSERVATION_ATOMICADD
        atomicAddFunc(GlobalAMx_d, AngMomx_d[id * nv + lev]);
        atomicAddFunc(GlobalAMy_d, AngMomy_d[id * nv + lev]);
        atomicAddFunc(GlobalAMz_d, AngMomz_d[id * nv + lev]);
#endif
    }
}

__global__ void CalcMass(double *Mass_d,
                         double *GlobalMass_d,
                         double *Rho_d,
                         double  A,
                         double *Altitudeh_d,
                         double *lonlat_d,
                         double *areasT,
                         int     num,
                         bool    DeepModel) {

    int id  = blockIdx.x * blockDim.x + threadIdx.x;
    int nv  = gridDim.y;
    int lev = blockIdx.y;

    if (id < num) {
        //calculate control volume
        double zup, zlow, Vol;
        if (DeepModel) {
            zup  = Altitudeh_d[lev + 1] + A;
            zlow = Altitudeh_d[lev] + A;
            Vol  = 0.5 * areasT[id] / A * (zup * zup - zlow * zlow);
        }
        else {
            zup  = Altitudeh_d[lev + 1];
            zlow = Altitudeh_d[lev];
            Vol  = areasT[id] * (zup - zlow);
        }

        //mass in control volume = density*volume
        Mass_d[id * nv + lev] = Rho_d[id * nv + lev] * Vol;

#ifdef GLOBAL_CONSERVATION_ATOMICADD
        atomicAddFunc(GlobalMass_d, Mass_d[id * nv + lev]);
#endif
    }
}
