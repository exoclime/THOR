// ==============================================================================
// This file is part of THOR.
//
//     THOR is free software : you can redistribute it and / or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation, either version 3 of the License, or
//     (at your option) any later version.
//
//     THOR is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
//     GNU General Public License for more details.
//
//     You find a copy of the GNU General Public License in the main
//     THOR directory under <license.txt>.If not, see
//     <http://www.gnu.org/licenses/>.
// ==============================================================================
//
//
//
//
// Description: Sets up mid latitude jet to test steady state capability of core
//
//
// Method:
//
// Known limitations: None.
//
// Known issues: None.
//
// If you use this code please cite the following references:
//
//       [1] Mendonca, J.M., Grimm, S.L., Grosheintz, L., & Heng, K., ApJ, 829, 115, 2016
//
//       [2] Held, I. M., & Suarez, M. J. 1994, Bullentin of the American
//           Meteorological Society
//
// Current Code Owner: Joao Mendonca, EEG. joao.mendonca@csh.unibe.ch
//
// History:
// Version Date       Comment
// ======= ====       =======
//
// 1.0     16/08/2017 Released version  (JM)
//
////////////////////////////////////////////////////////////////////////

__device__ double dTemp_dphi(double U,
                             double H,
                             double Rd,
                             double Omega,
                             double A,
                             double lat,
                             double Fz,
                             double dFzdz) {
// Calculate the meridional temperature gradient at latitude 'lat'
      double dU0dz, U0, fCori;

      dU0dz = U*pow(sin(M_PI*pow(sin(lat),2)),3)*dFzdz;
      U0 = U*pow(sin(M_PI*pow(sin(lat),2)),3)*Fz;
      fCori = 2*Omega*sin(lat);

      return -H/Rd*(A*fCori+2*U0*tan(lat))*dU0dz;
}

__global__ void setup_jet(double *Mh_d         ,
                            double *pressure_d   ,
                            double *Rho_d        ,
                            double *temperature_d,
                            double  Cp           ,
                            double  Rd           ,
                            double  Omega        ,
                            double  A            ,
                            double *Altitude_d   ,
                            double *lonlat_d     ,
                            int     num          ){

    int id = blockIdx.x * blockDim.x + threadIdx.x;
    int nv = gridDim.y;
    int lev = blockIdx.y;

    if (id < num){
        double z0 = 1823e3, z1 = 2486e3, dz0 = 414e3;
        double z = Altitude_d[id*nv+lev];
        double pref = 1e5;
        double U = 500, U0;
        double Fz, dFzdz;
        double Tbase = 1500, dT0dl, dT0dl_tmp, T0;
        double lat = lonlat_d[id*2 + 1];
        double lon = lonlat_d[id*2];
        double H = 580e3;
        int nlat;

        Fz = 0.5*(1-pow(tanh((z-z0)/dz0),3))*sin(M_PI*z/z1);
        dFzdz = 0.5*(-3/dz0*pow(tanh((z-z0)/dz0),2))*(1-pow(tanh((z-z0)/dz0),2))*sin(M_PI*z/z1)\
                + 0.5*(1-pow(tanh((z-z0)/dz0),3))*cos(M_PI*z/z1)*M_PI/z1;

        if (lat>=0) {
          U0 = U*pow(sin(M_PI*pow(sin(lat),2)),3)*Fz;
          dT0dl = dTemp_dphi(U,H,Rd,Omega,A,lat,Fz,dFzdz);
          T0 = Tbase;  //constant of integration
          for (nlat=1;nlat<=50;nlat++) { //do trapezoid rule. like a crazy person
            dT0dl_tmp = dTemp_dphi(U,H,Rd,Omega,A,nlat*lat/50,Fz,dFzdz); // ltmp = nlat*lat/50
            T0 += dT0dl_tmp + dT0dl;
            dT0dl = dT0dl_tmp;
          }
          T0 *= lat/50/2;
        }

//      Update temperature
        temperature_d[id*nv + lev] = T0;
        pressure_d[id*nv+lev] = pref*exp(-z/H);
        Rho_d[id*nv+lev] = pressure_d[id*nv+lev]/(Rd*temperature_d[id*nv+lev]);

//      Update momenta
        Mh_d[id*3*nv+lev*3+0] = U0*(-sin(lon))*Rho_d[id*nv+lev];
        Mh_d[id*3*nv+lev*3+1] = U0*(cos(lon))*Rho_d[id*nv+lev];
        Mh_d[id*3*nv+lev*3+2] = 0.0;
    }
}
