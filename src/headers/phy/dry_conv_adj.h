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
// Description: dry convective adjustment scheme
//
//
//
// Known limitations: None
//
// Known issues: None.
//
//
// If you use this code please cite the following reference:
//
//       [1] Mendonca, J.M., Grimm, S.L., Grosheintz, L., & Heng, K., ApJ, 829, 115, 2016
//
// Current Code Owners: Joao Mendonca (joao.mendonca@space.dtu.dk)
//                      Russell Deitrick (russell.deitrick@csh.unibe.ch)
//                      Urs Schroffenegger (urs.schroffenegger@csh.unibe.ch)
//
// History:
// Version Date       Comment
// ======= ====       =======
// 2.0     30/11/2018 Released version (RD & US)
// 1.0     16/08/2017 Released version  (JM)
//
////////////////////////////////////////////////////////////////////////

#include <math.h>


__global__ void dry_conv_adj(double *Pressure_d,    // Pressure [Pa]
                             double *Pressureh_d,   // Mid-point pressure [Pa]
                             double *Temperature_d, // Temperature [K]
                             double *profx_Qheat_d,
                             double *pt_d,        // Potential temperature [K]
                             double *Rho_d,       // Density [m^3/kg]
                             double *Cp_d,        // Specific heat capacity [J/kg/K]
                             double *Rd_d,        // Gas constant [J/kg/K]
                             double  Gravit,      // Gravity [m/s^2]
                             double *Altitude_d,  // Altitudes of the layers
                             double *Altitudeh_d, // Altitudes of the interfaces
                             double  time_step,
                             int  conv_adj_iter, // number of iterations of entire algorithm allowed
                             bool soft_adjust,
                             int  num, // Number of columns
                             int  nv) {          // Vertical levels
    //
    //  Description: Mixes entropy vertically on statically unstable columns
    //

    //

    int         id = blockIdx.x * blockDim.x + threadIdx.x;

    // Local arrays
    // double      theta[nv];
    // double      ph[nv + 1];

    // stability threshold
    double      stable = 0.0;

    double      ps, psm;
    double      pp, ptop;

    double      xi, xip, xim, a, b;

    if (id < num) {

                int  iter   = 0;
        bool repeat = true; //will repeat entire
        while ((repeat == true) && (iter < conv_adj_iter)) {
            // for (iter = 0; iter < ITERMAX; iter++) {
            // calculate pressure at the interfaces
            for (int lev = 0; lev <= nv; lev++) {
                if (lev == 0) {
                    // extrapolate to lower boundary
                    psm = Pressure_d[id * nv + 1]
                          - Rho_d[id * nv + 0] * Gravit * (-Altitude_d[0] - Altitude_d[1]);
                    ps                             = 0.5 * (Pressure_d[id * nv + 0] + psm);
                    Pressureh_d[id * (nv + 1) + 0] = ps;
                }
                else if (lev == nv) {
                    // extrapolate to top boundary
                    pp = Pressure_d[id * nv + nv - 2]
                         - Rho_d[id * nv + nv - 1] * Gravit
                               * (2 * Altitudeh_d[nv] - Altitude_d[nv - 1] - Altitude_d[nv - 2]);
                    if (pp < 0)
                        pp = 0; //prevents pressure from going negative
                    ptop                             = 0.5 * (Pressure_d[id * nv + nv - 1] + pp);
                    Pressureh_d[id * (nv + 1) + lev] = ptop;
                }
                else {
                    // interpolation between layers
                    xi  = Altitudeh_d[lev];
                    xim = Altitude_d[lev - 1];
                    xip = Altitude_d[lev];
                    a   = (xi - xip) / (xim - xip);
                    b   = (xi - xim) / (xip - xim);
                    Pressureh_d[id * (nv + 1) + lev] =
                        Pressure_d[id * nv + lev - 1] * a + Pressure_d[id * nv + lev] * b;
                }
            }

            // Compute Potential Temperature
            for (int lev = 0; lev < nv; lev++) {
                pt_d[id * nv + lev] =
                    Temperature_d[id * nv + lev]
                    * pow(Pressureh_d[id * (nv + 1) + 0] / Pressure_d[id * nv + lev],
                          Rd_d[id * nv + lev] / Cp_d[id * nv + lev]);
            }

            bool done_col = false;
            while (done_col == false) { // Unstable  column?
                int top = 0;
                int bot = nv - 1;

                for (int lev = 0; lev < nv - 1; lev++) {
                    // sweep upward, find lowest unstable layer
                    if (pt_d[id * nv + lev + 1] - pt_d[id * nv + lev] < stable) {
                        if (bot > lev)
                            bot = lev;
                    }
                }

                for (int lev = bot; lev < nv - 1; lev++) {
                    // sweep upward from unstable layer, find top
                    if (pt_d[id * nv + lev + 1] - pt_d[id * nv + lev] > stable) {
                        top = lev;
                        break;
                    }
                    else {
                        top = nv - 1;
                    }
                }

                if (bot < nv - 1) {
                    int    extend = 1;
                    double thnew;

                    while (extend == 1) {
                        double h   = 0.0; //Enthalpy;
                        double sum = 0.0;
                        extend     = 0;

                        for (int lev = bot; lev <= top; lev++) {
                            // calc adiabatic pressure, integrate upward for new pot. temp.
                            double pu = Pressureh_d[id * (nv + 1) + lev + 1];
                            double pl = Pressureh_d[id * (nv + 1) + lev];
                            double pi =
                                pow(Pressure_d[id * nv + lev] / Pressureh_d[id * (nv + 1) + 0],
                                    Rd_d[id * nv + lev]
                                        / Cp_d[id * nv
                                               + lev]); // adiabatic pressure wrt bottom of column
                            double deltap = pl - pu;

                            h   = h + pt_d[id * nv + lev] * pi * deltap;
                            sum = sum + pi * deltap;
                        }
                        thnew = h / sum;

                        // if (bot <= 0 && top >= nv - 1) {
                        //     // no need to extend again
                        //     extend = 0;
                        // }

                        if (bot > 0) {
                            // repeat if new pot. temp. is less than lower boundary p.t.
                            if ((thnew - pt_d[id * nv + bot - 1]) < stable) {
                                bot    = bot - 1;
                                extend = 1;
                            }
                        }

                        if (top < nv - 1) {
                            // repeat if new pot. temp. is greater p.t. above
                            if ((pt_d[id * nv + top + 1] - thnew) < stable) {
                                top    = top + 1;
                                extend = 1;
                            }
                        }
                    }

                    for (int lev = bot; lev <= top; lev++) {
                        pt_d[id * nv + lev] = thnew; // set new potential temperature
                    }
                }
                else {
                    done_col = true; //no unstable layers
                }
            }

            repeat = false;
            iter += 1;

            if (soft_adjust) {
                double Ttmp, Ptmp;

                for (int lev = 0; lev < nv; lev++) {
                    Ttmp = pt_d[id * nv + lev]
                           * pow(Pressure_d[id * nv + lev] / Pressureh_d[id * (nv + 1) + 0],
                                 Rd_d[id * nv + lev] / Cp_d[id * nv + lev]);
                    Ptmp = Ttmp * Rd_d[id * nv + lev] * Rho_d[id * nv + lev];
                    //reset pt value to beginning of time step
                    pt_d[id * nv + lev] =
                        Temperature_d[id * nv + lev]
                        * pow(Pressure_d[id * nv + lev] / Pressureh_d[id * (nv + 1) + 0],
                              -Rd_d[id * nv + lev] / Cp_d[id * nv + lev]);

                    profx_Qheat_d[id * nv + lev] +=
                        (Cp_d[id * nv + lev] - Rd_d[id * nv + lev]) / Rd_d[id * nv + lev]
                        * (Ptmp - Pressure_d[id * nv + lev]) / time_step;
                    //does not repeat
                }
            }
            // Compute Temperature & pressure from potential temperature
            else {
                for (int lev = 0; lev < nv; lev++) {
                    Temperature_d[id * nv + lev] =
                        pt_d[id * nv + lev]
                        * pow(Pressure_d[id * nv + lev] / Pressureh_d[id * (nv + 1) + 0],
                              Rd_d[id * nv + lev] / Cp_d[id * nv + lev]);
                    Pressure_d[id * nv + lev] =
                        Temperature_d[id * nv + lev] * Rd_d[id * nv + lev] * Rho_d[id * nv + lev];
                    //check pt again
                    pt_d[id * nv + lev] =
                        Temperature_d[id * nv + lev]
                        * pow(Pressure_d[id * nv + lev] / Pressureh_d[id * (nv + 1) + 0],
                              -Rd_d[id * nv + lev] / Cp_d[id * nv + lev]);
                    if (lev > 0) {
                        if (pt_d[id * nv + lev] - pt_d[id * nv + lev - 1] < stable)
                            repeat = true;
                    }
                }
            }
        }
        //printf("id = %d, iter = %d\n", id, iter);
        
        
        

    }
}


__device__ void bezier_altitude_interpolation_dry_conv(int id, int nlay, int iter, double* xi, double* yi, double x, double &y) {

    double dx, dx1, dy, dy1;
    double w, yc, t;
    //xc = (xi(1) + xi(2))/2.0_dp ! Control point (no needed here, implicitly included)
    dx  =   xi[iter]                    -       xi[iter + 1];
    dx1 =   xi[iter - 1]                -       xi[iter];
    dy  =   yi[id * nlay + iter]        -       yi[id * nlay + iter + 1];
    dy1 =   yi[id * nlay + iter - 1]    -       yi[id * nlay + iter];

    if (x > xi[iter + 1] && x < xi[iter])
    {
        // left hand side interpolation
        w   =   dx1 / (dx + dx1);

        yc  =   yi[id * nlay + iter] -
                dx / 2.0 * 
                (
                    w * dy / dx + 
                    (1.0 - w) * dy1 / dx1
                );

        t   =   (x - xi[iter + 1]) / dx;

        y   =   pow(1.0 - t, 2) * yi[id * nlay + iter + 1] + 
                2.0 * t * (1.0 - t) * yc + 
                pow(t, 2) * yi[id * nlay + iter];
    } else
    {
        // right hand side interpolation
        w   =   dx / (dx + dx1);

        yc  =   yi[id * nlay + iter] + 
                dx1 / 2.0 * 
                (
                    w * dy1 / dx1 +
                    (1.0 - w) * dy / dx
                );

        t   =   (x - xi[iter]) / (dx1);

        y   =   pow(1.0 - t, 2) * yi[id * nlay + iter] + 
                2.0 * t * (1.0 - t) * yc + 
                pow(t, 2) * yi[id * nlay + iter - 1];

    }
}


__global__ void ray_dry_conv_adj(double *Pressure_d,    // Pressure [Pa]
                             double *Pressureh_d,   // Mid-point pressure [Pa]
                             double *dT_conv_d,
                             double *Temperature_d, // Temperature [K]
                             double *profx_Qheat_d,
                             double *pt_d,          // Potential temperature [K]
                             double *Rho_d,         // Density [m^3/kg]
                             double *Cp_d,          // Specific heat capacity [J/kg/K]
                             double *Rd_d,          // Gas constant [J/kg/K]
                             double  Gravit,        // Gravity [m/s^2]
                             double *Altitude_d,    // Altitudes of the layers
                             double *Altitudeh_d,   // Altitudes of the interfaces
                             double timestep,       // time step [s]
                             int conv_adj_iter, // number of iterations of entire algorithm allowed
                             bool soft_adjust,
                             int num,           // Number of columns
                             int nv) {          // Vertical levels
    //
    //  Description: Mixes entropy vertically on statically unstable columns
    //

    //

    int         id = blockIdx.x * blockDim.x + threadIdx.x;

    // Local arrays
    // double      theta[nv];
    // double      ph[nv + 1];

    // stability threshold
    double      stable = 0.0;

    double      ps, psm;
    double      pp, ptop;

    double      xi, xip, xim, a, b;

    if (id < num) {
        int  iter   = 0;
        bool repeat = true; //will repeat entire

        bool ray_mode = false;

        if (ray_mode == true)
        {

            // constants & parameters

            int itermax1 = 4;
            int itermax2 = 5;
            const double small = 1e-6;
            double const euler = 2.71828182845904523536028;
            double dT_factor =  timestep;

            // work variables
            int i;
            bool did_adj;
            double pfact, Tbar;
            double d_p;
            double d_p_lower;
            double d_p_pfact;
            double condi;
            double d_T;
            double d_T_lower;

            for (i = 0; i <nv; i++){
                    dT_conv_d[id * nv + i] =  Temperature_d[id * nv + i];
            }

            for (int lev = 0; lev <= nv; lev++) {
                        if (lev == 0) {
                            // extrapolate to lower boundary
                            psm = Pressure_d[id * nv + 1]
                                - Rho_d[id * nv + 0] * Gravit * (-Altitude_d[0] - Altitude_d[1]);
                            ps                             = 0.5 * (Pressure_d[id * nv + 0] + psm);
                            Pressureh_d[id * (nv + 1) + 0] = ps;
                        }
                        else if (lev == nv) {
                            // extrapolate to top boundary
                            pp = Pressure_d[id * nv + nv - 2]
                                - Rho_d[id * nv + nv - 1] * Gravit
                                    * (2 * Altitudeh_d[nv] - Altitude_d[nv - 1] - Altitude_d[nv - 2]);
                            if (pp < 0)
                                pp = 0; //prevents pressure from going negative
                            ptop                             = 0.5 * (Pressure_d[id * nv + nv - 1] + pp);
                            Pressureh_d[id * (nv + 1) + lev] = ptop;
                        }
                        else {
                            // interpolation between layers
                            /*
                            xi  = Altitudeh_d[lev];
                            xim = Altitude_d[lev - 1];
                            xip = Altitude_d[lev];
                            a   = (xi - xip) / (xim - xip);
                            b   = (xi - xim) / (xip - xim);
                            Pressureh_d[id * (nv + 1) + lev] =
                                Pressure_d[id * nv + lev - 1] * a + Pressure_d[id * nv + lev] * b;
                            */

                            bezier_altitude_interpolation_dry_conv( id,
                                                                    nv, 
                                                                    lev, 
                                                                    Altitude_d, 
                                                                    Pressure_d, 
                                                                    Altitudeh_d[lev], 
                                                                    Pressureh_d[id * (nv + 1) + lev]);
                        }

                    }
            
            // start operations
            while ((repeat == true) && (iter < itermax1))
            {
                
                for (int iter2 = 0; iter2 < itermax2; iter2++)
                {
                    did_adj = false;

                    // Downward pass
                    for (i = nv - 1; i > 0; i--){
                        
                        /*
                        d_p = Rho_d[id * nv + i] * Gravit * (Altitudeh_d[i+1] - Altitudeh_d[i]);
                        d_p_lower = Rho_d[id * nv + i-1] * Gravit * (Altitudeh_d[i] - Altitudeh_d[i-1]);
                        d_p_pfact = Rho_d[id * nv + i-1] * Gravit * (Altitude_d[i] - Altitude_d[i-1]);

                        pfact = pow( ( (Pressure_d[id * nv + i - 1] - d_p_pfact) / Pressure_d[id * nv + i - 1]) ,
                                Rd_d[id * nv + i-1] / Cp_d[id * nv + i-1]);
                        condi = (dT_conv_d[id * nv + i - 1]  * pfact - small);
                        */

                        d_p = Pressureh_d[id * (nv + 1) + i] - Pressureh_d[id * (nv + 1) + i + 1];
                        d_p_lower = Pressureh_d[id * (nv + 1) + i - 1] - Pressureh_d[id * (nv + 1) + i];

                        pfact = pow( ( Pressure_d[id * nv + i] / Pressure_d[id * nv + i - 1]) ,
                                Rd_d[id * nv + i-1] / Cp_d[id * nv + i-1]);
                        condi = (dT_conv_d[id * nv + i - 1]  * pfact - small);

                        if (dT_conv_d[id * nv + i] < condi) {

                            Tbar = (d_p * dT_conv_d[id * nv + i] +
                                d_p_lower * dT_conv_d[id * nv + i - 1]) /
                                (d_p_lower + d_p); 

                            dT_conv_d[id * nv + i - 1] =  (d_p_lower + d_p) *Tbar / (d_p_lower + pfact * d_p); 

                            dT_conv_d[id * nv + i] = dT_conv_d[id * nv + i - 1] * pfact;

                            //Temperature_d[id * nv + i] =  Temperature_d[id * nv + i] + ((d_T - Temperature_d[id * nv + i]) / timestep);

                            did_adj = true;
                        }
                    }

                    // Upward pass
                    for (i = 1; i < nv; i++) { 
                        /*
                        d_p = Rho_d[id * nv + i] * Gravit * (Altitudeh_d[i+1] - Altitudeh_d[i]);
                        d_p_lower = Rho_d[id * nv + i-1] * Gravit * (Altitudeh_d[i] - Altitudeh_d[i-1]);
                        d_p_pfact = Rho_d[id * nv + i-1] * Gravit * (Altitude_d[i] - Altitude_d[i-1]);

                        pfact = pow( ( (Pressure_d[id * nv + i - 1] - d_p_pfact) / Pressure_d[id * nv + i - 1]) ,
                                Rd_d[id * nv + i-1] / Cp_d[id * nv + i-1]);
                        condi = (dT_conv_d[id * nv + i - 1] * pfact - small);
                        */

                        d_p = Pressureh_d[id * (nv + 1) + i] - Pressureh_d[id * (nv + 1) + i + 1];
                        d_p_lower = Pressureh_d[id * (nv + 1) + i - 1] - Pressureh_d[id * (nv + 1) + i];

                        pfact = pow( ( Pressure_d[id * nv + i] / Pressure_d[id * nv + i - 1]) ,
                                Rd_d[id * nv + i-1] / Cp_d[id * nv + i-1]);
                        condi = (dT_conv_d[id * nv + i - 1]  * pfact - small);

                        if (dT_conv_d[id * nv + i] < condi) {

                            
                            Tbar = (d_p * dT_conv_d[id * nv + i] +
                                d_p_lower * dT_conv_d[id * nv + i - 1]) /
                                (d_p_lower + d_p); 

                            dT_conv_d[id * nv + i - 1] =  (d_p_lower + d_p) *Tbar / (d_p_lower + pfact * d_p); 

                            dT_conv_d[id * nv + i] = dT_conv_d[id * nv + i - 1] * pfact;

                            //Temperature_d[id * nv + i] =  Temperature_d[id * nv + i] + ((d_T - Temperature_d[id * nv + i]) / timestep);

                            did_adj = true;
                        }
                    }

                    // ! If no adjustment required, exit the loop
                    if (did_adj == false)
                    {
                        break;
                    }
                }

                repeat = false;
                iter += 1;
                if (soft_adjust) {
                    double Ttmp, Ptmp;

                    for (int lev = 0; lev < nv; lev++) {
                        //reset pt value to beginning of time step
                        
                        pt_d[id * nv + lev] =
                            Temperature_d[id * nv + lev]
                            * pow(Pressureh_d[id * (nv + 1) + 0] / Pressure_d[id * nv + lev],
                                  Rd_d[id * nv + lev] / Cp_d[id * nv + lev]);


                        profx_Qheat_d[id * nv + lev] += (dT_conv_d[id * nv + i] - Temperature_d[id * nv + i])  /
                            (timestep);

                        //does not repeat
                    }
                } else
                {
                    for (i = 0; i <nv; i++){ 
                        Temperature_d[id * nv + i] = Temperature_d[id * nv + i] +
                            dT_factor * (dT_conv_d[id * nv + i] - Temperature_d[id * nv + i])  /
                            (itermax1 * timestep);
                        
                    }

                    for (i = 1; i <nv; i++){
                        Pressure_d[id * nv + i] = Rho_d[id * nv + i] * Rd_d[id * nv + i] * Temperature_d[id * nv + i];
                        // converve the mass might be better than hypsometric equation

                        
                        // interpolation between layers
                        /*
                        xim = Pressure_d[i - 1];
                        xip = Pressure_d[i];
                        a   = (xip) / (xim + xip);
                        b   = (xim) / (xip + xim);
                        

                        xi  = Altitudeh_d[i];
                        xim = Altitude_d[i - 1];
                        xip = Altitude_d[i];
                        a   = (xi - xip) / (xim - xip);
                        b   = (xi - xim) / (xip - xim);
                        
                        Pressure_d[id * nv + i] = Pressure_d[id * nv + i] + 
                            (   Pressure_d[id * nv + i - 1] *
                                    pow(euler,
                                        Gravit * (Altitude_d[i] - Altitude_d[i-1]) /
                                        (
                                            0.5*(Temperature_d[id * nv + i] + Temperature_d[id * nv + i - 1]) *
                                            0.5*(Rd_d[id * nv + i] + Rd_d[id * nv + i -1])
                                            //(a * Temperature_d[id * nv + i] + b * Temperature_d[id * nv + i - 1]) *
                                            //(a * Rd_d[id * nv + i] + b * Rd_d[id * nv + i -1])
                                        )
                                    ) -
                                Pressure_d[id * nv + i]
                            ) /
                            (timestep * itermax1) ;

                        */
                        
                        
                        
                    }
                    for (int lev = 0; lev <= nv; lev++) {
                        if (lev == 0) {
                            // extrapolate to lower boundary
                            psm = Pressure_d[id * nv + 1]
                                - Rho_d[id * nv + 0] * Gravit * (-Altitude_d[0] - Altitude_d[1]);
                            ps                             = 0.5 * (Pressure_d[id * nv + 0] + psm);
                            Pressureh_d[id * (nv + 1) + 0] = ps;
                        }
                        else if (lev == nv) {
                            // extrapolate to top boundary
                            pp = Pressure_d[id * nv + nv - 2]
                                - Rho_d[id * nv + nv - 1] * Gravit
                                    * (2 * Altitudeh_d[nv] - Altitude_d[nv - 1] - Altitude_d[nv - 2]);
                            if (pp < 0)
                                pp = 0; //prevents pressure from going negative
                            ptop                             = 0.5 * (Pressure_d[id * nv + nv - 1] + pp);
                            Pressureh_d[id * (nv + 1) + lev] = ptop;
                        }
                        else {
                            // interpolation between layers
                            /*
                            xi  = Altitudeh_d[lev];
                            xim = Altitude_d[lev - 1];
                            xip = Altitude_d[lev];
                            a   = (xi - xip) / (xim - xip);
                            b   = (xi - xim) / (xip - xim);
                            Pressureh_d[id * (nv + 1) + lev] =
                                Pressure_d[id * nv + lev - 1] * a + Pressure_d[id * nv + lev] * b;
                            */
                            
                            bezier_altitude_interpolation_dry_conv( id,
                                                                    nv, 
                                                                    lev, 
                                                                    Altitude_d, 
                                                                    Pressure_d, 
                                                                    Altitudeh_d[lev], 
                                                                    Pressureh_d[id * (nv + 1) + lev]);
                        }

                        // Compute Potential Temperature
                        pt_d[id * nv + lev] =
                            Temperature_d[id * nv + lev]
                            * pow(Pressureh_d[id * (nv + 1) + 0] / Pressure_d[id * nv + lev],
                                Rd_d[id * nv + lev] / Cp_d[id * nv + lev]);

                    }

                    repeat = true;

                    for (i = 0; i <nv; i++){
                            dT_conv_d[id * nv + i] =  Temperature_d[id * nv + i];
                    }
                    
                }
                
            }
  
        } else
        {
             while ((repeat == true) && (iter < conv_adj_iter)) {
                // for (iter = 0; iter < ITERMAX; iter++) {
                // calculate pressure at the interfaces
                for (int lev = 0; lev <= nv; lev++) {
                    if (lev == 0) {
                        // extrapolate to lower boundary
                        psm = Pressure_d[id * nv + 1]
                            - Rho_d[id * nv + 0] * Gravit * (-Altitude_d[0] - Altitude_d[1]);
                        ps                             = 0.5 * (Pressure_d[id * nv + 0] + psm);
                        Pressureh_d[id * (nv + 1) + 0] = ps;
                    }
                    else if (lev == nv) {
                        // extrapolate to top boundary
                        pp = Pressure_d[id * nv + nv - 2]
                            - Rho_d[id * nv + nv - 1] * Gravit
                                * (2 * Altitudeh_d[nv] - Altitude_d[nv - 1] - Altitude_d[nv - 2]);
                        if (pp < 0)
                            pp = 0; //prevents pressure from going negative
                        ptop                             = 0.5 * (Pressure_d[id * nv + nv - 1] + pp);
                        Pressureh_d[id * (nv + 1) + lev] = ptop;
                    }
                    else {
                        // interpolation between layers
                        xi  = Altitudeh_d[lev];
                        xim = Altitude_d[lev - 1];
                        xip = Altitude_d[lev];
                        a   = (xi - xip) / (xim - xip);
                        b   = (xi - xim) / (xip - xim);
                        Pressureh_d[id * (nv + 1) + lev] =
                            Pressure_d[id * nv + lev - 1] * a + Pressure_d[id * nv + lev] * b;
                    }
                }

                // Compute Potential Temperature
                for (int lev = 0; lev < nv; lev++) {
                    pt_d[id * nv + lev] =
                        Temperature_d[id * nv + lev]
                        * pow(Pressureh_d[id * (nv + 1) + 0] / Pressure_d[id * nv + lev],
                            Rd_d[id * nv + lev] / Cp_d[id * nv + lev]);
                }

                bool done_col = false;
                while (done_col == false) { // Unstable  column?
                    int top = 0;
                    int bot = nv - 1;

                    for (int lev = 0; lev < nv - 1; lev++) {
                        // sweep upward, find lowest unstable layer
                        if (pt_d[id * nv + lev + 1] - pt_d[id * nv + lev] < stable) {
                            if (bot > lev)
                                bot = lev;
                        }
                    }

                    for (int lev = bot; lev < nv - 1; lev++) {
                        // sweep upward from unstable layer, find top
                        if (pt_d[id * nv + lev + 1] - pt_d[id * nv + lev] > stable) {
                            top = lev;
                            break;
                        }
                        else {
                            top = nv - 1;
                        }
                    }

                    if (bot < nv - 1) {
                        int    extend = 1;
                        double thnew;

                        while (extend == 1) {
                            double h   = 0.0; //Enthalpy;
                            double sum = 0.0;
                            extend     = 0;

                            for (int lev = bot; lev <= top; lev++) {
                                // calc adiabatic pressure, integrate upward for new pot. temp.
                                double pu = Pressureh_d[id * (nv + 1) + lev + 1];
                                double pl = Pressureh_d[id * (nv + 1) + lev];
                                double pi =
                                    pow(Pressure_d[id * nv + lev] / Pressureh_d[id * (nv + 1) + 0],
                                        Rd_d[id * nv + lev]
                                            / Cp_d[id * nv
                                                + lev]); // adiabatic pressure wrt bottom of column
                                double deltap = pl - pu;

                                h   = h + pt_d[id * nv + lev] * pi * deltap;
                                sum = sum + pi * deltap;
                            }
                            thnew = h / sum;

                            // if (bot <= 0 && top >= nv - 1) {
                            //     // no need to extend again
                            //     extend = 0;
                            // }

                            if (bot > 0) {
                                // repeat if new pot. temp. is less than lower boundary p.t.
                                if ((thnew - pt_d[id * nv + bot - 1]) < stable) {
                                    bot    = bot - 1;
                                    extend = 1;
                                }
                            }

                            if (top < nv - 1) {
                                // repeat if new pot. temp. is greater p.t. above
                                if ((pt_d[id * nv + top + 1] - thnew) < stable) {
                                    top    = top + 1;
                                    extend = 1;
                                }
                            }
                        }

                        for (int lev = bot; lev <= top; lev++) {
                            pt_d[id * nv + lev] = thnew; // set new potential temperature
                        }
                    }
                    else {
                        done_col = true; //no unstable layers
                    }
                }

                repeat = false;
                iter += 1;

                if (soft_adjust) {
                    double Ttmp, Ptmp;

                    for (int lev = 0; lev < nv; lev++) {
                        Ttmp = pt_d[id * nv + lev]
                            * pow(Pressure_d[id * nv + lev] / Pressureh_d[id * (nv + 1) + 0],
                                    Rd_d[id * nv + lev] / Cp_d[id * nv + lev]);
                        Ptmp = Ttmp * Rd_d[id * nv + lev] * Rho_d[id * nv + lev];
                        //reset pt value to beginning of time step
                        pt_d[id * nv + lev] =
                            Temperature_d[id * nv + lev]
                            * pow(Pressure_d[id * nv + lev] / Pressureh_d[id * (nv + 1) + 0],
                                -Rd_d[id * nv + lev] / Cp_d[id * nv + lev]);

                        profx_Qheat_d[id * nv + lev] +=
                            (Cp_d[id * nv + lev] - Rd_d[id * nv + lev]) / Rd_d[id * nv + lev]
                            * (Ptmp - Pressure_d[id * nv + lev]) / timestep;
                        //does not repeat
                    }
                }
                // Compute Temperature & pressure from potential temperature
                else {
                    for (int lev = 0; lev < nv; lev++) {
                        Temperature_d[id * nv + lev] =
                            pt_d[id * nv + lev]
                            * pow(Pressure_d[id * nv + lev] / Pressureh_d[id * (nv + 1) + 0],
                                Rd_d[id * nv + lev] / Cp_d[id * nv + lev]);
                        Pressure_d[id * nv + lev] =
                            Temperature_d[id * nv + lev] * Rd_d[id * nv + lev] * Rho_d[id * nv + lev];
                        //check pt again
                        pt_d[id * nv + lev] =
                            Temperature_d[id * nv + lev]
                            * pow(Pressure_d[id * nv + lev] / Pressureh_d[id * (nv + 1) + 0],
                                -Rd_d[id * nv + lev] / Cp_d[id * nv + lev]);
                        if (lev > 0) {
                            if (pt_d[id * nv + lev] - pt_d[id * nv + lev - 1] < stable)
                                repeat = true;
                        }
                    }
                }
            }
            //printf("id = %d, iter = %d\n", id, iter);
        }
        
        

        
        //printf("id = %d, iter = %d\n", id, iter);
    }
}
