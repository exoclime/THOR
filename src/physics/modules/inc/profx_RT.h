//
//
//
//
// Description:
//
//
// Method:
//   Calls to
//
// Known limitations:
//
//
// Known issues:
//
//
// Current Code Owner:
//
// History:
// Version Date       Comment
// ======= ====       =======
//
// 1.0
//
////////////////////////////////////////////////////////////////////////


__host__ double solve_kepler(double mean_anomaly,
                             double ecc         ) {

// Solve Kepler's equation (see Murray & Dermott 1999)
// Get eccentric anomaly from mean anomaly and eccentricity

  double ecc_anomaly, fi, fi_1, fi_2, fi_3, di_1, di_2, di_3;

  ecc_anomaly = mean_anomaly + sign(sin(mean_anomaly))*0.85*e;
  di_3 = 1.0;

  while (di_3 > 1e-15) {
    fi = ecc_anomaly - ecc*sin(ecc_anomaly) - mean_anomaly;
    fi_1 = 1.0 - ecc*cos(ecc_anomaly);
    fi_2 = ecc*sin(ecc_anomaly);
    fi_3 = ecc*cos(ecc_anomaly);
    di_1 = -fi / fi_1;
    di_2 = -fi / (fi_1 + 0.5*di_1*fi_2);
    di_3 = -fi / (fi_1 + 0.5*di_2*fi_2 + 1./6.*di_2*di_2*fi_3);
    ecc_anomaly = ecc_anomaly + di_3;
  }
  return ecc_anomaly;
}

__host__ double calc_r_orb(double ecc_anomaly,
                           double ecc        ) {

// Calc relative distance between planet and star (units of semi-major axis)

  double r = 1.0 - ecc*cos(ecc_anomaly);
  return r;
}

__host__ double ecc2true_anomaly(double ecc_anomaly,
                                 double ecc        ){

// Convert from eccentric to true anomaly

  double tanf2, true_anomaly;
  tanf2 = sqrt((1.0+ecc)/(1.0-ecc))*tan(ecc_anomaly/2.0);
  true_anomaly = 2.0*atan(tanf2);
  if (true_anomaly<0.0) true_anomaly += 2*M_PI;
  return true_anomaly;
}

__host__ double true2ecc_anomaly(double true_anomaly,
                                 double ecc         ){

// Convert from true to eccentric anomaly

  double cosE, ecc_anomaly;
  while (true_anomaly<0.0) true_anomaly += 2*M_PI;
  while (true_anomaly>=2*M_PI) true_anomaly -= 2*M_PI;
  cosE = (cos(true_anomaly) + ecc) / (1.0 + ecc*cos(true_anomaly));
  if (true_anomaly < M_PI){
    ecc_anomaly = acos(cosE);
  } else {
    ecc_anomaly = 2*M_PI - acos(cosE);
  }
  return ecc_anomaly;
}

// need to calc mean_anomaly_i at start
// ecc, obliquity, mean_motion, longp, true_long_i and alpha_i are user inputs
__host__ void update_spin_orbit() {

// Update the insolation related parameters for spin and orbit

  mean_anomaly = fmod((mean_motion*time + mean_anomaly_i),(2*M_PI));

  ecc_anomaly = fmod(solve_kepler(mean_anomaly, ecc),(2*M_PI));

  r_orb = calc_r_orb(ecc_anomaly, ecc);
  true_long = fmod((ecc2true_anomaly(ecc_anomaly, ecc) + longp),(2*M_PI));

  sin_decl = sin(obliquity)*sin(true_long);
  cos_decl = sqrt(1.0 - sin_decl*sin_decl);
  alpha = -Omega*time + true_long - true_long_i + alpha_i;
}

__device__ double calc_insol(double *lonlat_d,  //latitude/longitude grid
                             double alpha    ,  //current RA of star (relative to zero long on planet)
                             double sin_decl ,  //declination of star
                             double cos_decl ,
                             double r_orb    ){ //host star planet distance (units of semi-major axis)

// Calculate the insolation (scaling) at a point on the surface

  double insol_q;

  insol_q = pow(r,-2)*( sin(lonlat_d[id*nv+1])*sin_decl + \
                      cos(lonlat_d[id*nv+1])*cos_decl*cos(lonlat_d[id*nv]-alpha) );
  return insol_q;     //insolation scaling factor
}


__device__ void radcsw(double *phtemp     ,
                       double  coszrs     ,
                       double *dtemp      ,
                       double *tau_d      ,
                       double *fnet_up_d  ,
                       double *fnet_dn_d  ,
                       double  incflx     ,
                       double  alb        ,
                       double  tausw      ,
                       double  ps0        ,
                       double  Cp         ,
                       double  gravit     ,
                       int     id         ,
                       int     nv         ){

//  Calculate upward, downward, and net flux.
//  Downward Directed Radiation

    double gocp;
    double tau = (tausw/ps0)*(phtemp[id*(nv+1)+nv]);
    double flux_top = incflx*coszrs*(1.0-alb);

	// Extra layer to avoid over heating at the top.
    fnet_dn_d[id*(nv+1) + nv] = flux_top*exp(-(1.0/coszrs)*tau);
    // if (isnan(fnet_dn_d[id*(nv+1)+nv]) || isinf(fnet_dn_d[id*(nv+1)+nv])) {
    //   printf("%d, %d\n",id,nv);
    // }

    // Normal integration
	for(int lev = nv;lev >= 1;lev--) fnet_dn_d[id*(nv+1) + lev-1] = fnet_dn_d[id*(nv+1) + lev]*exp(-(1.0/coszrs)*tau_d[id*nv*2 + (lev-1)*2]);
	for(int lev = 0;lev <= nv;lev++) fnet_up_d[id*(nv+1) + lev]   = 0.0;

	// Update temperature rates.
  // if (id==0) printf("shortwave down = \n");
  for(int lev = 0;lev < nv;lev++){
        gocp  = gravit/Cp;
        dtemp[id*nv+lev]     = gocp * ((fnet_up_d[id*(nv+1) + lev] - fnet_dn_d[id*(nv+1) + lev]) - (fnet_up_d[id*(nv+1) + lev+1] - fnet_dn_d[id*(nv+1) + lev+1])) / (phtemp[id*(nv+1)+lev] - phtemp[id*(nv+1)+lev +1]);

        // if (isnan(dtemp[id*nv+lev])) {
        //   printf("%d, %d\n",id,lev);
        // }
        // if (id == 0){
        //   printf("%f, ",fnet_dn_d[id*(nv+1)+lev-1]);
        // }
    }
}

__device__ int factorial_num(int number) {

    int temp;
    int n1;
    temp = 1;
    n1   = number;
    while(n1 != 0){
        temp = temp*n1;
        n1 = n1 - 1;
    }
    return temp;
}

__device__ double source_func_lin(double bb      ,
                                  double bl      ,
                                  double bt      ,
                                  double tau     ,
                                  double diff_fac){

    double e1 = 0.0;
    double e2 = 0.0;
    double e  = 0.0;

    if(tau >= 1e-10){
        e1 = bb - bl + (bl-(diff_fac/(tau/2.0))*(bb-bl))*(1.0-exp(-(tau/2.0)/diff_fac));
        e2 = bl - bt + (bt-(diff_fac/(tau/2.0))*(bl-bt))*(1.0-exp(-(tau/2.0)/diff_fac));
        e  = e1 + e2 *  exp(-(tau/2.0)/diff_fac);
    }
    else{
        for(int i = 0; i < 5; i++){
            int fac = factorial_num(i+2);
            e1 = e1 + (pow(-1.0,i+2.0))*((bl + (i+1)*bl)/fac)*pow((tau/2.0)/diff_fac,i+1.0);
            e2 = e2 + (pow(-1.0,i+2.0))*((bt + (i+1)*bb)/fac)*pow((tau/2.0)/diff_fac,i+1.0);
        }
        e = e1 + e2*exp(-(tau/2.0)/diff_fac);
    }

    return e;
}

__device__ void radclw(double *phtemp         ,
                       double *ttemp          ,
                       double *thtemp         ,
                       double *dtemp      ,
                       double *tau_d      ,
                       double *fnet_up_d  ,
                       double *fnet_dn_d  ,
					   double diff_fac    ,
					   double tlow          ,
                       double Cp          ,
                       double gravit      ,
                       int    id          ,
                       int    nv          ){

    double gocp = gravit/Cp;
    double tb, tl, tt;
    double bb, bl, bt;

    double bc = 5.677036E-8;  //Stefan–Boltzmann constant W⋅m−2⋅K−4

//
//  Calculate upward, downward, and net flux.
//  Downward Directed Radiation
//
   fnet_dn_d[id*(nv+1) + nv] = 0.0; // Upper boundary
    for(int lev = nv-1;lev >= 0; lev--){
        double ed = 0.0;
        if(tau_d[id*nv*2 + 2*lev + 1] < 0.0) tau_d[id*nv*2 + 2*lev + 1] = 0.0;

        tb = thtemp[id*(nv+1)+lev];
        tl = ttemp[id*nv+lev];
        tt = thtemp[id*(nv+1)+lev+1];

        bb = bc*tb*tb*tb*tb;
        bl = bc*tl*tl*tl*tl;
        bt = bc*tt*tt*tt*tt;

        ed = source_func_lin(bb, bl, bt, tau_d[id*nv*2 + 2*lev + 1], diff_fac);

        fnet_dn_d[id*(nv+1) + lev] = ed + fnet_dn_d[id*(nv+1) + lev+1] * exp(-(1./diff_fac)*tau_d[id*nv*2 + 2*lev + 1]);
    }
//
//  Upward Directed Radiation
//
    fnet_up_d[id*(nv+1) + 0] = bc*tlow*tlow*tlow*tlow; // Lower boundary;
    if (fnet_up_d[id*(nv+1)+0] < fnet_dn_d[id*(nv+1)+0]) fnet_up_d[id*(nv+1)+0] = fnet_dn_d[id*(nv+1)+0];
    for(int lev = 1;lev <= nv; lev++){

        double eu = 0.0;

        if(tau_d[id*nv*2 + 2*(lev-1) + 1] < 0.0) tau_d[id*nv*2 + 2*(lev-1) + 1] = 0.0;

        tb = thtemp[id*(nv+1)+lev-1];
        tl = ttemp[id*nv+lev-1];
        tt = thtemp[id*(nv+1)+lev];

        bb = bc*tb*tb*tb*tb;
        bl = bc*tl*tl*tl*tl;
        bt = bc*tt*tt*tt*tt;

        eu = source_func_lin(bt, bl, bb, tau_d[id*nv*2 + 2*(lev-1) + 1], diff_fac);

        fnet_up_d[id*(nv+1) + lev] = eu + fnet_up_d[id*(nv+1) + lev-1] * exp(-(1./diff_fac)*tau_d[id*nv*2 + 2*(lev-1) + 1]);
    }

    // if(id==0) printf("\nlongwave down = \n");
    for(int lev = 0;lev < nv;lev++) {
        dtemp[id*nv+lev] = dtemp[id*nv+lev] + gocp * ((fnet_up_d[id*(nv+1) + lev] - fnet_dn_d[id*(nv+1) + lev]) - (fnet_up_d[id*(nv+1) + lev+1] - fnet_dn_d[id*(nv+1) + lev+1])) /
                                    (phtemp[id*(nv+1)+lev] - phtemp[id*(nv+1)+lev + 1]);

        // if (isnan(dtemp[id*nv+lev])) {
        //   printf("%d, %d\n",id,lev);
        // }
    }
    // if(id==0) printf("\nlongwave up = \n");
    // for(int lev = 0;lev < nv;lev++) {
    //     // if(id==0) printf("%f, ",fnet_up_d[id*(nv+1)+lev]);
    // }

}


__device__ void computetau (double *tau_d,
                            double *phtemp   ,
                            double  cosz ,
                            double  tausw,
                            double  taulw,
                            double  ps0  ,
                            int id       ,
                            int nv       ){

    for (int lev = 0; lev < nv; lev++ ){
        tau_d[id*2*nv + lev*2]     =(tausw/ps0)*(phtemp[id*(nv+1)+lev]-phtemp[id*(nv+1)+lev+1]);
        tau_d[id*2*nv + lev*2 + 1] =(taulw/ps0)*(phtemp[id*(nv+1)+lev]-phtemp[id*(nv+1)+lev+1]) + (taulw/pow(ps0, 2.0))*(pow(phtemp[id*(nv+1)+lev],2.0)-pow(phtemp[id*(nv+1)+lev+1],2.0));
    }
}

__global__ void rtm_dual_band (double *pressure_d   ,
                               double *Rho_d        ,
                               double *temperature_d,
                               double *fnet_up_d    ,
                               double *fnet_dn_d    ,
                               double *tau_d        ,
                               double  gravit       ,
                               double  Cp           ,
                               double *lonlat_d     ,
                               double *Altitude_d   ,
                               double *Altitudeh_d  ,
                               double *phtemp       ,
                               double *dtemp        ,
                               double *ttemp        ,
                               double *thtemp       ,
                               double timestep      ,
                               double tstar         ,
                               double planet_star_dist,
                               double radius_star   ,
                               double diff_fac      ,
                               double tlow          ,
                               double alb           ,
                               double tausw         ,
                               double taulw         ,
                               double incflx        ,
                               double ps0           ,
                               int num              ,
                               int nv               ,
                               int nvi              ,
                               double A             ){


//
//  Description:
//
//
//
//  Input: .
//
//  Output:
//

    int id = blockIdx.x * blockDim.x + threadIdx.x;
    //int id;
    //for (id=0;id<num;id++){

    double coszrs;

    // double *ph, *th, *t;
    double ps, psm;
    double pp, ptop;


    // RTM parameters
    // double bc = 5.677036E-8; // Stefan–Boltzmann constant [W m−2 K−4]

	//Star
	// double tstar = 4520; //10170;                                   // Star effective temperature [K]
  //   double planet_star_dist = 0.015*149597870.7; //0.03462*149597870.7;          // Planet-star distance [km]
  //   double radius_star = 0.667*695508; //2.362*695508;                      // Star radius [km]
  //   double resc_flx = pow(radius_star/planet_star_dist,2.0);//
  //   double incflx = resc_flx*bc*tstar*tstar*tstar*tstar;                // Incoming stellar flux [W m-2]
	// //Planet
	// double diff_fac = 0.5;         // Diffusivity factor: 0.5-1.0
	// double tlow       = 970; // 3500;        // Lower boundary temperature: upward flux coming from the planet's interior
  //   double alb      = 0.18; //0.1;         // Bond albedo
  //   double tausw    = 532.0;       // Absorption coefficient for the shortwaves
  //   double taulw    = 1064.0;      // Absorption coefficient for the longwaves
  //   double ps0      = 10000000;    // Reference bottom pressure

    double xi, xip, xim, a, b;

    if(id < num){

      for(int lev = 0; lev < nv; lev++) {
        dtemp[id*nv+lev] = 0.0;
      }
		// Calculate pressures and temperatures at interfaces
        for (int lev = 0; lev <= nv; lev++ ){
            if(lev == 0){
                psm = pressure_d[id*nv + 1] - Rho_d[id*nv + 0] * gravit * (-Altitude_d[0] - Altitude_d[1]);
                ps  = 0.5*(pressure_d[id*nv + 0] + psm);

                phtemp[id*nvi+0] = ps;
                ttemp[id*nv+0]  = temperature_d[id*nv + 0];
                thtemp[id*nvi+0] = ttemp[id*nv+0];
            }
            else if(lev == nv){
                // pp = pressure_d[id*nv + nv-2] - Rho_d[id*nv + nv-1] * gravit * (2*Altitudeh_d[nv]-Altitude_d[nv-1]-Altitude_d[nv-2]);
                pp = pressure_d[id*nv+nv-2] + (pressure_d[id*nv+nv-1]-pressure_d[id*nv+nv-2])/(Altitude_d[nv-1]-Altitude_d[nv-2])*(2*Altitudeh_d[nv]-Altitude_d[nv-1]-Altitude_d[nv-2]);
                if (pp < 0) pp = 0;  //prevents pressure at the top from becoming negative
                ptop  = 0.5*(pressure_d[id*nv + nv-1] + pp);

                phtemp[id*nvi+nv] = ptop;
                thtemp[id*nvi+nv] = temperature_d[id*nv + nv-1];
            }
            else{
                ttemp[id*nv+lev] = temperature_d[id*nv + lev];
                xi  = Altitudeh_d[lev  ] ;
                xim = Altitude_d[lev-1 ] ;
                xip = Altitude_d[lev   ] ;
                a   = (xi - xip)/(xim -xip);
                b   = (xi - xim)/(xip -xim);

                phtemp[id*nvi+lev] = pressure_d[id*nv + lev-1]*a + pressure_d[id*nv + lev]*b;
                thtemp[id*nvi+lev] = temperature_d[id*nv + lev-1]*a + ttemp[id*nv+lev]*b;
            }
            // if (id==431){
            //   printf("%f, ",pressure_d[lev]);
            // }
            // if (phtemp[id*nvi+lev]<=0){
            //   printf("%d, %d",id,lev);
            // }
        }

        // if (phtemp[id*nvi+nv]<=0){
        //   for (int lev = 0; lev < nv; lev++ ){
        //     printf("%f, ",pressure_d[id*nv+lev]);
        //   }
        //   printf("stop");
        // }
		// Cosine of the zenith angle
        coszrs = cos(lonlat_d[id*2 + 1])*cos(lonlat_d[id*2 + 0]);

		// Compute opacities
        computetau (tau_d ,
                    phtemp    ,
                    coszrs,
                    tausw ,
                    taulw ,
                    ps0   ,
                    id    ,
                    nv    );

        // for(int lev = 0; lev <=nv; lev++){
        //       fnet_up_d[id*nvi + lev] = 0.0;
        //       if(id==0) {
        //         printf("%f, ",tau_d[2*lev]);
        //       }
        // }
        // for(int lev = 0; lev <=nv; lev++) {
        //       fnet_dn_d[id*nvi + lev] = 0.0;
        //       if(id==0) printf("%f, ",tau_d[2*lev+1]);
        // }
        //
        // if (id == 0){
        //   printf("stop\n");
        // }

        if(coszrs > 0.0){

               radcsw(phtemp         ,
                      coszrs     ,
                      dtemp      ,
                      tau_d      ,
                      fnet_up_d  ,
                      fnet_dn_d  ,
                      incflx     ,
                      alb        ,
                      tausw      ,
                      ps0        ,
                      Cp         ,
                      gravit     ,
                      id         ,
                      nv         );

        }

        // if (id == 0){
        //   printf("stop\n");
        // }

        for(int lev = 0; lev <=nv; lev++)fnet_up_d[id*nvi + lev] = 0.0;
        for(int lev = 0; lev <=nv; lev++)fnet_dn_d[id*nvi + lev] = 0.0;

        radclw(phtemp       ,
               ttemp        ,
               thtemp       ,
               dtemp    ,
               tau_d    ,
               fnet_up_d,
               fnet_dn_d,
               diff_fac ,
               tlow     ,
               Cp       ,
               gravit   ,
               id       ,
               nv       );

        // if(id==0) printf("\ndT = \n");
        for(int lev = 0; lev < nv; lev++) {
          temperature_d[id*nv + lev] = ttemp[id*nv+lev] + dtemp[id*nv+lev]*timestep;
          // if (id==0) printf("%e, ",dtemp[id*nv+lev]);
          if (isnan(temperature_d[id*nv+lev])) {
            printf("%d, %d\n",id,lev);
          }
          if (temperature_d[id*nv+lev]<0) {
            printf("%d, %d\n",id,lev);
          }
        }
    }
}
