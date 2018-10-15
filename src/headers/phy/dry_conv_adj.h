#include <math.h>

template<int nv>
__global__ void dry_conv_adj (double *Pressure_d    , // Pressure [Pa] 
                              double *Temperature_d , // Temperature [K]
                              double *Rho_d         , // Density [m^3/kg]
                              double Cp             , // Specific heat capacity [J/kg/K]
                              double Rd             , // Gas constant [J/kg/K]
                              double Gravit         , // Gravity [m/s^2]
                              double *Altitude_d    , // Altitudes of the layers
                              double *Altitudeh_d   , // Altitudes of the interfaces
                              int num               ){ // Number of columns
	
//
//  Description: Mixes entropy vertically on statically unstable columns
//

	int id = blockIdx.x * blockDim.x + threadIdx.x;

	// Local arrays
	double theta[nv];
	double ph[nv+1];
	
	// stability threshold
	double stable = 0.0;

	double ps, psm;
	double pp, ptop; 	

	double xi, xip, xim, a, b;

	if(id < num){
		for (int lev = 0; lev <= nv; lev++ ){
			if (lev == 0){
				psm = Pressure_d[id*nv + 1] - Rho_d[id*nv + 0] * Gravit * (-Altitude_d[0] - Altitude_d[1]);
				ps  = 0.5*(Pressure_d[id*nv + 0] + psm);
				ph[0] = ps;		
			}
			else if (lev == nv){

				pp = Pressure_d[id*nv + nv-2] - Rho_d[id*nv + nv-1] * Gravit * (2*Altitudeh_d[nv]-Altitude_d[nv-1]-Altitude_d[nv-2]);	
                if(pp < 0) pp = 0;
				ptop  = 0.5*(Pressure_d[id*nv + nv-1] + pp);
				ph[nv] = ptop;		
			}
			else{
				xi  = Altitudeh_d[lev  ] ;
				xim = Altitude_d[lev-1 ] ;
				xip = Altitude_d[lev   ] ;
				a   = (xi - xip)/(xim -xip);
				b   = (xi - xim)/(xip -xim);
				ph[lev] = Pressure_d[id*nv + lev-1]*a + Pressure_d[id*nv + lev]*b;
			}
		}

		// Compute Potential Temperature		
		for(int lev = 0; lev < nv; lev++){
			theta[lev] = Temperature_d[id*nv +lev]*pow(ph[0]/Pressure_d[id*nv + lev], Rd/Cp);
		}

		bool  done_col = false;
		while(done_col == false){ // Unstable  column?
			int top = 0 ;
			int bot = nv-1; 

			for(int lev = 0; lev < nv-1; lev++){
				if(theta[lev+1]-theta[lev] < stable){
					if(bot > lev) bot = lev;
				}
			}

			if(bot > nv - 2) done_col = true;

			for(int lev = bot; lev < nv-1; lev++){			
				if(theta[lev+1]-theta[lev] > stable){
					top = lev;
					goto stop1;
				}
				else{
					top = nv -1;
				}
			}
			stop1:
			if(bot < nv -1){
				double h   =  0.0;//Enthalpy;
				double sum =  0.0;

				for(int lev = bot; lev <= top; lev++){				
					double pu = ph[lev+1];
					double pl = ph[lev  ];
					double pi = pow(Pressure_d[id*nv + lev]/ph[0],Rd/Cp);
					double deltap = pl - pu;

					h = h + theta[lev]*pi*deltap;
					sum = sum + pi * deltap;
				}
				double thnew = h/sum;

				int extend = 0;

				if(bot > 0){
					if((thnew-theta[bot-1]) < stable){
						bot = bot - 1;
						extend = 1      ;
					}
				}
				
				if(top < nv-1){
					if((theta[top+1] - thnew) < stable){
						top = top + 1;
						extend = 1.0;
					}

				}

				if(extend != 0) goto stop1;
				
				for(int lev = bot; lev <= top; lev++){
					theta[lev] = thnew;
				}
			}
		}

		// Compute Temperature & pressure
		for(int lev = 0; lev < nv; lev++){
			Temperature_d[id*nv +lev] = theta[lev]*pow(Pressure_d[id*nv + lev]/ph[0], Rd/Cp);
			Pressure_d[id*nv + lev] = Temperature_d[id*nv + lev]*Rd*Rho_d[id*nv + lev];
		}
	}
	
}