__device__ double BilinearInterpolation_dev(double q11, 
									        double q12, 
									        double q21, 
									        double q22, 
									        double x1 , 
									        double x2 , 
									        double y1 , 
									        double y2 , 
									        double x  , 
									        double y  ) {
	
    double x2x1, y2y1, x2x, y2y, yy1, xx1;
    
    x2x1 = x2 - x1;
    y2y1 = y2 - y1;
    x2x = x2 - x;
    y2y = y2 - y;
    yy1 = y - y1;
    xx1 = x - x1;
    
    return 1.0 / (x2x1 * y2y1) * (q11 * x2x * y2y + q21 * xx1 * y2y + q12 * x2x * yy1 + q22 * xx1 * yy1);
}


__device__ int locate_min_i_dev (double *array_d,
							     int     N      ,
							     double  val    ){

	int id = -1;
	if(val >= array_d[N-1]){
		id = N-2;
	}
	else if(val < array_d[0]){
		id = 0;
	}
	else{
		for(int j = 1; j < N; j++){
			if(val >= array_d[j-1] && val < array_d[j]){
				id = j-1;
				break;
			}
		}
	}	
	return id;
}

__device__ int locate_max_i_dev (double *array_d,
						         int     N      ,
					 	         double  val    ){

	int id = -1;
	if(val >= array_d[N-1]){
		id = N-1;
	}
	if(val < array_d[0]){
		id = 1;
	}
	else{
		for(int j = 1; j < N; j++){
			if(val >= array_d[j-1] && val < array_d[j]){
				id = j;
				break;
			}
		}
	}	
	return id;
	
}


__global__ void Tracers_relax_vulcan(double *tracer_d     ,
									 double *tauch4_d     ,
									 double *tauco_d      ,
									 double *tauh2o_d     ,
									 double *tauco2_d     ,
									 double *taunh3_d     ,
									 double *ch4eq_d      ,
									 double *coeq_d       ,
									 double *h2oeq_d      ,
									 double *co2eq_d      ,
									 double *nh3eq_d      ,
									 double *P_che_d      ,
									 double *T_che_d      ,
									 double *temperature_d,
									 double *pressure_d   ,
									 double *Rho_d        ,
									 double dt            ,
									 int ntr              ,
					                 int num              ){

	int id = blockIdx.x * blockDim.x + threadIdx.x;
	int nv = gridDim.y;
	int lev = blockIdx.y;
	int itr = blockIdx.z;
	
	double T, lP, P ;
	double tau, treq, tr, rho;
	
	int NT = 55 ;
	int NP = 135;
	
	int TMIN , TMAX ;
	int PMIN , PMAX ;
	
	double q11, q12, q21, q22; 
	double x1 ,	x2 , y1 , y2 ; 
	double x  ,	y;
	
	if(id < num && itr !=3){
		P  = pressure_d[id*nv + lev]/100000;
		lP = log(P);
	    T  = temperature_d[id*nv + lev];
			
		if(lP < P_che_d[0])   lP = P_che_d[0]   ;
		if(lP > P_che_d[NP-1])lP = P_che_d[NP-1];
		if(T  < T_che_d[0])   T  = T_che_d[0]   ;
		if(T  > T_che_d[NT-1])T  = T_che_d[NT-1];		    
		    
		// Find the 4 nearst points.
		TMIN = locate_min_i_dev (T_che_d,
							     NT     ,
						 	     T      );	
		TMAX = locate_max_i_dev (T_che_d,
							     NT     ,
						 	     T      );	
			
		PMIN = locate_min_i_dev  (P_che_d,
							 	  NP     ,
							 	  lP     );	
		PMAX = locate_max_i_dev  (P_che_d,
							 	  NP     ,
							 	  lP     );
			
		// Interpolate timescale.
		if(itr == 0){
			q11 = ch4eq_d[PMIN*NT + TMIN];
			q12 = ch4eq_d[PMIN*NT + TMAX];
			q21 = ch4eq_d[PMAX*NT + TMIN];
			q22 = ch4eq_d[PMAX*NT + TMAX];
		}
		else if(itr == 1){
			q11 = coeq_d[PMIN*NT + TMIN];
			q12 = coeq_d[PMIN*NT + TMAX];
			q21 = coeq_d[PMAX*NT + TMIN];
			q22 = coeq_d[PMAX*NT + TMAX];
		}	
		else if(itr == 2){
			q11 = h2oeq_d[PMIN*NT + TMIN];
			q12 = h2oeq_d[PMIN*NT + TMAX];
			q21 = h2oeq_d[PMAX*NT + TMIN];
			q22 = h2oeq_d[PMAX*NT + TMAX];
		}
		else if(itr == 3){
			q11 = co2eq_d[PMIN*NT + TMIN];
			q12 = co2eq_d[PMIN*NT + TMAX];
			q21 = co2eq_d[PMAX*NT + TMIN];
			q22 = co2eq_d[PMAX*NT + TMAX];
		}
		else if(itr == 4){
			q11 = nh3eq_d[PMIN*NT + TMIN];
			q12 = nh3eq_d[PMIN*NT + TMAX];
			q21 = nh3eq_d[PMAX*NT + TMIN];
			q22 = nh3eq_d[PMAX*NT + TMAX];
		}
			
		x1  = T_che_d[TMIN];
		x2  = T_che_d[TMAX];
		y1  = P_che_d[PMIN];
		y2  = P_che_d[PMAX];
		x   = T ;
		y   = lP;
					
		treq= BilinearInterpolation_dev(q11, 
									    q12, 
									    q21, 
									    q22, 
									    x1 , 
									    x2 , 
									    y1 , 
									    y2 , 
									    x  , 
									    y  );
		
		// Interpolate chemical equilibrium.
		if(itr == 0){
			q11 = tauch4_d[PMIN*NT + TMIN];
			q12 = tauch4_d[PMIN*NT + TMAX];
			q21 = tauch4_d[PMAX*NT + TMIN];
			q22 = tauch4_d[PMAX*NT + TMAX];
		}
		else if(itr == 1){
			q11 = tauco_d[PMIN*NT + TMIN];
			q12 = tauco_d[PMIN*NT + TMAX];
			q21 = tauco_d[PMAX*NT + TMIN];
			q22 = tauco_d[PMAX*NT + TMAX];
		}
		else if(itr == 2){
			q11 = tauh2o_d[PMIN*NT + TMIN];
			q12 = tauh2o_d[PMIN*NT + TMAX];
			q21 = tauh2o_d[PMAX*NT + TMIN];
			q22 = tauh2o_d[PMAX*NT + TMAX];	
		}
		else if(itr == 3){
			q11 = tauco2_d[PMIN*NT + TMIN];
			q12 = tauco2_d[PMIN*NT + TMAX];
			q21 = tauco2_d[PMAX*NT + TMIN];
			q22 = tauco2_d[PMAX*NT + TMAX];
		}
		else if(itr == 4){
			q11 = taunh3_d[PMIN*NT + TMIN];
			q12 = taunh3_d[PMIN*NT + TMAX];
			q21 = taunh3_d[PMAX*NT + TMIN];
			q22 = taunh3_d[PMAX*NT + TMAX];
		}
		
		x1  = T_che_d[TMIN];
		x2  = T_che_d[TMAX];
		y1  = P_che_d[PMIN];
		y2  = P_che_d[PMAX];
		x   = T ;
		y   = lP;
			
		tau = BilinearInterpolation_dev(q11, 
									    q12, 
									    q21, 
									    q22, 
									    x1 , 
									    x2 , 
									    y1 , 
									    y2 , 
									    x  , 
									    y  );
	
		rho   =   Rho_d[id*nv + lev];
		tr    =   tracer_d[id*nv*ntr + lev*ntr + itr]/rho;
		tr = (treq*(dt/tau) + tr)/(1+(dt/tau)); // Implicit
//      Fix concentrations
		if(tr > 1.0) tr = 1.0;
		if(tr < 0.0) tr = 0.0;		
				
		tracer_d[id*nv*ntr + lev*ntr + itr] = tr * rho;
	}
}

__global__ void Tracers_relax_vulcan_co2(double *tracer_d     ,
									 double *tauch4_d     ,
									 double *tauco_d      ,
									 double *tauh2o_d     ,
									 double *tauco2_d     ,
									 double *taunh3_d     ,
									 double *ch4eq_d      ,
									 double *coeq_d       ,
									 double *h2oeq_d      ,
									 double *co2eq_d      ,
									 double *nh3eq_d      ,
									 double *P_che_d      ,
									 double *T_che_d      ,
									 double *temperature_d,
									 double *pressure_d   ,
									 double *Rho_d        ,
									 double dt            ,
									 int ntr              ,
					                 int num              ){

	int id = blockIdx.x * blockDim.x + threadIdx.x;
	int nv = gridDim.y;
	int lev = blockIdx.y;
	int itr = blockIdx.z;
	
	double T, lP, P ;
	double tau, treq, treqco, treqh2o, tr, rho;
	
	int NT = 55 ;
	int NP = 135;
	
	int TMIN , TMAX ;
	int PMIN , PMAX ;
	
	double q11, q12, q21, q22; 
	double x1 ,	x2 , y1 , y2 ; 
	double x  ,	y;
	
	if(id < num && itr == 3){
		
		P  = pressure_d[id*nv + lev]/100000;
		lP = log(P);
	    T  = temperature_d[id*nv + lev];
			
		if(lP < P_che_d[0])   lP = P_che_d[0]   ;
		if(lP > P_che_d[NP-1])lP = P_che_d[NP-1];
		if(T  < T_che_d[0])   T  = T_che_d[0]   ;
		if(T  > T_che_d[NT-1])T  = T_che_d[NT-1];		    
		    
		// Find the 4 nearst points.
		TMIN = locate_min_i_dev (T_che_d,
							     NT     ,
						 	     T      );	
		TMAX = locate_max_i_dev (T_che_d,
							     NT     ,
						 	     T      );	
			
		PMIN = locate_min_i_dev  (P_che_d,
							 	  NP     ,
							 	  lP     );	
		PMAX = locate_max_i_dev  (P_che_d,
							 	  NP     ,
							 	  lP     );
			
		// Interpolate timescale.
		x1  = T_che_d[TMIN];
		x2  = T_che_d[TMAX];
		y1  = P_che_d[PMIN];
		y2  = P_che_d[PMAX];
		x   = T ;
		y   = lP;
		q11 = coeq_d[PMIN*NT + TMIN];
		q12 = coeq_d[PMIN*NT + TMAX];
		q21 = coeq_d[PMAX*NT + TMIN];
		q22 = coeq_d[PMAX*NT + TMAX];
	
		treqco= BilinearInterpolation_dev(q11, 
									    q12, 
									    q21, 
									    q22, 
									    x1 , 
									    x2 , 
									    y1 , 
									    y2 , 
									    x  , 
									    y  );
		
		q11 = h2oeq_d[PMIN*NT + TMIN];
		q12 = h2oeq_d[PMIN*NT + TMAX];
		q21 = h2oeq_d[PMAX*NT + TMIN];
		q22 = h2oeq_d[PMAX*NT + TMAX];

		treqh2o= BilinearInterpolation_dev(q11, 
										    q12, 
										    q21, 
										    q22, 
										    x1 , 
										    x2 , 
										    y1 , 
										    y2 , 
										    x  , 
										    y  );			
			
		q11 = co2eq_d[PMIN*NT + TMIN];
		q12 = co2eq_d[PMIN*NT + TMAX];
		q21 = co2eq_d[PMAX*NT + TMIN];
		q22 = co2eq_d[PMAX*NT + TMAX];
		
		treq= BilinearInterpolation_dev(q11, 
									    q12, 
									    q21, 
									    q22, 
									    x1 , 
									    x2 , 
									    y1 , 
									    y2 , 
									    x  , 
									    y  );
		
		rho   =   Rho_d[id*nv + lev];
		double cov = tracer_d[id*nv*ntr + lev*ntr + 1]/rho;
		double h2ov= tracer_d[id*nv*ntr + lev*ntr + 2]/rho;
		double factor = (cov*h2ov)/(treqco*treqh2o);
		
		// Interpolate chemical equilibrium.
		q11 = tauco2_d[PMIN*NT + TMIN];
		q12 = tauco2_d[PMIN*NT + TMAX];
		q21 = tauco2_d[PMAX*NT + TMIN];
		q22 = tauco2_d[PMAX*NT + TMAX];
		
		x1  = T_che_d[TMIN];
		x2  = T_che_d[TMAX];
		y1  = P_che_d[PMIN];
		y2  = P_che_d[PMAX];
		x   = T ;
		y   = lP;
			
		tau = BilinearInterpolation_dev(q11, 
									    q12, 
									    q21, 
									    q22, 
									    x1 , 
									    x2 , 
									    y1 , 
									    y2 , 
									    x  , 
									    y  );
	
		tr    =   tracer_d[id*nv*ntr + lev*ntr + itr]/rho;
		tr = (factor*treq*(dt/tau) + tr)/(1+(dt/tau)); // Implicit
//      Fix concentrations
		if(tr > 1.0) tr = 1.0;
		if(tr < 0.0) tr = 0.0;		
				
		tracer_d[id*nv*ntr + lev*ntr + itr] = tr * rho;
	}
}

