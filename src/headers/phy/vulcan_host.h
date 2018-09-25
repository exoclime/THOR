__host__ double BilinearInterpolation_host(double q11, 
							               double q12, 
								           double q21, 
								           double q22, 
								           double x1 , 
								           double x2 , 
								           double y1 , 
								           double y2 , 
								           double x  , 
								           double y  ) {
	
    float x2x1, y2y1, x2x, y2y, yy1, xx1;
    
    x2x1 = x2 - x1;
    y2y1 = y2 - y1;
    x2x  = x2 - x ;
    y2y  = y2 - y ;
    yy1  = y  - y1;
    xx1  = x  - x1;
    
    return 1.0 / (x2x1 * y2y1) * (q11 * x2x * y2y + q21 * xx1 * y2y + q12 * x2x * yy1 + q22 * xx1 * yy1);
}


__host__ int locate_min_i_host (double *array_d,
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

__host__ int locate_max_i_host (double *array_d,
						   	    int     N      ,
						   	    double  val    ){

	int id = -1;
	if(val >= array_d[N-1]){
		id = N-1;
	}
	else if(val < array_d[0]){
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

__host__ double Compute_tracer_host (double *treq       ,
					       	   	     double *P_che      ,
					       	   	     double *T_che      ,
					       	   	     double temperature ,
					       	   	     double pressure    ){

//
//  Description:
//
//  
//
//  Input: .
//
//  Output:
//


	double T, lP, P ;
	
	int NT = 55;
	int NP = 135;
	
	int TMIN , TMAX ;
	int PMIN , PMAX ;
	
	double q11, q12, q21, q22; 
	double x1 ,	x2 , y1 , y2 ; 
	double x  ,	y;
	double tr;

		
	P  = pressure/100000.0;
	lP = log(P);
	T  = temperature;
			
	if(lP < P_che[0])   lP = P_che[0]   ;
	if(lP > P_che[NP-1])lP = P_che[NP-1];
	if(T  < T_che[0])   T  = T_che[0]   ;
	if(T  > T_che[NT-1])T  = T_che[NT-1];
		    	    
	// Find the 4 nearst points.
	TMIN = locate_min_i_host (T_che,
							  NT     ,
						 	  T      );	
	TMAX = locate_max_i_host (T_che,
							  NT     ,
						 	  T      );	
			
	PMIN = locate_min_i_host (P_che,
							  NP     ,
						 	  lP     );	
	PMAX = locate_max_i_host (P_che,
							  NP     ,
						 	  lP     );	
						
	// Interpolate timescale.
	q11 = treq[PMIN*NT + TMIN];
	q12 = treq[PMIN*NT + TMAX];
	q21 = treq[PMAX*NT + TMIN];
	q22 = treq[PMAX*NT + TMAX];
	x1  = T_che[TMIN];
	x2  = T_che[TMAX];
	y1  = P_che[PMIN];
	y2  = P_che[PMAX];
	x   = T ;
	y   = lP;
	
	tr= BilinearInterpolation_host(q11, 
								   q12, 
								   q21, 
								   q22, 
								   x1 , 
								   x2 , 
								   y1 , 
								   y2 , 
								   x  , 
								   y  );
		
	return tr;
}