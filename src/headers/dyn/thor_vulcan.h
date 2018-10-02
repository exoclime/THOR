//
// 
//
//
// Description:
//        Contains subroutines call in thor_drive:
//							-
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
// Current Code Owner: ...
//
// History:
// Version Date       Comment
// ======= ====       =======
//
// 1.0    ... 
//
////////////////////////////////////////////////////////////////////////
template <int NX, int NY>
__global__ void Tracer_Eq (double *tracers_d  ,
			   	   	       double *tracerk_d  ,
			   	   	   	   double *Rho_d      ,
						   double *Rhok_d     ,
						   double *Mh_d       ,
						   double *Mhk_d      ,
						   double *Wh_d       ,
						   double *Whk_d      ,
						   double *difftr_d   ,
						   double *div_d      ,
						   double *Altitude_d ,
						   double *Altitudeh_d,
						   double  A          ,
						   double  dt         ,     
						   int *maps_d        ,
					 	   int     ntr        ,
						   int nl_region      ,
						   int DeepModel      ){

		int x   = threadIdx.x;
		int y   = threadIdx.y;
		int ib  = blockIdx.x;
		int nv  = gridDim.y;
		int lev = blockIdx.y;
		int itr = blockIdx.z;
		
		int pt1, pt2, pt3, pt4, pt5, pt6;
		double div0, div1, div2, div3, div4, div5, div6;
		int nhl = nl_region + 2;
		int nhl2 = nhl*nhl;

		int ir = (y + 1)*nhl + x + 1;   // Region index
		int iri, ir2, twot;

		__shared__ double nflxtr_s[NX*NY];
		__shared__ double v1_s[3 * (NX + 2)*(NY + 2)];
		__shared__ double tr_s[(NX + 2)*(NY + 2)];

		int pent_ind = 0;
		int ig, id;

		double xi, xim1, xip1;
		double a, b;
		double alt, rscale;
		double dtr;
		double trht, trhl;
		double dz;
		double whl, wht;
		double r, tr;
		double r2p, r2m, r2l;
		double altht, althl;
		double dtr_dalt;

		// Load shared memory
		ig = maps_d[ib*nhl2 + ir];
		id = ig;
		if (x == 0 && y == 0) if (maps_d[ib * nhl2] == -1) pent_ind = 1;
		v1_s[ir * 3 + 0] = Mh_d[ig * 3 * nv + lev * 3 + 0] + Mhk_d[ig * 3 * nv + lev * 3 + 0];
		v1_s[ir * 3 + 1] = Mh_d[ig * 3 * nv + lev * 3 + 1] + Mhk_d[ig * 3 * nv + lev * 3 + 1];
		v1_s[ir * 3 + 2] = Mh_d[ig * 3 * nv + lev * 3 + 2] + Mhk_d[ig * 3 * nv + lev * 3 + 2];

		tr_s[ir] = tracerk_d[ig*nv*ntr + lev*ntr + itr]/Rhok_d[ig*nv + lev];
		
		///////////////////////////////
		//////////// Halo /////////////
		///////////////////////////////
		if (x == 0) {
			ir2 = (y + 1) * nhl + x;
			ig = maps_d[ib * nhl2 + ir2];
			v1_s[ir2 * 3 + 0] = Mh_d[ig * 3 * nv + lev * 3 + 0] + Mhk_d[ig * 3 * nv + lev * 3 + 0];
			v1_s[ir2 * 3 + 1] = Mh_d[ig * 3 * nv + lev * 3 + 1] + Mhk_d[ig * 3 * nv + lev * 3 + 1];
			v1_s[ir2 * 3 + 2] = Mh_d[ig * 3 * nv + lev * 3 + 2] + Mhk_d[ig * 3 * nv + lev * 3 + 2];
			tr_s[ir2] = tracerk_d[ig*nv*ntr + lev*ntr + itr]/Rhok_d[ig*nv + lev];
		}
		if (x == nhl - 3){
			ir2 = (y + 1) * nhl + x + 2;
			ig = maps_d[ib * nhl2 + ir2];
			v1_s[ir2 * 3 + 0] = Mh_d[ig * 3 * nv + lev * 3 + 0] + Mhk_d[ig * 3 * nv + lev * 3 + 0];
			v1_s[ir2 * 3 + 1] = Mh_d[ig * 3 * nv + lev * 3 + 1] + Mhk_d[ig * 3 * nv + lev * 3 + 1];
			v1_s[ir2 * 3 + 2] = Mh_d[ig * 3 * nv + lev * 3 + 2] + Mhk_d[ig * 3 * nv + lev * 3 + 2];
			tr_s[ir2] = tracerk_d[ig*nv*ntr + lev*ntr + itr]/Rhok_d[ig*nv + lev];
		}
		if (y == 0){
			twot = 1;
			ir2 = y * nhl + (x + 1);
			if (x == 0) twot = 2;

			for (int k = 0; k < twot; k++){
				if (k == 1) ir2 = y * nhl + x;
				ig = maps_d[ib * nhl2 + ir2];

				if (ig >= 0){
					v1_s[ir2 * 3 + 0] = Mh_d[ig * 3 * nv + lev * 3 + 0] + Mhk_d[ig * 3 * nv + lev * 3 + 0];
					v1_s[ir2 * 3 + 1] = Mh_d[ig * 3 * nv + lev * 3 + 1] + Mhk_d[ig * 3 * nv + lev * 3 + 1];
					v1_s[ir2 * 3 + 2] = Mh_d[ig * 3 * nv + lev * 3 + 2] + Mhk_d[ig * 3 * nv + lev * 3 + 2];
					tr_s[ir2] = tracerk_d[ig*nv*ntr + lev*ntr + itr]/Rhok_d[ig*nv + lev];
				}
				else{
					v1_s[ir2 * 3 + 0] = 0.0;
					v1_s[ir2 * 3 + 1] = 0.0;
					v1_s[ir2 * 3 + 2] = 0.0;
					tr_s[ir2] = 0.0;
				}
			}
		}
		if (y == nhl - 3) {
			twot = 1;
			ir2 = (y + 2) * nhl + (x + 1);
			if (x == nhl - 3) twot = 2;
			for (int k = 0; k < twot; k++){
				if (k == 1) ir2 = (y + 2) * nhl + (x + 2);
				ig = maps_d[ib * nhl2 + ir2];
				v1_s[ir2 * 3 + 0] = Mh_d[ig * 3 * nv + lev * 3 + 0] + Mhk_d[ig * 3 * nv + lev * 3 + 0];
				v1_s[ir2 * 3 + 1] = Mh_d[ig * 3 * nv + lev * 3 + 1] + Mhk_d[ig * 3 * nv + lev * 3 + 1];
				v1_s[ir2 * 3 + 2] = Mh_d[ig * 3 * nv + lev * 3 + 2] + Mhk_d[ig * 3 * nv + lev * 3 + 2];
				tr_s[ir2] = tracerk_d[ig*nv*ntr + lev*ntr + itr]/Rhok_d[ig*nv + lev];
			}
		}
		__syncthreads();
		//////////////////////////////////////////////

		iri = (y  )*nl_region + x;

		pt1 = (y + 2)*nhl + x + 1;
		pt2 = (y + 2)*nhl + x + 2;
		pt3 = (y + 1)*nhl + x + 2;
		pt4 = (y    )*nhl + x + 1;
	 	pt5 = (pent_ind)*((y + 1)*nhl + x) + (!pent_ind)*((y   )*nhl + x);
		pt6 = (y + 1)*nhl + x;

		altht = Altitudeh_d[lev + 1];
		althl = Altitudeh_d[lev];

		if (DeepModel == 1){
			alt = Altitude_d[lev];
			r2p = pow(altht + A, 2.0);
			r2m = pow(alt + A, 2.0);
			r2l = pow(althl + A, 2.0);
			rscale = A / (alt + A);
		}
		else{
			r2p = 1.0;
			r2m = 1.0;
			r2l = 1.0;
			rscale = 1.0;
		}

		nflxtr_s[iri] = 0.0;

		for (int k = 0; k < 3; k++){

			div0 = div_d[id * 7 * 3 + 3 * 0 + k];
			div1 = div_d[id * 7 * 3 + 3 * 1 + k];
			div2 = div_d[id * 7 * 3 + 3 * 2 + k];
			div3 = div_d[id * 7 * 3 + 3 * 3 + k];
			div4 = div_d[id * 7 * 3 + 3 * 4 + k];
			div5 = div_d[id * 7 * 3 + 3 * 5 + k];
			div6 = div_d[id * 7 * 3 + 3 * 6 + k];

			nflxtr_s[iri] += rscale*(div0 * v1_s[ir  * 3 + k]* tr_s[ir]  +
									 div1 * v1_s[pt1 * 3 + k]* tr_s[pt1] +
									 div2 * v1_s[pt2 * 3 + k]* tr_s[pt2] +
									 div3 * v1_s[pt3 * 3 + k]* tr_s[pt3] +
									 div4 * v1_s[pt4 * 3 + k]* tr_s[pt4] +
									 div5 * v1_s[pt5 * 3 + k]* tr_s[pt5] +
									 div6 * v1_s[pt6 * 3 + k]* tr_s[pt6]);
		}
		
		if(lev==0){
			trhl = 0.0;
			xi  = Altitudeh_d[lev +1 ] ;
			xim1= Altitude_d[lev ] ;
			xip1= Altitude_d[lev+1   ] ;
			a = (xi - xip1)/(xim1 - xip1);
			b = (xi - xim1)/(xip1 - xim1);				
			trht  = tracerk_d[id*nv*ntr + lev*ntr + itr]     / Rhok_d[id*nv + lev]   *a + 
					tracerk_d[id*nv*ntr + (lev+1)*ntr + itr] / Rhok_d[id*nv + lev+1] *b;
		}
		else if(lev==nv-1){
			xi  = Altitudeh_d[lev ] ;
			xim1= Altitude_d[lev-1 ] ;
			xip1= Altitude_d[lev   ] ;
			a = (xi - xip1)/(xim1 - xip1);
			b = (xi - xim1)/(xip1 - xim1);				
			trhl  = tracerk_d[id*nv*ntr + (lev-1)*ntr + itr] / Rhok_d[id*nv + lev - 1] * a + 
					tracerk_d[id*nv*ntr + lev*ntr + itr] / Rhok_d[id*nv + lev] * b;		
			trht = 0.0;
		}
		else{
			xi  = Altitudeh_d[lev ] ;
			xim1= Altitude_d[lev-1 ] ;
			xip1= Altitude_d[lev   ] ;
			a = (xi - xip1)/(xim1 - xip1);
			b = (xi - xim1)/(xip1 - xim1);				
			trhl  = tracerk_d[id*nv*ntr + (lev-1)*ntr + itr] / Rhok_d[id*nv + lev-1] * a + 
					tracerk_d[id*nv*ntr + lev*ntr + itr] / Rhok_d[id*nv + lev] * b;
			
			xi  = Altitudeh_d[lev +1 ] ;
			xim1= Altitude_d[lev ] ;
			xip1= Altitude_d[lev+1   ] ;
			a = (xi - xip1)/(xim1 - xip1);
			b = (xi - xim1)/(xip1 - xim1);				
			trht  = tracerk_d[id*nv*ntr + lev*ntr + itr] / Rhok_d[id*nv + lev] * a + 
					tracerk_d[id*nv*ntr + (lev+1)*ntr + itr] / Rhok_d[id*nv + lev +1] * b;
		}
			
		if (lev == 0){
			whl= 0.0;
			wht= Wh_d[id*(nv + 1) + lev + 1] + Whk_d[id*(nv + 1) + lev + 1];
		}
		else{
			whl= Wh_d[id*(nv + 1) + lev    ] + Whk_d[id*(nv + 1) + lev    ];
			wht= Wh_d[id*(nv + 1) + lev + 1] + Whk_d[id*(nv + 1) + lev + 1];
		}

		dz = altht - althl;
		dtr_dalt = (wht*trht*r2p - whl*trhl*r2l)/(dz*r2m);
		dtr = - ( nflxtr_s[iri] + dtr_dalt)*dt;
		r   = Rhok_d[id*nv + lev] + Rho_d[id*nv + lev];			
		tr  = tr_s[ir]*r + dtr;
	
		tracers_d[id*nv*ntr + lev*ntr + itr] = tr - tr_s[ir] * Rhok_d[id*nv + lev];

		tracers_d[id*nv*ntr + lev*ntr + itr] = tracers_d[id*nv*ntr + lev*ntr + itr] + difftr_d[id*nv*ntr + lev*ntr + itr]*dt;	
}
	
template <int NN>
__global__ void Tracer_Eq_Poles(double *tracers_d  ,
   	       	   	   	   	   	   	double *tracerk_d  ,
								double *Rho_d      ,
								double *Rhok_d     ,
								double *Mh_d       ,
								double *Mhk_d      ,
								double *Wh_d       ,
								double *Whk_d      ,
								double *difftr_d   ,
								double *div_d      ,
								double *Altitude_d ,
								double *Altitudeh_d,
								double  A          ,
								double  dt         ,
								int *point_local_d ,
							 	int    ntr         ,
								int    num         ,
								int    nv          ,
								int DeepModel      ){

	int id = blockIdx.x * blockDim.x + threadIdx.x;
	id += num - 2; 	// Poles
	
	int itr = blockIdx.z;
	
	

	__shared__ double div_p[3 * 7];
	__shared__ double v1_p[3 * NN];
	__shared__ double tr_p[NN*5];
	__shared__ int local_p[NN];
	
	double nflxtr_p;

	double alt, rscale;
	double wht, whl;
	double trht, trhl;
	double dz;
	double r;
	double altht, althl;
	double xi, xim1, xip1, a, b;
	double dtr_dalt, dtr, tr;
	double r2p, r2m, r2l;

	if (id < num){
		for (int i = 0; i < 5; i++)local_p[i] = point_local_d[id * 6 + i];
		for (int i = 0; i < 7; i++) for (int k = 0; k < 3; k++) div_p[i * 3 + k] = div_d[id * 7 * 3 + i * 3 + k];
	
		for (int lev = 0; lev < nv; lev++){
			v1_p[0] = Mh_d[id * 3 * nv + lev * 3 + 0] + Mhk_d[id * 3 * nv + lev * 3 + 0];
			v1_p[1] = Mh_d[id * 3 * nv + lev * 3 + 1] + Mhk_d[id * 3 * nv + lev * 3 + 1];
			v1_p[2] = Mh_d[id * 3 * nv + lev * 3 + 2] + Mhk_d[id * 3 * nv + lev * 3 + 2];
			tr_p[0] = tracerk_d[id*nv*ntr + lev*ntr + itr]/ Rhok_d[id*nv + lev];
			for (int i = 1; i < 6; i++){
				v1_p[i * 3 + 0] = Mh_d[local_p[i - 1] * 3 * nv + lev * 3 + 0] + Mhk_d[local_p[i - 1] * 3 * nv + lev * 3 + 0];
				v1_p[i * 3 + 1] = Mh_d[local_p[i - 1] * 3 * nv + lev * 3 + 1] + Mhk_d[local_p[i - 1] * 3 * nv + lev * 3 + 1];
				v1_p[i * 3 + 2] = Mh_d[local_p[i - 1] * 3 * nv + lev * 3 + 2] + Mhk_d[local_p[i - 1] * 3 * nv + lev * 3 + 2];
				tr_p[i] = tracerk_d[local_p[i - 1] * nv * ntr + lev * ntr + itr]/ Rhok_d[local_p[i - 1]*nv + lev];
			}

			if (lev == 0){
				altht = Altitudeh_d[lev + 1];
				althl = Altitudeh_d[lev];
			}

			alt = Altitude_d[lev];

			if (DeepModel == 1){
				r2p = pow(altht + A, 2.0);
				r2m = pow(alt + A, 2.0);
				r2l = pow(althl + A, 2.0);
				rscale = A / (alt + A);
			}
			else{
				r2p = 1.0;
				r2m = 1.0;
				r2l = 1.0;
				rscale = 1.0;
			}

			nflxtr_p = 0.0;

			for (int k = 0; k < 3; k++){

				nflxtr_p += rscale*(div_p[3 * 0 + k] * v1_p[0 * 3 + k] * tr_p[0] +
									div_p[3 * 1 + k] * v1_p[1 * 3 + k] * tr_p[1] +
									div_p[3 * 2 + k] * v1_p[2 * 3 + k] * tr_p[2] +
									div_p[3 * 3 + k] * v1_p[3 * 3 + k] * tr_p[3] +
									div_p[3 * 4 + k] * v1_p[4 * 3 + k] * tr_p[4] +
									div_p[3 * 5 + k] * v1_p[5 * 3 + k] * tr_p[5]);
			}
			
			if (lev == 0){
				whl = 0.0;
				wht = Wh_d[id*(nv + 1) + lev + 1] + Whk_d[id*(nv + 1) + lev + 1];
			}

			if(lev==0){
				trhl = 0.0;
				xi  = Altitudeh_d[lev +1 ] ;
				xim1= Altitude_d[lev ] ;
				xip1= Altitude_d[lev+1   ] ;
				a = (xi - xip1)/(xim1 - xip1);
				b = (xi - xim1)/(xip1 - xim1);				
				trht  = tracerk_d[id*nv*ntr + lev*ntr + itr] / Rhok_d[id*nv + lev] * a + 
						tracerk_d[id*nv*ntr + (lev+1)*ntr + itr] / Rhok_d[id*nv + lev] * b;
			}
			else if(lev==nv-1){
				xi  = Altitudeh_d[lev ] ;
				xim1= Altitude_d[lev-1 ] ;
				xip1= Altitude_d[lev   ] ;
				a = (xi - xip1)/(xim1 - xip1);
				b = (xi - xim1)/(xip1 - xim1);				
				trhl  = tracerk_d[id*nv*ntr + (lev-1)*ntr + itr] / Rhok_d[id*nv + lev-1] * a + 
						tracerk_d[id*nv*ntr + lev*ntr + itr] / Rhok_d[id*nv + lev] * b;		
				trht = 0.0;
			}
			else{
				xi  = Altitudeh_d[lev ] ;
				xim1= Altitude_d[lev-1 ] ;
				xip1= Altitude_d[lev   ] ;
				a = (xi - xip1)/(xim1 - xip1);
				b = (xi - xim1)/(xip1 - xim1);				
				trhl  = tracerk_d[id*nv*ntr + (lev-1)*ntr + itr] / Rhok_d[id*nv + lev-1] * a + 
						tracerk_d[id*nv*ntr + lev*ntr + itr] / Rhok_d[id*nv + lev] * b;
					
				xi  = Altitudeh_d[lev +1 ] ;
				xim1= Altitude_d[lev ] ;
				xip1= Altitude_d[lev+1   ] ;
				a = (xi - xip1)/(xim1 - xip1);
				b = (xi - xim1)/(xip1 - xim1);				
				trht  = tracerk_d[id*nv*ntr + lev*ntr + itr] / Rhok_d[id*nv + lev] * a + 
						tracerk_d[id*nv*ntr + (lev+1)*ntr + itr] / Rhok_d[id*nv + lev+1] * b;
			}
						
			dz = altht - althl;
			dtr_dalt = (wht*trht*r2p - whl*trhl*r2l)/(dz*r2m);
			dtr = -(nflxtr_p + dtr_dalt)*dt;
			r   = Rhok_d[id*nv + lev] + Rho_d[id*nv + lev];			
			tr  = tr_p[0]*r + dtr;
			
			tracers_d[id*nv*ntr + lev*ntr + itr] = tr - tr_p[0] * Rhok_d[id*nv + lev];

			tracers_d[id*nv*ntr + lev*ntr + itr] = tracers_d[id*nv*ntr + lev*ntr + itr] + difftr_d[id*nv*ntr + lev*ntr + itr]*dt;	

			if (lev != nv - 1){
				althl = altht;
				altht = Altitudeh_d[lev + 2];
				whl = wht;
				wht = Wh_d[id*(nv + 1) + lev + 2] + Whk_d[id*(nv + 1) + lev + 2];
			}
		}
	}
}


template <int NX, int NY>
__global__ void Tracer_Eq_Diffusion (double * difftr_d     , //
							  	     double * diff_d       , //
							  	     double * tracer_d     , //
							  	     double * Rho_d        , //
							  	     double * areasTr_d    , //
							  	     double * nvecoa_d     , //
							  	     double * nvecti_d     , //
							  	     double * nvecte_d     , //
							  	     double * K_d          , //
							  	     double * Altitude_d   , //
							  	     double   A            , //
							  	     int    * maps_d       , //
							 	     int      ntr          , //
							  	     int      nl_region    , //
							  	     int      laststep     , //
							  	     int      DeepModel    ){//

	int x = threadIdx.x;
	int y = threadIdx.y;
	int ib = blockIdx.x;
	int nv = gridDim.y;
	int lev = blockIdx.y;
	int var = blockIdx.z;

	int nhl = nl_region + 2;
	int nhl2 = nhl*nhl;
	int pt1, pt2, pt3, twot;

	double alt;
	double rscale;
	double sdiff, vdiff;
	double AT1, AT2;
	double o3 = 1.0 / 3.0;
	double o6 = 1.0 / 6.0;
	double lap, lap1, lap2, lap3, lap4, lap5, lap6;
	double lapx1, lapx2,
		   lapy1, lapy2,
		   lapz1, lapz2;

	int ir = (y + 1) * nhl + (x + 1);   // Region index
	int ir2, id, jp1, jp2;

	/////////////////////////////////////////
	__shared__ double a_s[(NX + 2)*(NY + 2)];
	__shared__ double Rho_s[(NX + 2)*(NY + 2)];
	/////////////////////////////////////////

	bool pent_ind;
	int ig       ;

	ig = maps_d[ib * nhl2 + ir];
	id = ig;
	pent_ind = 0;
	if (x == 0 && y == 0) if(maps_d[ib * nhl2] == - 1) pent_ind = 1;
	int n_faces = (pent_ind ? 5 : 6);
		
	Rho_s[ir] = Rho_d[ig * nv + lev];

	///////////////////////////////
	//////////// Halo /////////////
	///////////////////////////////
	if (x == 0) {
		ir2 = (y + 1) * nhl + x;
		ig = maps_d[ib * nhl2 + ir2];
		Rho_s[ir2] = Rho_d[ig * nv + lev];
	}
	if (x == nhl - 3){
		ir2 = (y + 1) * nhl + x + 2;
		ig = maps_d[ib * nhl2 + ir2];
		Rho_s[ir2] = Rho_d[ig * nv + lev];
	}
	if (y == 0){
		twot = 1;
		ir2 = y * nhl + (x + 1);
		if (x == 0) twot = 2;
		for (int k = 0; k < twot; k++){
			if (k == 1) ir2 = y * nhl + x;
			ig = maps_d[ib * nhl2 + ir2];
			if (ig >= 0) Rho_s[ir2] = Rho_d[ig * nv + lev];
			else         Rho_s[ir2] = 0.0;
		}
	}
	if (y == nhl - 3) {
		twot = 1;
		ir2 = (y + 2) * nhl + (x + 1);
		if (x == nhl - 3) twot = 2;
		for (int k = 0; k < twot; k++){
			if (k == 1) ir2 = (y + 2) * nhl + (x + 2);
			ig = maps_d[ib * nhl2 + ir2];
			Rho_s[ir2] = Rho_d[ig * nv + lev];
		}
	}
	__syncthreads();
	//////////////////////////////////////////////
	if (laststep) a_s[ir] = diff_d[id * nv * 6 + lev * 6 + var];
	else a_s[ir] = tracer_d[id*nv*ntr + lev*ntr + var];
	
	if (!laststep) a_s[ir] = a_s[ir] / Rho_s[ir];
	
	///////////////////////////////
	//////////// Halo /////////////
	///////////////////////////////
	if (x == 0) {
		ir2 = (y + 1) * nhl + x;
		ig = maps_d[ib * nhl2 + ir2];
		if (laststep) a_s[ir2]   = diff_d[ig * nv * 6 + lev * 6 + var];
		else a_s[ir2] = tracer_d[ig*nv*ntr + lev*ntr + var];

		if (!laststep) a_s[ir2] = a_s[ir2] / Rho_s[ir2];
	}
	if (x == nhl - 3){
		ir2 = (y + 1) * nhl + x + 2;
		ig = maps_d[ib * nhl2 + ir2];
		if (laststep) a_s[ir2]   = diff_d[ig * nv * 6 + lev * 6 + var];
		else a_s[ir2] = tracer_d[ig*nv*ntr + lev*ntr + var];

		if (!laststep) a_s[ir2] = a_s[ir2] / Rho_s[ir2];
	}
	if (y == 0){
		twot = 1;
		ir2 = y * nhl + (x + 1);
		if (x == 0) twot = 2;
		for (int k = 0; k < twot; k++){
			if (k == 1) ir2 = y * nhl + x;
			ig = maps_d[ib * nhl2 + ir2];
			if (ig >= 0){
				if (laststep) a_s[ir2] = diff_d[ig * nv * 6 + lev * 6 + var];
				else a_s[ir2] = tracer_d[ig*nv*ntr + lev*ntr + var];

				if (!laststep) a_s[ir2] = a_s[ir2] / Rho_s[ir2];
			}
			else a_s[ir2] = 0.0;
		}
	}
	if (y == nhl - 3) {
		twot = 1;
		ir2 = (y + 2) * nhl + (x + 1);
		if (x == nhl - 3) twot = 2;
		for (int k = 0; k < twot; k++){
			if (k == 1) ir2 = (y + 2) * nhl + (x + 2);
			ig = maps_d[ib * nhl2 + ir2];
			if (laststep) a_s[ir2] = diff_d[ig * nv * 6 + lev * 6 + var];
			else a_s[ir2] = tracer_d[ig*nv*ntr + lev*ntr + var];

			if (!laststep) a_s[ir2] = a_s[ir2] / Rho_s[ir2];
		}
	}
	__syncthreads();
	//////////////////////////////////////////////

	if (DeepModel == 1){
		alt = Altitude_d[lev];
		rscale = A / (alt + A);
	}
	else  rscale = 1.0;

	if (laststep)sdiff = -K_d[lev]     ;

	lap = 0.0;

	for (int j = 0; j < 6; j++){

		jp1 = (j + 1) % n_faces;
		jp2 = (j + 2) % n_faces;

		if (j == 0){
			pt1 = (y + 2)*nhl + x + 1;
			pt2 = (y + 2)*nhl + x + 2;
			pt3 = (y + 1)*nhl + x + 2;
		}
		else if (j == 1){
			pt1 = (y + 2)*nhl + x + 2;
			pt2 = (y + 1)*nhl + x + 2;
			pt3 = (y    )*nhl + x + 1;
		}
		else if (j == 2){
			pt1 = (y + 1)*nhl + x + 2;
			pt2 = (y    )*nhl + x + 1;
			pt3 = pent_ind * ((y + 1)*nhl + x) + !pent_ind * ((y    )*nhl + x);
		}
		else if (j == 3){
			pt1 = (y    )*nhl + x + 1;
			pt2 = pent_ind * ((y + 1)*nhl + x    ) + !pent_ind * ((y    )*nhl + x);
			pt3 = pent_ind * ((y + 2)*nhl + x + 1) + !pent_ind * ((y + 1)*nhl + x);
		}
		else if (j == 4){
			pt1 = pent_ind * ((y + 1)*nhl + x)     + !pent_ind * ((y    )*nhl + x    );
			pt2 = pent_ind * ((y + 2)*nhl + x + 1) + !pent_ind * ((y + 1)*nhl + x    );
			pt3 = pent_ind * ((y + 2)*nhl + x + 2) + !pent_ind * ((y + 2)*nhl + x + 1);
		}
		else{
			pt1 = (y + 1)*nhl + x;
			pt2 = (y + 2)*nhl + x + 1;
			pt3 = (y + 2)*nhl + x + 2;
		}

		if (laststep) vdiff = 0.5 * rscale * (2.0 * Rho_s[ir] + Rho_s[pt1] + 2.0 * Rho_s[pt2] + Rho_s[pt3]) * o6 * sdiff;
		else          vdiff = 0.5 * rscale;

		if (j == 0){
			AT1 = rscale / areasTr_d[id * 6 + j  ];
			lap1 = o6*(a_s[ir]  + a_s[pt1]) - o3*a_s[pt2];
			lap2 = o6*(a_s[pt1] + a_s[pt2]) - o3*a_s[ir] ;
			lap3 = o6*(a_s[ir]  + a_s[pt2]) - o3*a_s[pt1];
		}
		AT2 = rscale / areasTr_d[id * 6 + jp1];
		lap4 = o6*(a_s[ir]  + a_s[pt2]) - o3*a_s[pt3];
		lap5 = o6*(a_s[pt2] + a_s[pt3]) - o3*a_s[ir] ;
		lap6 = o6*(a_s[ir]  + a_s[pt3]) - o3*a_s[pt2];

		if (j == 0){
			lapx1   = (- lap1 * nvecti_d[id * 6 * 3 + j   * 3 + 0]
				    +    lap2 * nvecte_d[id * 6 * 3 + j   * 3 + 0]
				    +	 lap3 * nvecti_d[id * 6 * 3 + jp1 * 3 + 0]) * AT1;

			lapx2   = (- lap4 * nvecti_d[id * 6 * 3 + jp1 * 3 + 0]
					+	 lap5 * nvecte_d[id * 6 * 3 + jp1 * 3 + 0]
					+	 lap6 * nvecti_d[id * 6 * 3 + jp2 * 3 + 0]) * AT2;

			lapy1   = (- lap1 * nvecti_d[id * 6 * 3 + j   * 3 + 1]
				    +    lap2 * nvecte_d[id * 6 * 3 + j   * 3 + 1]
				    +	 lap3 * nvecti_d[id * 6 * 3 + jp1 * 3 + 1]) * AT1;

			lapy2   = (- lap4 * nvecti_d[id * 6 * 3 + jp1 * 3 + 1]
					+	 lap5 * nvecte_d[id * 6 * 3 + jp1 * 3 + 1]
					+	 lap6 * nvecti_d[id * 6 * 3 + jp2 * 3 + 1]) * AT2;

			lapz1   = (- lap1 * nvecti_d[id * 6 * 3 + j   * 3 + 2]
				    +    lap2 * nvecte_d[id * 6 * 3 + j   * 3 + 2]
				    +	 lap3 * nvecti_d[id * 6 * 3 + jp1 * 3 + 2]) * AT1;

			lapz2   = (- lap4 * nvecti_d[id * 6 * 3 + jp1 * 3 + 2]
					+	 lap5 * nvecte_d[id * 6 * 3 + jp1 * 3 + 2]
					+	 lap6 * nvecti_d[id * 6 * 3 + jp2 * 3 + 2]) * AT2;
		}
		else{
			lapx1   = lapx2;

			lapx2   = (- lap4 * nvecti_d[id * 6 * 3 + jp1 * 3 + 0]
					+	 lap5 * nvecte_d[id * 6 * 3 + jp1 * 3 + 0]
					+	 lap6 * nvecti_d[id * 6 * 3 + jp2 * 3 + 0]) * AT2;

			lapy1   = lapy2;

			lapy2   = (- lap4 * nvecti_d[id * 6 * 3 + jp1 * 3 + 1]
					+	 lap5 * nvecte_d[id * 6 * 3 + jp1 * 3 + 1]
					+	 lap6 * nvecti_d[id * 6 * 3 + jp2 * 3 + 1]) * AT2;

			lapz1   = lapz2;

			lapz2   = (- lap4 * nvecti_d[id * 6 * 3 + jp1 * 3 + 2]
					+	 lap5 * nvecte_d[id * 6 * 3 + jp1 * 3 + 2]
					+	 lap6 * nvecti_d[id * 6 * 3 + jp2 * 3 + 2]) * AT2;
		}

		lap += ((lapx1 + lapx2) * nvecoa_d[id * 6 * 3 + j * 3 + 0] + (lapy1 + lapy2) * nvecoa_d[id * 6 * 3 + j * 3 + 1] + (lapz1 + lapz2) * nvecoa_d[id * 6 * 3 + j * 3 + 2])*vdiff;

		if (pent_ind && j == 4) break;
	}

	if (laststep)difftr_d[id*nv*ntr + lev*ntr + var] = lap;
	else diff_d[id * nv * 6 + lev * 6 + var] = lap;
}


template <int NN>
__global__ void Tracer_Eq_Diffusion_Poles (double * difftr_d     ,
									 	   double * diff_d       ,
									 	   double * tracer_d     ,
									 	   double * Rho_d        ,
									 	   double * areasTr_d    ,
									 	   double * nvecoa_d     ,
									 	   double * nvecti_d     ,
									 	   double * nvecte_d     ,
									 	   double * K_d          ,
									 	   double * Altitude_d   ,
									 	   double * Altitudeh_d  ,
									 	   double   A            ,
									 	   int    * local_d      ,
									 	   int      ntr          ,
									 	   int      num          ,
									 	   int      laststep     ,
									 	   int      DeepModel    ){

	int id = blockIdx.x * blockDim.x + threadIdx.x;
	id += num - 2; 	// Poles

	int nv = gridDim.y;
	int lev = blockIdx.y;
	int var = blockIdx.z;

	double alt, sdiff, vdiff;
	double rscale;
	int jp1, jp2;
	int kp1, kp2;
	double o3 = 1.0 / 3.0;
	double o6 = 1.0 / 6.0;
	double AT1, AT2;
	double lap, lap1, lap2, lap3, lap4, lap5, lap6;
	double lapx1, lapx2,
		   lapy1, lapy2,
		   lapz1, lapz2;
/////////////////////////////////////////
	__shared__ double a_p[NN+1];
	__shared__ double Rho_p[NN+1];

	__shared__ double areasTr_p[NN];
	__shared__ double nvecoa_p[NN * 3];
	__shared__ double nvecti_p[NN * 3];
	__shared__ double nvecte_p[NN * 3];
	__shared__ int local_p[NN];
/////////////////////////////////////////
	
	for (int i = 0; i < 5; i++) local_p[i] = local_d[id * 6 + i];
	for (int i = 0; i < 5; i++) areasTr_p[i]  = areasTr_d[id * 6 + i];
	for (int i = 0; i < 5; i++) for (int k = 0; k < 3; k++) nvecoa_p[i * 3 + k] = nvecoa_d[id * 6 * 3 + i * 3 + k];
	for (int i = 0; i < 5; i++) for (int k = 0; k < 3; k++) nvecti_p[i * 3 + k] = nvecti_d[id * 6 * 3 + i * 3 + k];
	for (int i = 0; i < 5; i++) for (int k = 0; k < 3; k++) nvecte_p[i * 3 + k] = nvecte_d[id * 6 * 3 + i * 3 + k];


	alt = Altitude_d[lev];

	Rho_p[0] = Rho_d[id * nv + lev];
	for (int i = 1; i < 6; i++) Rho_p[i] = Rho_d[local_p[i - 1] * nv + lev];

	if (laststep) sdiff = -K_d[lev];

	if (laststep){
		a_p[0] = diff_d[id * nv *6 + lev * 6 + var];
		for (int i = 1; i < 6; i++) a_p[i] = diff_d[local_p[i - 1] * nv * 6 + lev * 6 + var];
	}
	else{
		a_p[0] = tracer_d[id*nv*ntr + lev*ntr + 0] / Rho_p[0];
		for (int i = 1; i < 6; i++) a_p[i] = tracer_d[local_p[i - 1]*nv*ntr + lev*ntr + var] / Rho_p[i];		

		if (DeepModel == 1)rscale = A / (alt + A);
		else  rscale = 1.0;

		lap = 0.0;
			
		for (int k = 0; k < 5; k++){
			int j = k + 1;
			jp1 = (j + 1) % 5;
			jp2 = (j + 2) % 5;
			kp1 = (k + 1) % 5;
			kp2 = (k + 2) % 5;

			if (laststep) vdiff = 0.5 * sdiff * rscale * (2.0 * Rho_p[0] + Rho_p[j] + 2.0 * Rho_p[jp1] + Rho_p[jp2])*o6;
			else          vdiff = 0.5 * rscale;

			if (k == 0){
				lap1 = (o6*(a_p[0] + a_p[j])   - o3*a_p[jp1]);
				lap2 = (o6*(a_p[j] + a_p[jp1]) - o3*a_p[0])  ;
				lap3 = (o6*(a_p[0] + a_p[jp1]) - o3*a_p[j])  ;
				AT1 = rscale / areasTr_p[k];
			}
			AT2 = rscale / areasTr_p[kp1];
			lap4 = (o6*(a_p[0] + a_p[jp1])   - o3*a_p[jp2]);
			lap5 = (o6*(a_p[jp1] + a_p[jp2]) - o3*a_p[0])  ;
			lap6 = (o6*(a_p[0] + a_p[jp2])   - o3*a_p[jp1]);
			
			if (k == 0){
				lapx1  = (- lap1 * nvecti_p[k * 3 + 0]
					   +    lap2 * nvecte_p[k * 3 + 0]
					   +    lap3 * nvecti_p[kp1 * 3 + 0]) * AT1;
				lapx2  = (- lap4 * nvecti_p[kp1 * 3 + 0]
					   +    lap5 * nvecte_p[kp1 * 3 + 0]
					   +    lap6 * nvecti_p[kp2 * 3 + 0]) * AT2;
				
				lapy1  = (- lap1 * nvecti_p[k * 3 + 1]
					   +	lap2 * nvecte_p[k * 3 + 1]
					   +	lap3 * nvecti_p[kp1 * 3 + 1]) * AT1;
				lapy2  = (- lap4 * nvecti_p[kp1 * 3 + 1]
					  +	    lap5 * nvecte_p[kp1 * 3 + 1]
					  +		lap6 * nvecti_p[kp2 * 3 + 1]) * AT2;
				
				lapz1  = (- lap1 * nvecti_p[k * 3 + 2]
					  +     lap2 * nvecte_p[k * 3 + 2]
					  + 	lap3 * nvecti_p[kp1 * 3 + 2]) * AT1;
				lapz2  = (- lap4 * nvecti_p[kp1 * 3 + 2]
					  +	    lap5 * nvecte_p[kp1 * 3 + 2]
					  +		lap6 * nvecti_p[kp2 * 3 + 2]) * AT2;
			}
			else{
				lapx1  = lapx2;
				lapx2  = (- lap4 * nvecti_p[kp1 * 3 + 0]
					   +    lap5 * nvecte_p[kp1 * 3 + 0]
					   +    lap6 * nvecti_p[kp2 * 3 + 0]) * AT2;
				
				lapy1  = lapy2;
				lapy2  = (- lap4 * nvecti_p[kp1 * 3 + 1]
					  +	    lap5 * nvecte_p[kp1 * 3 + 1]
					  +		lap6 * nvecti_p[kp2 * 3 + 1]) * AT2;
				
				lapz1  = lapz2;
				lapz2  = (- lap4 * nvecti_p[kp1 * 3 + 2]
					  +	    lap5 * nvecte_p[kp1 * 3 + 2]
					  +		lap6 * nvecti_p[kp2 * 3 + 2]) * AT2;

			}
			lap += ((lapx1 + lapx2) * nvecoa_p[k * 3 + 2] + (lapy1 + lapy2) * nvecoa_p[k * 3 + 1] + (lapz1 + lapz2) * nvecoa_p[k * 3 + 2])*vdiff;
		}
		if (laststep)  difftr_d[id*nv*ntr + lev*ntr + var] =  lap;
		else diff_d[id * nv * 6 + lev * 6 + var] = lap;
	}
}
