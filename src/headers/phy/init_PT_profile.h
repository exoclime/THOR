#pragma once
#include <math.h>

#include <iostream>
#include <fstream>
#include <string>

void text_file_to_array(std::string name ,double *array, int Nlen){
        
        std::ifstream inFile;
        inFile.open(name);
        if (!inFile)
        {
            printf("\nError opening the file: %s \n", name);
            
        }
        for (int i = 0; i < Nlen; i++)
        {
            inFile >> array[i];
        }  
        inFile.close();
    
}

void bilinear_interp(double xval, double yval, double x1, double x2, double y1, double y2,
    double a11, double a12, double a21, double a22, double &aval) {
    
    double lxval, lyval, lx1, lx2, ly1, ly2, la11, la21, la12, la22;
    double norm;

    lxval = xval;
    lyval = yval;
    lx1 = x1;
    lx2 = x2;
    ly1 = y1;
    ly2 = y2;
    la11 = a11;
    la21 = a21;
    la12 = a12;
    la22 = a22;

    norm = 1.0 / (lx2 - lx1) / (ly2 - ly1);

    aval = la11 * (lx2 - lxval) * (ly2 - lyval) * norm +
        la21 * (lxval - lx1) * (ly2 - lyval) * norm +
        la12 * (lx2 - lxval) * (lyval - ly1) * norm + 
        la22 * (lxval - lx1) * (lyval - ly1) * norm;

    aval = aval;

}

void bilinear_log_interp(double xval, double yval, double x1, double x2, double y1, double y2,
    double a11, double a12, double a21, double a22, double &aval) {
    
    double lxval, lyval, lx1, lx2, ly1, ly2, la11, la21, la12, la22;
    double norm;

    lxval = log10(xval);
    lyval = log10(yval);
    lx1 = log10(x1);
    lx2 = log10(x2);
    ly1 = log10(y1);
    ly2 = log10(y2);
    la11 = log10(a11);
    la21 = log10(a21);
    la12 = log10(a12);
    la22 = log10(a22);

    if (isnan( lx1))
        {                
            printf("variable lx1 in bilinear_log_interp is NaN \n");                  
        }

    if (x2 > x1)
        {                
            printf("x2 > x1 \n");                  
        }

    if (y2 > y1)
        {                
            printf("y2 > y1 \n");                  
        }

    if (x2 < x1)
        {                
            printf("x2 < x1 \n");                  
        }

    if (y2 < y1)
        {                
            printf("y2 < y1 \n");                  
        }

    if (x2 == x1)
        {                
            printf("x2 == x1 \n");                  
        }

    if (y2 == y1)
        {                
            printf("y2 == y1 \n");                  
        }

    norm = 1.0 / (lx2 - lx1) / (ly2 - ly1);

    if (isnan( la11 * (lx2 - lxval) * (ly2 - lyval) * norm))
        {                
            printf("operations (la11 * (lx2 - lxval) * (ly2 - lyval) * norm) inside bilinear_log_interp for kRoss gets NaN \n");                  
        }

    

    aval = la11 * (lx2 - lxval) * (ly2 - lyval) * norm +
        la21 * (lxval - lx1) * (ly2 - lyval) * norm +
        la12 * (lx2 - lxval) * (lyval - ly1) * norm + 
        la22 * (lxval - lx1) * (lyval - ly1) * norm;

    aval = pow(10, aval);

}


void bilinear_interpolation_polynomial_fit(double x, double y, double x1, double x2, double y1, double y2,
    double z11, double z12, double z21, double z22, double &z){
        
        double a00, a10, a01, a11;

        a00 = 1.0/((x2-x1)*(y2-y1)) * (x2*y2*z11 - x2*y1*z12 - x1*y2*z21 + x1*y1*z22);
        a10 = 1.0/((x2-x1)*(y2-y1)) * (-y2*z11 + y1*z12 + y2*z21 - y1*z22);
        a01 = 1.0/((x2-x1)*(y2-y1)) * (-x2*z11 + x2*z12 + x1*z21 - x1*z22);
        a11 = 1.0/((x2-x1)*(y2-y1)) * (z11 - z12 - z21 + z22);

        z = a00 + a10*x + a01*y + a11*x*y;

        if (isnan( (x2*y2*z11 - x2*y1*z12 - x1*y2*z21 + x1*y1*z22)))
        {                
            printf("operations (x2*y2*z11 - x2*y1*z12 - x1*y2*z21 + x1*y1*z22) inside bilinear interpolation for kRoss gets NaN \n");                  
        }
        if (isnan( (1/((x2-x1)*(y2-y1)) * (x2*y2*z11 - x2*y1*z12 - x1*y2*z21 + x1*y1*z22))))
        {                
            printf("operations (1/((x2-x1)*(y2-y1)) * (x2*y2*z11 - x2*y1*z12 - x1*y2*z21 + x1*y1*z22)) inside bilinear interpolation for kRoss greater than 1e50 \n");                  
        }

        if ((1/((x2-x1)*(y2-y1)) * (x2*y2*z11 - x2*y1*z12 - x1*y2*z21 + x1*y1*z22)) < 1e-20  || (1/((x2-x1)*(y2-y1)) * (x2*y2*z11 - x2*y1*z12 - x1*y2*z21 + x1*y1*z22)) > -1e-20 )
        {                
            printf("operations (1/((x2-x1)*(y2-y1)) * (x2*y2*z11 - x2*y1*z12 - x1*y2*z21 + x1*y1*z22)) inside bilinear interpolation for kRoss close to 0 within 1e-20 \n");                  
        }

        if ( (x2*y2*z11 - x2*y1*z12 - x1*y2*z21 + x1*y1*z22) < 1e-20 && (x2*y2*z11 - x2*y1*z12 - x1*y2*z21 + x1*y1*z22) > -1e-20 )
        {                
            printf("operations (x2*y2*z11 - x2*y1*z12 - x1*y2*z21 + x1*y1*z22) inside bilinear interpolation for kRoss close to 0 within 1e-20 \n");                  
        }

        if ( (x2*y2*z11 - x2*y1*z12 - x1*y2*z21 + x1*y1*z22) < -1e50 || (x2*y2*z11 - x2*y1*z12 - x1*y2*z21 + x1*y1*z22) > 1e50 )
        {                
            printf("operations (x2*y2*z11 - x2*y1*z12 - x1*y2*z21 + x1*y1*z22) inside bilinear interpolation for kRoss close greater than 1e50 \n");                  
        }

        if ( 1/((x2-x1)*(y2-y1)) < 1e-20 && 1/((x2-x1)*(y2-y1)) > -1e-20)
        {                
            printf("operations 1/((x2-x1)*(y2-y1)) inside bilinear interpolation for kRoss greater close to 0 within 1e-20\n");                  
        }
        if ( 1/((x2-x1)*(y2-y1)) < -1e50 || 1/((x2-x1)*(y2-y1)) > 1e50 )
        {                
            printf("operations 1/((x2-x1)*(y2-y1)) inside bilinear interpolation for kRoss greater than 1e50\n");                  
        }

        if (  isnan( 1/((x2-x1)*(y2-y1)))  )
        {                
            printf("operation 1/((x2-x1)*(y2-y1)) inside bilinear interpolation for kRoss  is NaN \n");                  
        }

        if (x1>1e10)
        {                
            printf("variable x1 greater than 1e10 \n");                  
        }

        if (x2>1e10)
        {                
            printf("variable x2 greater than 1e10 \n");                  
        }

        if (y1>1e10)
        {                
            printf("variable y1 greater than 1e10 \n");                  
        }

        if (y2>1e10)
        {                
            printf("variable y2 greater than 1e10 \n");                  
        }

        if (x1 < 1e-6)
        {                
            printf("variable x1 smaller than < 1e-6 \n");                  
        }

        if (x2 < 1e-6)
        {                
            printf("variable x2 smaller than < 1e-6 \n");                  
        }

        if (y1 < 1e-6)
        {                
            printf("variable y1 smaller than < 1e-6 \n");                  
        }

        if (y2 < 1e-6)
        {                
            printf("variable y2 smaller than < 1e-6 \n");                  
        }
        
        
        if (isnan( a00))
        {                
            printf("variable a00 inside bilinear interpolation for kRoss  is NaN \n");                  
        }
        if (isnan( a10))
        {                
            printf("variable a10 inside bilinear interpolation for kRoss  is NaN \n");                  
        }
        if (isnan( a01))
        {                
            printf("variable a01 inside bilinear interpolation for kRoss  is NaN \n");                  
        }
        if (isnan( a11))
        {                
            printf("variable a11 inside bilinear interpolation for kRoss  is NaN \n");                  
        }
        if (isnan( z))
        {                
            printf("variable z inside bilinear interpolation for kRoss  is NaN \n");                  
        }
    
}


// Calculates the IR band Rosseland mean opacity (local T) according to a
 // bilinear_interpolation_polynomial_fit to the opacities of Freedman et al. (2014)

void k_Ross_Freedman_bilinear_interpolation_polynomial_fit(double Tin, double Pin, double *OpaTableTemperature,
    double *OpaTablePressure, double *OpaTableKappa, double &k_IR){
        
    double x, y, x1, x2, y1, y2, z11, z12, z21, z22;
 
    int iter = 0;
    int jump_for_higher_temp = 15;
    double increasing_factor = 1; //1e+16;
    double dyncm_2_to_Pa = 0.1;
    int len = sizeof(OpaTableTemperature)/sizeof(*OpaTableTemperature);
    
    // exclude values off the table and insure that values within the table are used
    if (Tin <= OpaTableTemperature[0]) {
        x = OpaTableTemperature[0] + 0.11 ;      
    } else if (Tin >= OpaTableTemperature[len - 1]) {
        x = OpaTableTemperature[len - 1] - 0.1 ;      
    } else {
        x = Tin;
    }
    
    if (Pin <= OpaTablePressure[0] * dyncm_2_to_Pa) {
        y = OpaTablePressure[0] / 10 + 0.1 ;      
    } else if (Pin >= OpaTablePressure[len - 1]  * dyncm_2_to_Pa) {
        y = OpaTablePressure[len - 1]  * dyncm_2_to_Pa - 0.1 ;      
    } else {
        y = Pin;
    }
    
    // iterating through the table and get iter-position for z22
    while (x > OpaTableTemperature[iter]) {
        iter++;
    }
    
    while (y > OpaTablePressure[iter]  * dyncm_2_to_Pa) {
        iter++;
    }

    // setting all inputs variables for the interpolation
    // increased opacities for operations by 
    
    z11 = OpaTableKappa[iter - (jump_for_higher_temp) - 1] * increasing_factor;
    z12 = OpaTableKappa[iter - (jump_for_higher_temp)] * increasing_factor;
    z21 = OpaTableKappa[iter - 1] * increasing_factor;
    z22 = OpaTableKappa[iter] * increasing_factor;

    if (iter%jump_for_higher_temp==0) {
        printf("case A \n");
        x1 = OpaTableTemperature[iter - jump_for_higher_temp];
        x2 = OpaTableTemperature[iter + 1];
        y1 = OpaTablePressure[iter];
        y2 = OpaTablePressure[iter + 1] ;

        z11 = OpaTableKappa[iter - (jump_for_higher_temp)] * increasing_factor;
        z12 = OpaTableKappa[iter - (jump_for_higher_temp) + 1] * increasing_factor;
        z21 = OpaTableKappa[iter] * increasing_factor;
        z22 = OpaTableKappa[iter + 1] * increasing_factor;
    } else{
        printf("case B \n");
        x1 = OpaTableTemperature[iter - jump_for_higher_temp - 1];
        x2 = OpaTableTemperature[iter];
        y1 = OpaTablePressure[iter-1];
        y2 = OpaTablePressure[iter];

        z11 = OpaTableKappa[iter - (jump_for_higher_temp) - 1] * increasing_factor;
        z12 = OpaTableKappa[iter - (jump_for_higher_temp)] * increasing_factor;
        z21 = OpaTableKappa[iter - 1] * increasing_factor;
        z22 = OpaTableKappa[iter] * increasing_factor;    
    }

    // interpolate values from the table
    bilinear_log_interp(x, y, x1, x2, y1, y2, z11,  z12,  z21,  z22, k_IR);

    if (isnan( k_IR))
    {                
        printf("output variable k_IR for bilinear interpolation kRoss  is NaN at level   \n");                  
    }
    if (isnan( len))
    {                
        printf("input length for arrays for bilinear interpolation kRoss  is NaN at level   \n");                  
    }
    if (isnan( x))
    {                
        printf("input variable x for bilinear interpolation kRoss  is NaN at level   \n");                  
    }
    if (isnan( y))
    {                
        printf("input variable y for bilinear interpolation kRoss  is NaN at level   \n");                  
    }

    if (isnan( x1))
    {                
        printf("input variable x1 for bilinear interpolation kRoss  is NaN at level   \n");                  
    }
    if (isnan( y1))
    {                
        printf("input variable y1 for bilinear interpolation kRoss  is NaN at level   \n");                  
    }

    if (isnan( x2))
    {                
        printf("input variable x2 for bilinear interpolation kRoss  is NaN at level   \n");                  
    }
    if (isnan( y2))
    {                
        printf("input variable y2 for bilinear interpolation kRoss  is NaN at level   \n");                  
    }

    if (isnan( z11))
    {                
        printf("input variable z11 for bilinear interpolation kRoss  is NaN at level   \n");                  
    }
    if (isnan( z12))
    {                
        printf("input variable z12 for bilinear interpolation kRoss  is NaN at level   \n");                  
    }

    if (isnan( z21))
    {                
        printf("input variable z21 for bilinear interpolation kRoss  is NaN at level   \n");                  
    }
    if (isnan( z22))
    {                
        printf("input variable z22 for bilinear interpolation kRoss  is NaN at level   \n");                  
    }

    if (k_IR == 0)
    {                
        printf("output variable k_IR for bilinear interpolation is 0 at level   \n");                  
    }

    if (k_IR > 1.0)
    {                
        printf("output variable k_IR for bilinear interpolation greater than 1 at level   \n");                  
    }


    //  converted from [cm2 g-1] to [m2 kg-1] and redo temporal 1e+5 higher values
    k_IR = 10 * k_IR / increasing_factor;
    
}

void create_pressure_layers(int i, int nlay, double *(&pl), double P_Ref){

    double pe[nlay];

    if (nlay>52)
    {
        
        printf("error too many vertical layers - max. 52 layers");
    }

    double a[] = {0.05,
        0.24134615384615382,
        0.43269230769230765,
        0.6240384615384615,
        0.8153846153846154,
        1.0067307692307692,
        1.198076923076923,
        1.3894230769230769,
        1.5807692307692307,
        1.7721153846153845,
        1.9634615384615384,
        2.154807692307692,
        2.346153846153846,
        2.5374999999999996,
        2.7288461538461535,
        2.9201923076923073,
        3.111538461538461,
        3.302884615384615,
        3.494230769230769,
        3.6855769230769226,
        3.8769230769230765,
        4.06826923076923,
        4.259615384615384,
        4.450961538461538,
        4.642307692307692,
        4.833653846153846,
        5.0249999999999995,
        5.216346153846153,
        5.407692307692307,
        5.599038461538461,
        5.790384615384615,
        5.981730769230769,
        6.1730769230769225,
        6.364423076923076,
        6.55576923076923,
        6.747115384615384,
        6.938461538461538,
        7.129807692307692,
        7.3211538461538455,
        7.512499999999999,
        7.703846153846153,
        7.895192307692307,
        8.086538461538462,
        8.277884615384615,
        8.46923076923077,
        8.660576923076924,
        8.851923076923077,
        9.04326923076923,
        9.234615384615385,
        9.42596153846154,
        9.617307692307692,
        9.808653846153845,
        10.0};

    double b[] = {2.2727272727272727e-09,
        3.3324801199555632e-09,
        4.886386449955581e-09,
        7.164865709274777e-09,
        1.0505779916856185e-08,
        1.5404533195722723e-08,
        2.2587532278054278e-08,
        3.311990100120702e-08,
        4.8563421131024246e-08,
        7.120811960951402e-08,
        1.044118429103743e-07,
        1.5309817194616815e-07,
        2.244865103413433e-07,
        3.2916260648071946e-07,
        4.826482506250917e-07,
        7.077029080613574e-07,
        1.0376985836576538e-06,
        1.5215683562398778e-06,
        2.231062371262057e-06,
        3.2713872393891584e-06,
        4.7968064935739005e-06,
        7.033515402807869e-06,
        1.0313182111433727e-05,
        1.5122128718327065e-05,
        2.2173445063102982e-05,
        3.251272854003574e-05,
        4.767312946227964e-05,
        6.99026927232105e-05,
        0.00010249770688584301,
        0.00015029149103676587,
        0.00022037109867454092,
        0.0003231282143520963,
        0.0004738000742310444,
        0.0006947289044117072,
        0.0010186749155926304,
        0.001493674118160335,
        0.0021901612939630483,
        0.0032114143475163253,
        0.004708868766816817,
        0.006904573083274554,
        0.010124115116188805,
        0.014844901437009984,
        0.02176694912547564,
        0.031916687102401954,
        0.04679915911598253,
        0.06862119764924554,
        0.10061866186840697,
        0.14753626376411874,
        0.21633113302523627,
        0.31720444805898806,
        0.4651141075319441,
        0.6819927474188744,
        1.0000000000000009};

    // Contruct pressure array in pa
    for (int lev = 0; lev < (nlay+1); lev++)
    {
        pe[lev] = a[nlay+1-lev] + b[nlay+1-lev]*P_Ref;
        
    }

    //! Pressure layers

    for (int lev = 0; lev < nlay; lev++)
    {
        pl[lev] = (pe[ lev] - pe[lev + 1] ) / ( logl(pe[lev] ) - logl( pe[lev+1] ));
        
    }

}

// Calculates the IR band Rosseland mean opacity (local T) according to the
    // Freedman et al. (2014) fit and coefficents

void k_Ross_Freedman(double Tin, double Pin, double met, double& k_IR) {
    // dependcies
    //// pow from math
    //// log10 from math        
    //// atan from math
    //// onedivpi -> namespace constants::onedivpi

    // Input:
    // T - Local gas temperature [K]
    // P - Local gas pressure [pa]
    // met - Local metallicity [M/H] (log10 from solar, solar [M/H] = 0.0)

    // Call by reference (Input&Output):
    // k_IR - IR band Rosseland mean opacity [m2 kg-1]

    // work variables
    double k_lowP;
    double k_hiP;
    double T;
    double P;
    double Tl10;
    double Pl10;
    const double pi = atan((double)(1)) * 4;
    const double onedivpi = 1.0 / pi;

    // Coefficent parameters for Freedman et al. (2014) table fit
    double c1 = 10.602;
    double c2 = 2.882;
    double c3 = 6.09e-15;
    double c4 = 2.954;
    double c5 = -2.526;
    double c6 = 0.843;
    double c7 = -5.490;
    double c8_l = -14.051, c8_h = 82.241;
    double c9_l = 3.055, c9_h = -55.456;
    double c10_l = 0.024, c10_h = 8.754;
    double c11_l = 1.877, c11_h = 0.7048;
    double c12_l = -0.445, c12_h = -0.0414;
    double c13_l = 0.8321, c13_h = 0.8321;

    // start operations

    T = Tin;
    P = Pin * ((double)10.0); // Convert to dyne cm-2


    Tl10 = log10(T);
    Pl10 = log10(P);

    // Low pressure expression
    k_lowP = c1 * atan(Tl10 - c2) -
        (c3 / (Pl10 + c4)) * exp(pow((Tl10 - c5), 2.0)) + c6 * met + c7;

    // De log10
    k_lowP = pow(((double)10.0), k_lowP);

    // Temperature split for coefficents = 800 K
    if (T <= 800.0)
    {
        k_hiP = c8_l + c9_l * Tl10 + c10_l * pow(Tl10, 2.0) +
            Pl10 * (c11_l + c12_l * Tl10) +
            c13_l * met * (0.5 + onedivpi * atan((Tl10 - ((double)2.5)) / (double)0.2));
    }
    else
    {
        k_hiP = c8_h + c9_h * Tl10 +
            c10_h * pow(Tl10, 2.0) + Pl10 * (c11_h + c12_h * Tl10) +
            c13_h * met * (0.5 + onedivpi * atan((Tl10 - ((double)2.5)) / (double)0.2));
    }

    // De log10
    k_hiP = pow(((double)10.0), k_hiP);

    // Total Rosseland mean opacity - converted to m2 kg-1
    k_IR = (k_lowP + k_hiP) / ((double)10.0);   

    // Avoid divergence in fit for large values
    if (k_IR > 1.0e10)
    {
        k_IR = 1.0e10;
    }
}

///////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

// Calculates 3 band grey visual gamma values and 2 picket fence IR gamma values
    // according to the coefficents and equations in:
    // Parmentier & Menou (2014) and Parmentier et al. (2015)
    // NOTE: This does not calculate the opacity - call k_Ross_Freedman for that
void gam_Parmentier(double Teff, int table_num, double(&gam_V)[3], double(&Beta_V)[3],
    double(&Beta)[2], double& gam_1, double& gam_2, double& gam_P, double& tau_lim) {
    // dependcies
    //// pow from math
    //// log10 from math        


    // Input:
    // Teff - Effective temperature [K] (See Parmentier papers for various ways to calculate this)
    // for non-irradiated atmosphere Teff = Tint
    // table_num - Table selection from Parmentier et al. (2015): 1 = w. TiO/VO, 2 = w.o. TiO/VO

    // Call by reference (Input&Output):
    // gam_V(3) - gamma ratio for 3 visual bands (gam_V = kV_Ross/kIR_Ross)
    // beta_V(3) - fraction of total incident stellar flux in band (1/3 for Parmentier values)
    // Beta - equilvalent bandwidth for picket fence IR model
    // gam_1 - gamma ratio for IR band 1 (gam_1 = kIR_1/kIR_Ross)
    // gam_2 - gamma ratio for IR band 2 (gam_2 = kIR_2/kIR_Ross)
    // gam_P - gamma ratio for Planck mean (gam_P = kIR_Planck/kIR_Ross)
    // tau_lim - tau limit variable (usually for IC system)

    // work variables
    double  R = 0;
    double aP = 0;
    double bP = 0;
    double cP = 0;
    double aV1 = 0, bV1 = 0, aV2 = 0, bV2 = 0, aV3 = 0, bV3 = 0;
    double aB = 0, bB = 0;
    double l10T = 0, l10T2 = 0, RT = 0;
    int i;

    // start operations

    // Log 10 T_eff variables
    l10T = log10(Teff);
    l10T2 = pow(l10T, 2.0);

    if (table_num == 1) {
        // First table in Parmentier et al. (2015) w. TiO/VO
        // Start large if statements with visual band and Beta coefficents
        if (Teff <= 200.0)
        {
            aV1 = -5.51; bV1 = 2.48;
            aV2 = -7.37; bV2 = 2.53;
            aV3 = -3.03; bV3 = -0.20;
            aB = 0.84; bB = 0.0;
        }
        else if (Teff > 200.0 && Teff <= 300.0)
        {
            aV1 = 1.23; bV1 = -0.45;
            aV2 = 13.99; bV2 = -6.75;
            aV3 = -13.87; bV3 = 4.51;
            aB = 0.84; bB = 0.0;
        }
        else if (Teff > 300.0 && Teff <= 600.0)
        {
            aV1 = 8.65; bV1 = -3.45;
            aV2 = -15.18; bV2 = 5.02;
            aV3 = -11.95; bV3 = 3.74;
            aB = 0.84; bB = 0.0;
        }
        else if (Teff > 600.0 && Teff <= 1400.0)
        {
            aV1 = -12.96; bV1 = 4.33;
            aV2 = -10.41; bV2 = 3.31;
            aV3 = -6.97; bV3 = 1.94;
            aB = 0.84; bB = 0.0;
        }
        else if (Teff > 1400.0 && Teff < 2000.0)
        {
            aV1 = -23.75; bV1 = 7.76;
            aV2 = -19.95; bV2 = 6.34;
            aV3 = -3.65; bV3 = 0.89;
            aB = 0.84; bB = 0.0;
        }
        else if (Teff >= 2000.0)
        {
            aV1 = 12.65; bV1 = -3.27;
            aV2 = 13.56; bV2 = -3.81;
            aV3 = -6.02; bV3 = 1.61;
            aB = 6.21; bB = -1.63;
        }

        // gam_P coefficents
        aP = -2.36;
        bP = 13.92;
        cP = -19.38;
    }
    else if (table_num == 2)
    {
        // ! Appendix table from Parmentier et al. (2015) - without TiO and VO
        if (Teff <= 200.0)
        {
            aV1 = -5.51; bV1 = 2.48;
            aV2 = -7.37; bV2 = 2.53;
            aV3 = -3.03; bV3 = -0.20;
            aB = 0.84; bB = 0.0;
        }
        else if (Teff > 200.0 && Teff <= 300.0)
        {
            aV1 = 1.23; bV1 = -0.45;
            aV2 = 13.99; bV2 = -6.75;
            aV3 = -13.87; bV3 = 4.51;
            aB = 0.84; bB = 0.0;
        }
        else if (Teff > 300.0 && Teff <= 600.0)
        {
            aV1 = 8.65; bV1 = -3.45;
            aV2 = -15.18; bV2 = 5.02;
            aV3 = -11.95; bV3 = 3.74;
            aB = 0.84; bB = 0.0;
        }
        else if (Teff > 600.0 && Teff <= 1400.0)
        {
            aV1 = -12.96; bV1 = 4.33;
            aV2 = -10.41; bV2 = 3.31;
            aV3 = -6.97; bV3 = 1.94;
            aB = 0.84; bB = 0.0;
        }
        else if (Teff > 1400.0 && Teff < 2000.0)
        {
            aV1 = -1.68; bV1 = 0.75;
            aV2 = 6.96; bV2 = -2.21;
            aV3 = 0.02; bV3 = -0.28;
            aB = 3.0; bB = -0.69;
        }
        else if (Teff >= 2000.0)
        {
            aV1 = 10.37; bV1 = -2.91;
            aV2 = -2.4; bV2 = 0.62;
            aV3 = -16.54; bV3 = 4.74;
            aB = 3.0; bB = -0.69;
        }

        // gam_P coefficents
        if (Teff <= 1400.0)
        {
            aP = -2.36;
            bP = 13.92;
            cP = -19.38;
        }
        else
        {
            aP = -12.45;
            bP = 82.25;
            cP = -134.42;
        }
    }

    // Calculation of all values
    // Visual band gamma
    gam_V[0] = pow(((double)10.0), (aV1 + bV1 * l10T));
    gam_V[1] = pow(((double)10.0), (aV2 + bV2 * l10T));
    gam_V[2] = pow(((double)10.0), (aV3 + bV3 * l10T));



    // Visual band fractions
    for (i = 0; i < 3; i++)
    {
        Beta_V[i] = ((double)1.0) / ((double)3.0);
    }

    // gamma_Planck - if < 1 then make it grey approximation (k_Planck = k_Ross, gam_P = 1)
    gam_P = pow(((double)10.0), (aP * l10T2 + bP * l10T + cP));
    if (gam_P < 1.0000001)
    {
        gam_P = 1.0000001;
    }

    // equivalent bandwidth value
    Beta[0] = aB + bB * l10T;
    Beta[1] = (1.0) - Beta[0];

    // IR band kappa1/kappa2 ratio - Eq. 96 from Parmentier & Menou (2014)
    RT = (gam_P - 1.0) / (2.0 * Beta[0] * Beta[1]);
    R = 1.0 + RT + sqrt(pow(RT, 2.0) + RT);

    // gam_1 and gam_2 values - Eq. 92, 93 from Parmentier & Menou (2014)
    gam_1 = Beta[0] + R - Beta[0] * R;
    gam_2 = gam_1 / R;

    // Calculate tau_lim parameter
    tau_lim = ((double)1.0) / (gam_1 * gam_2) * sqrt(gam_P / ((double)3.0));


}

///////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////







    // Calculates the Bond Albedo according to Parmentier et al. (2015) expression
void Bond_Parmentier_host(double Teff0, double grav, double& AB) {
    // dependcies
    //// pow from math
    //// log10 from math        


    // Input:
    // Teff0 - Atmospheric profile effective temperature [K] with zero albedo
    // grav - Surface gravity of planet [m s-2]

    // Call by reference (Input&Output):
    // AB - Bond albedo

    // work variables
    double a, b;

    // start operations

    if (Teff0 <= 250.0)
    {
        a = ((double)-0.335) * pow(grav, ((double)0.070));
        b = 0.0;
    }
    else if (Teff0 > 250.0 && Teff0 <= 750.0)
    {
        a = -0.335 * pow(grav, ((double)0.070)) + 2.149 * pow(grav, ((double)0.135));
        b = -0.896 * pow(grav, ((double)0.135));
    }
    else if (Teff0 > 750.0 && Teff0 < 1250.0)
    {
        a = -0.335 * pow(grav, ((double)0.070)) - 0.428 * pow(grav, ((double)0.135));
        b = 0.0;
    }
    else if (Teff0 >= 1250.0)
    {
        a = 16.947 - ((double)3.174) * pow(grav, ((double)0.070)) - 4.051 *
            pow(grav, ((double)0.135));
        b = -5.472 + ((double)0.917) * pow(grav, ((double)0.070)) + 1.170 *
            pow(grav, ((double)0.135));
    }

    // Final Bond Albedo expression
    AB = pow(((double)10.0), (a + b * log10(Teff0)));

}



///////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

// This subroutine follows Parmentier & Guillot (2014, 2015) non-grey picket fence scheme
void Parmentier_IC(int id, const int nlay, double* pl, double Tint, double mu, double Tirr,
    double grav, double* (&Tl), int  table_num, double met) {
    // dependcies
    //// pow -> math    
    //// sqrt -> math
    //// exp -> math

    // Input:
    // 

    // Call by reference (Input & Output):
    // 

    // work variables
    int i, j, k;
    double Teff0, Teff, Tmu, Bond, Tskin;
    double gam_V[3] = { 0 }, Beta_V[3] = { 0 };
    double Beta[2];
    double gam_1, gam_2, gam_P, tau_lim;
    double a0, a1, b0, A, B, At1, At2;
    double a2[3], a3[3], b1[3], b2[3], b3[3], Av1[3], Av2[3];
    double C[3], D[3], E[3];
    double kRoss[nlay];
    double tau[nlay + 1];

    
    
    double summy;

    // start operations

    // Effective temperature parameter
    Tmu = pow((mu * pow(Tirr, 4.0)), (1.0 / 4.0));

    // Find Bond albedo of planet - Bond albedo is given by mu = 1/sqrt(3)
    Teff0 = pow(((pow(Tint, 4.0) + (1.0 / sqrt(((double)3.0))) * pow(Tirr, 4.0))), (1.0 / 4.0));
    Bond_Parmentier_host(Teff0, grav, Bond);

    



    Teff = pow((pow(Tint, 4.0) + (((double)1.0) - Bond) * mu * pow(Tirr, 4.0)), (1.0 / 4.0));


    // Find the V band gamma, beta and IR gamma and beta ratios for this profile
    // Passed mu, so make lat = acos(mu) and lon = 0
    gam_Parmentier(Teff, table_num, gam_V, Beta_V, Beta, gam_1, gam_2,
        gam_P, tau_lim);

    

    for (i = 0; i < 3; i++)
    {
        gam_V[i] = gam_V[i] / mu;
    }

    // Hard work starts here - first calculate all the required coefficents
    At1 = pow(gam_1, 2.0) * log(1.0 + 1.0 / (tau_lim * gam_1));
    At2 = pow(gam_2, 2.0) * log(1.0 + 1.0 / (tau_lim * gam_2));
    for (i = 0; i < 3; i++)
    {
        Av1[i] = pow(gam_1, 2.0) * log(1.0 + gam_V[i] / gam_1);
        Av2[i] = pow(gam_2, 2.0) * log(1.0 + gam_V[i] / gam_2);
    }

    a0 = 1.0 / gam_1 + 1.0 / gam_2;

    a1 = -1.0 / (((double)3.0) * pow(tau_lim, 2.0)) * (gam_P / (1.0 - gam_P) *
        (gam_1 + gam_2 - 2.0) / (gam_1 + gam_2) +
        (gam_1 + gam_2) * tau_lim - (At1 + At2) * pow(tau_lim, 2.0));

    for (i = 0; i < 3; i++)
    {
        a2[i] = pow(tau_lim, 2.0) / (gam_P * pow(gam_V[i], 2.0)) *
            ((3.0 * pow(gam_1, 2.0) - pow(gam_V[i], 2.0)) * (3.0 * pow(gam_2, 2.0) - pow(gam_V[i], 2.0)) *
                (gam_1 + gam_2) - 3.0 * gam_V[i] * (6.0 * pow(gam_1, 2.0) * pow(gam_2, 2.0) - pow(gam_V[i], 2.0) *
                    (pow(gam_1, 2.0) + pow(gam_2, 2.0)))) / (1.0 - pow(gam_V[i], 2.0) * pow(tau_lim, 2.0));

        a3[i] = -pow(tau_lim, 2.0) * (3.0 * pow(gam_1, 2.0) - pow(gam_V[i], 2.0)) *
            (3.0 * pow(gam_2, 2.0) - pow(gam_V[i], 2.0)) * (Av2[i] + Av1[i]) /
            (gam_P * pow(gam_V[i], 3.0) * (1.0 - pow(gam_V[i], 2.0) * pow(tau_lim, 2.0)));

        b1[i] = gam_1 * gam_2 * (3.0 * pow(gam_1, 2.0) - pow(gam_V[i], 2.0)) * (3.0 * pow(gam_2, 2.0) -
            pow(gam_V[i], 2.0)) * pow(tau_lim, 2) / (gam_P * pow(gam_V[i], 2.0) *
                (pow(gam_V[i], 2.0) * pow(tau_lim, 2.0) - 1.0));

        b2[i] = 3.0 * (gam_1 + gam_2) * pow(gam_V[i], 3.0) / ((3.0 * pow(gam_1, 2.0) - pow(gam_V[i], 2.0)) *
            (3.0 * pow(gam_2, 2.0) - pow(gam_V[i], 2.0)));

        b3[i] = (Av2[i] - Av1[i]) / (gam_V[i] * (gam_1 - gam_2));
    }

    b0 = 1.0 / (gam_1 * gam_2 / (gam_1 - gam_2) * (At1 - At2) / 3.0 - pow((gam_1 * gam_2), 2.0) /
        sqrt(3.0 * gam_P) - pow((gam_1 * gam_2), 3.0) /
        ((1.0 - gam_1) * (1.0 - gam_2) * (gam_1 + gam_2)));

    A = 1.0 / ((double)3.0) * (a0 + a1 * b0);
    B = -1.0 / ((double)3.0) * pow((gam_1 * gam_2), 2.0) / gam_P * b0;

    for (i = 0; i < 3; i++)
    {
        C[i] = -1.0 / ((double)3.0) * (b0 * b1[i] * (1.0 + b2[i] + b3[i]) * a1 + a2[i] + a3[i]);

        D[i] = 1.0 / ((double)3.0) * pow((gam_1 * gam_2), 2.0) / gam_P * b0 * b1[i] * (1.0 + b2[i] + b3[i]);

        E[i] = (3.0 - pow((gam_V[i] / gam_1), 2.0)) * (3.0 - pow((gam_V[i] / gam_2), 2.0)) /
            (9.0 * gam_V[i] * (pow((gam_V[i] * tau_lim), 2.0) - 1.0));
    }

    
    // T-p structure calculation - we follow exactly V. Parmentier's method
    // Estimate the skin temperature by setting tau = 0
    tau[nlay-1] = 0.0;
    summy = 0.0;
    for (i = 0; i < 3; i++)
    {
        summy += 3.0 * Beta_V[i] * pow(Tmu, 4.0) / 4.0 * (C[i] + D[i] * exp(-tau[nlay-1] / tau_lim) +
            E[i] * exp(-gam_V[i] * tau[nlay-1]));
    }

    Tskin = 3.0 * pow(Tint, 4) / 4.0 * (tau[nlay-1] + A + B * exp(-tau[nlay-1] / tau_lim)) + summy;
    Tskin = pow(Tskin, (1.0 / 4.0));

    
    // Estimate the opacity TOA at the skin temperature - assume this is = first layer optacity
    k_Ross_Freedman(Tskin, pl[id * nlay + nlay-1], met, kRoss[nlay-1]);

    
    
    


    // Recalculate the upmost tau with new kappa
    tau[nlay-1] = kRoss[nlay-1] / grav * pl[id * nlay + nlay-1];

    
    
    // More accurate layer T at uppermost layer
    summy = 0.0;
    for (i = 0; i < 3; i++)
    {
        summy += 3.0 * Beta_V[i] * pow(Tmu, 4.0) / 4.0 * (C[i] + D[i] * exp(-tau[nlay-1] / tau_lim) +
            E[i] * exp(-gam_V[i] * tau[nlay-1]));
    }
    
    Tl[id * nlay + nlay-1] = 3.0 * pow(Tint, 4) / 4.0 * (tau[nlay-1] + A + B * exp(-tau[nlay-1] / tau_lim)) + summy;
    Tl[id * nlay + nlay-1] = pow(Tl[id * nlay + nlay-1], (1.0 / 4.0));

    

    
    // Now we can loop in optical depth space to find the T-p profile
    for (i = nlay-2; i>-1; i--)
    {
        // Initial guess for layer
        k_Ross_Freedman(Tl[id * nlay + i+1], sqrt(pl[id * nlay + i+1] * pl[id * nlay + i]), met, kRoss[i]);
        
        tau[i] = tau[i+1] + kRoss[i] / grav *  (pl[id * nlay + i] - pl[id * nlay + i+1]) ;
    
        summy = 0.0;
        for (j = 0; j < 3; j++)
        {
            summy = +3.0 * Beta_V[j] * pow(Tmu, 4.0) / 4.0 * (C[j] + D[j] * exp(-tau[i] / tau_lim) +
                E[j] * exp(-gam_V[j] * tau[i]));
        }
        Tl[id * nlay + i] = 3.0 * pow(Tint, 4.0) / 4.0 * (tau[i] + A + B * exp(-tau[i] / tau_lim)) + summy;
        Tl[id * nlay + i] = pow(Tl[id * nlay + i], (1.0 / 4.0));

        // Convergence loop
        for (j = 0; j < 5; j++)
        {
            k_Ross_Freedman(sqrt(Tl[id * nlay + i+1] * Tl[id * nlay + i]), sqrt(pl[id * nlay + i+1] * pl[id * nlay + i]), met, kRoss[i]);
            
            tau[i] = tau[i+1] + kRoss[i] / grav *  (pl[id * nlay + i] - pl[id * nlay + i+1]);
            summy = 0.0;
            for (k = 0; k < 3; k++)
            {
                summy += 3.0 * Beta_V[k] * pow(Tmu, 4.0) / 4.0 * (C[k] + D[k] * exp(-tau[i] / tau_lim) +
                    E[k] * exp(-gam_V[k] * tau[i]));
            }
            Tl[id * nlay + i] = 3.0 * pow(Tint, 4.0) / 4.0 * (tau[i] + A + B * exp(-tau[i] / tau_lim)) + summy;
            Tl[id * nlay + i] = pow(Tl[id * nlay + i], (1.0 / 4.0));
        }
        
    }


}

///////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

// This subroutine follows Parmentier & Guillot (2014, 2015) non-grey picket fence scheme
void Parmentier_bilinear_interpolation_IC(int id, const int nlay, double* pl, double Tint, double mu, double Tirr,
    double *OpaTableTemperature, double *OpaTablePressure, double *OpaTableKappa, double grav, double* (&Tl), int  table_num, double met) {
    // dependcies
    //// pow -> math    
    //// sqrt -> math
    //// exp -> math

    // Input:
    // 

    // Call by reference (Input & Output):
    // 

    // work variables
    int i, j, k;
    double Teff0, Teff, Tmu, Bond, Tskin;
    double gam_V[3] = { 0 }, Beta_V[3] = { 0 };
    double Beta[2];
    double gam_1, gam_2, gam_P, tau_lim;
    double a0, a1, b0, A, B, At1, At2;
    double a2[3], a3[3], b1[3], b2[3], b3[3], Av1[3], Av2[3];
    double C[3], D[3], E[3];
    double kRoss[nlay];
    double tau[nlay + 1];

    
    
    double summy;

    // start operations

    // Effective temperature parameter
    Tmu = pow((mu * pow(Tirr, 4.0)), (1.0 / 4.0));

    // Find Bond albedo of planet - Bond albedo is given by mu = 1/sqrt(3)
    Teff0 = pow(((pow(Tint, 4.0) + (1.0 / sqrt(((double)3.0))) * pow(Tirr, 4.0))), (1.0 / 4.0));
    Bond_Parmentier_host(Teff0, grav, Bond);

    



    Teff = pow((pow(Tint, 4.0) + (((double)1.0) - Bond) * mu * pow(Tirr, 4.0)), (1.0 / 4.0));


    // Find the V band gamma, beta and IR gamma and beta ratios for this profile
    // Passed mu, so make lat = acos(mu) and lon = 0
    gam_Parmentier(Teff, table_num, gam_V, Beta_V, Beta, gam_1, gam_2,
        gam_P, tau_lim);

    

    for (i = 0; i < 3; i++)
    {
        gam_V[i] = gam_V[i] / mu;
    }

    // Hard work starts here - first calculate all the required coefficents
    At1 = pow(gam_1, 2.0) * log(1.0 + 1.0 / (tau_lim * gam_1));
    At2 = pow(gam_2, 2.0) * log(1.0 + 1.0 / (tau_lim * gam_2));
    for (i = 0; i < 3; i++)
    {
        Av1[i] = pow(gam_1, 2.0) * log(1.0 + gam_V[i] / gam_1);
        Av2[i] = pow(gam_2, 2.0) * log(1.0 + gam_V[i] / gam_2);
    }

    a0 = 1.0 / gam_1 + 1.0 / gam_2;

    a1 = -1.0 / (((double)3.0) * pow(tau_lim, 2.0)) * (gam_P / (1.0 - gam_P) *
        (gam_1 + gam_2 - 2.0) / (gam_1 + gam_2) +
        (gam_1 + gam_2) * tau_lim - (At1 + At2) * pow(tau_lim, 2.0));

    for (i = 0; i < 3; i++)
    {
        a2[i] = pow(tau_lim, 2.0) / (gam_P * pow(gam_V[i], 2.0)) *
            ((3.0 * pow(gam_1, 2.0) - pow(gam_V[i], 2.0)) * (3.0 * pow(gam_2, 2.0) - pow(gam_V[i], 2.0)) *
                (gam_1 + gam_2) - 3.0 * gam_V[i] * (6.0 * pow(gam_1, 2.0) * pow(gam_2, 2.0) - pow(gam_V[i], 2.0) *
                    (pow(gam_1, 2.0) + pow(gam_2, 2.0)))) / (1.0 - pow(gam_V[i], 2.0) * pow(tau_lim, 2.0));

        a3[i] = -pow(tau_lim, 2.0) * (3.0 * pow(gam_1, 2.0) - pow(gam_V[i], 2.0)) *
            (3.0 * pow(gam_2, 2.0) - pow(gam_V[i], 2.0)) * (Av2[i] + Av1[i]) /
            (gam_P * pow(gam_V[i], 3.0) * (1.0 - pow(gam_V[i], 2.0) * pow(tau_lim, 2.0)));

        b1[i] = gam_1 * gam_2 * (3.0 * pow(gam_1, 2.0) - pow(gam_V[i], 2.0)) * (3.0 * pow(gam_2, 2.0) -
            pow(gam_V[i], 2.0)) * pow(tau_lim, 2) / (gam_P * pow(gam_V[i], 2.0) *
                (pow(gam_V[i], 2.0) * pow(tau_lim, 2.0) - 1.0));

        b2[i] = 3.0 * (gam_1 + gam_2) * pow(gam_V[i], 3.0) / ((3.0 * pow(gam_1, 2.0) - pow(gam_V[i], 2.0)) *
            (3.0 * pow(gam_2, 2.0) - pow(gam_V[i], 2.0)));

        b3[i] = (Av2[i] - Av1[i]) / (gam_V[i] * (gam_1 - gam_2));
    }

    b0 = 1.0 / (gam_1 * gam_2 / (gam_1 - gam_2) * (At1 - At2) / 3.0 - pow((gam_1 * gam_2), 2.0) /
        sqrt(3.0 * gam_P) - pow((gam_1 * gam_2), 3.0) /
        ((1.0 - gam_1) * (1.0 - gam_2) * (gam_1 + gam_2)));

    A = 1.0 / ((double)3.0) * (a0 + a1 * b0);
    B = -1.0 / ((double)3.0) * pow((gam_1 * gam_2), 2.0) / gam_P * b0;

    for (i = 0; i < 3; i++)
    {
        C[i] = -1.0 / ((double)3.0) * (b0 * b1[i] * (1.0 + b2[i] + b3[i]) * a1 + a2[i] + a3[i]);

        D[i] = 1.0 / ((double)3.0) * pow((gam_1 * gam_2), 2.0) / gam_P * b0 * b1[i] * (1.0 + b2[i] + b3[i]);

        E[i] = (3.0 - pow((gam_V[i] / gam_1), 2.0)) * (3.0 - pow((gam_V[i] / gam_2), 2.0)) /
            (9.0 * gam_V[i] * (pow((gam_V[i] * tau_lim), 2.0) - 1.0));
    }

    
    // T-p structure calculation - we follow exactly V. Parmentier's method
    // Estimate the skin temperature by setting tau = 0
    tau[nlay-1] = 0.0;
    summy = 0.0;
    for (i = 0; i < 3; i++)
    {
        summy += 3.0 * Beta_V[i] * pow(Tmu, 4.0) / 4.0 * (C[i] + D[i] * exp(-tau[nlay-1] / tau_lim) +
            E[i] * exp(-gam_V[i] * tau[nlay-1]));
    }

    Tskin = 3.0 * pow(Tint, 4) / 4.0 * (tau[nlay-1] + A + B * exp(-tau[nlay-1] / tau_lim)) + summy;
    Tskin = pow(Tskin, (1.0 / 4.0));

    
    // Estimate the opacity TOA at the skin temperature - assume this is = first layer optacity
    k_Ross_Freedman_bilinear_interpolation_polynomial_fit(Tskin, pl[id * nlay + nlay-1], OpaTableTemperature, OpaTablePressure, OpaTableKappa, kRoss[nlay-1]);
    if (isnan( kRoss[nlay-1]))
    {                
        printf("after bilinear interpolation kRoss[nlay-1] is NaN at level %d  \n", nlay-1);                  
    }
    
    
    


    // Recalculate the upmost tau with new kappa
    tau[nlay-1] = kRoss[nlay-1] / grav * pl[id * nlay + nlay-1];

    
    
    // More accurate layer T at uppermost layer
    summy = 0.0;
    for (i = 0; i < 3; i++)
    {
        summy += 3.0 * Beta_V[i] * pow(Tmu, 4.0) / 4.0 * (C[i] + D[i] * exp(-tau[nlay-1] / tau_lim) +
            E[i] * exp(-gam_V[i] * tau[nlay-1]));
    }
    
    Tl[id * nlay + nlay-1] = 3.0 * pow(Tint, 4) / 4.0 * (tau[nlay-1] + A + B * exp(-tau[nlay-1] / tau_lim)) + summy;
    Tl[id * nlay + nlay-1] = pow(Tl[id * nlay + nlay-1], (1.0 / 4.0));

    

    
    // Now we can loop in optical depth space to find the T-p profile
    for (i = nlay-2; i>-1; i--)
    {
        // Initial guess for layer
        k_Ross_Freedman(Tl[id * nlay + i+1], sqrt(pl[id * nlay + i+1] * pl[id * nlay + i]), met, kRoss[i]);
        //k_Ross_Freedman_bilinear_interpolation_polynomial_fit(Tl[id * nlay + i+1], sqrt(pl[id * nlay + i+1] * pl[id * nlay + i]), OpaTableTemperature, OpaTablePressure, OpaTableKappa, kRoss[i]);
        

        tau[i] = tau[i+1] + kRoss[i] / grav *  (pl[id * nlay + i] - pl[id * nlay + i+1]) ;
    
        summy = 0.0;
        for (j = 0; j < 3; j++)
        {
            summy = +3.0 * Beta_V[j] * pow(Tmu, 4.0) / 4.0 * (C[j] + D[j] * exp(-tau[i] / tau_lim) +
                E[j] * exp(-gam_V[j] * tau[i]));
        }
        Tl[id * nlay + i] = 3.0 * pow(Tint, 4.0) / 4.0 * (tau[i] + A + B * exp(-tau[i] / tau_lim)) + summy;
        Tl[id * nlay + i] = pow(Tl[id * nlay + i], (1.0 / 4.0));

        // Convergence loop
        for (j = 0; j < 5; j++)
        {
           k_Ross_Freedman(sqrt(Tl[id * nlay + i+1] * Tl[id * nlay + i]), sqrt(pl[id * nlay + i+1] * pl[id * nlay + i]), met, kRoss[i]);
           //k_Ross_Freedman_bilinear_interpolation_polynomial_fit(sqrt(Tl[id * nlay + i+1] * Tl[id * nlay + i]), sqrt(pl[id * nlay + i+1] * pl[id * nlay + i]), OpaTableTemperature, OpaTablePressure, OpaTableKappa, kRoss[i]);
        
            
            tau[i] = tau[i+1] + kRoss[i] / grav *  (pl[id * nlay + i] - pl[id * nlay + i+1]);
            summy = 0.0;
            for (k = 0; k < 3; k++)
            {
                summy += 3.0 * Beta_V[k] * pow(Tmu, 4.0) / 4.0 * (C[k] + D[k] * exp(-tau[i] / tau_lim) +
                    E[k] * exp(-gam_V[k] * tau[i]));
            }
            Tl[id * nlay + i] = 3.0 * pow(Tint, 4.0) / 4.0 * (tau[i] + A + B * exp(-tau[i] / tau_lim)) + summy;
            Tl[id * nlay + i] = pow(Tl[id * nlay + i], (1.0 / 4.0));
        }
        
    }


}

///////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

// Subroutine that corrects for adiabatic region following Parmentier & Guillot (2015)
void adiabat_correction(int id, int nlay, double* (&Tl), double* pressure_h, double Gravit) {
    // dependcies
    //// main_parameters::nlay  -> "FMS_RC_para_&_const.cpp" 
    //// pow -> math    
    /// log10 -> math

    // Input:
    // 

    // Call by reference (Input & Output):
    // 

    // work variables
    int i, iRC, iRC1;
    double gradrad[nlay-1];
    double gradad[nlay-1];




    // start operations

    for (i = (nlay-1); i > 1; i--)
    {
        
        gradrad[i] = (log10(Tl[i]) - log10(Tl[i-1])) / (log10(pressure_h[id * nlay + i]) - log10(pressure_h[id * nlay + i-1]));
        gradad[i] = ((double)0.32) - ((double)0.10) * Tl[i] / ((double)3000.0);
    }

    gradrad[0] = 0.0;
    gradad[0] = 0.0;

    iRC = 1;
    iRC1 = 1;

    for (i = 1; i <= (nlay-1) ; i++)
    {
        if (iRC1 >= i -1)
        {
            if (gradrad[i] > ((double)0.7) * gradad[i])
            {
                iRC1 = i;
            }
            if (gradrad[i] > ((double)0.98) * gradad[i])
            {
                iRC = i;
            }
        }
    }

    if (iRC > -1 )
    {
        for (i = iRC; i > 0; i--)
        {
            gradad[i] = (double)0.32 - ((double)0.10) * Tl[id * nlay + i] / ((double)3000.0);
            if (gradad[i] < 0.0)
            {
                gradad[i] = 0.0;
            }
            
            Tl[id * nlay + i-1] = Tl[id * nlay + i] * pow((pressure_h[id * nlay + i-1] ) / pressure_h[id * nlay + i], gradad[i]);
        }
    }
}