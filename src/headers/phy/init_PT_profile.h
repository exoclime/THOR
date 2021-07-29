#pragma once
#include <math.h>

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
void Parmentier_IC(const int nlay, double* pl, double Tint, double mu, double Tirr,
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
    double kRoss[40];
    double tau[40 + 1]={0};

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
    k_Ross_Freedman(Tskin, pl[nlay-1], met, kRoss[nlay-1]);
    


    // Recalculate the upmost tau with new kappa
    tau[nlay-1] = kRoss[nlay-1] / grav * pl[nlay-1];

    // More accurate layer T at uppermost layer
    summy = 0.0;
    for (i = 0; i < 3; i++)
    {
        summy += 3.0 * Beta_V[i] * pow(Tmu, 4.0) / 4.0 * (C[i] + D[i] * exp(-tau[nlay-1] / tau_lim) +
            E[i] * exp(-gam_V[i] * tau[nlay-1]));
    }
    Tl[nlay-1] = 3.0 * pow(Tint, 4) / 4.0 * (tau[nlay-1] + A + B * exp(-tau[nlay-1] / tau_lim)) + summy;
    Tl[nlay-1] = pow(Tl[nlay-1], (1.0 / 4.0));

    // Now we can loop in optical depth space to find the T-p profile
    for (i = nlay-2; nlay>-1; i--)
    {
        // Initial guess for layer
        k_Ross_Freedman(Tl[i+1], sqrt(pl[i+1] * pl[i]), met, kRoss[i]);
        
        tau[i] = tau[i+1] + kRoss[i] / grav *  (pl[i] - pl[i+1]) ;
    
        summy = 0.0;
        for (j = 0; j < 3; j++)
        {
            summy = +3.0 * Beta_V[j] * pow(Tmu, 4.0) / 4.0 * (C[j] + D[j] * exp(-tau[i] / tau_lim) +
                E[j] * exp(-gam_V[j] * tau[i]));
        }
        Tl[i] = 3.0 * pow(Tint, 4.0) / 4.0 * (tau[i] + A + B * exp(-tau[i] / tau_lim)) + summy;
        Tl[i] = pow(Tl[i], (1.0 / 4.0));

        // Convergence loop
        for (j = 0; j < 5; j++)
        {
            k_Ross_Freedman(sqrt(Tl[i+1] * Tl[i]), sqrt(pl[i+1] * pl[i]), met, kRoss[i]);
            
            tau[i] = tau[i+1] + kRoss[i] / grav *  (pl[i] - pl[i+1]);
            summy = 0.0;
            for (k = 0; k < 3; k++)
            {
                summy += 3.0 * Beta_V[k] * pow(Tmu, 4.0) / 4.0 * (C[k] + D[k] * exp(-tau[i] / tau_lim) +
                    E[k] * exp(-gam_V[k] * tau[i]));
            }
            Tl[i] = 3.0 * pow(Tint, 4.0) / 4.0 * (tau[i] + A + B * exp(-tau[i] / tau_lim)) + summy;
            Tl[i] = pow(Tl[i], (1.0 / 4.0));
        }
    }


}

///////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

// Subroutine that corrects for adiabatic region following Parmentier & Guillot (2015)
void adiabat_correction(int nlay, double* (&Tl), double* pressure_h, double Gravit) {
    // dependcies
    //// main_parameters::nlay  -> "FMS_RC_para_&_const.cpp" 
    //// pow -> math    
    /// log10 -> math

    // Input:
    // 

    // Call by reference (Input & Output):
    // 

    // work variables
    int i, iRC, iRC1, d_p;
    double gradrad[nlay-1];
    double gradad[nlay-1];


    // start operations

    for (i = (nlay-2); i > 0; i--)
    {
        
        gradrad[i+1] = (log10(Tl[i+1]) - log10(Tl[i])) / (log10(pressure_h[i+1]) - log10(pressure_h[i]));
        gradad[i] = ((double)0.32) - ((double)0.10) * Tl[i] / ((double)3000.0);
    }

    gradrad[0] = 0.0;
    gradad[0] = 0.0;

    iRC = 1;
    iRC1 = 1;

    for (i = 0; i < (nlay-1) ; i++)
    {
        if (iRC1 <= i + 1)
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

    if (iRC > 0 )
    {
        for (i = iRC; i > 1; i--)
        {
            gradad[i] = (double)0.32 - ((double)0.10) * Tl[i] / ((double)3000.0);
            if (gradad[i] < 0.0)
            {
                gradad[i] = 0.0;
            }
            
            Tl[i] = Tl[i+1] * pow((pressure_h[i * nlay + i] ) / pressure_h[i * nlay + i+1], gradad[i+1]);
        }
    }
}