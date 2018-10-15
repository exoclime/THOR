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
// Build the model's grid
//
//
// Description:- The code first generates a standard icosahedral subdividing lower order grid.
//               The Platonic icosahedron is the lowest order grid.
//             - The spring dynamics method smoothes the grid distortions.
//             - Position of the centroids is corrected at the end.
//
// Method: - [1] - The standard grid is obtained by dividing the triangles of the icosahedral grid
//                 until the desired resolution is achieved.
//           [2] - Spring dynamics from Tomita et al. 2001.
//
//
// Known limitations:
//   - None.
//
// Known issues:
//   - None.
//
// Current Code Owner: Joao Mendonca, EEG. joao.mendonca@csh.unibe.ch
//
// If you use this code please cite the following reference:
//
// [1] - João M. Mendonça, Simon L. Grimm, Luc Grosheintz and Kevin Heng, 2016, Apj,
//       "THOR: A New and Flexible Global Circulation Model to Explore Planetary Atmospheres",
//       http://arxiv.org/abs/1607.05535
//
// [2] - Hirofumi Tomita, Motohiko Tsugawa, Masaki Satoh and Koji Goto Koji, 2001, Journal of Computational Physics
//       "Shallow Water Model on a Modified Icosahedral Geodesic Grid by Using Spring Dynamics",
//       http://adsabs.harvard.edu/abs/2001JCoPh.174..579T
//
//
// History:
// Version Date       Comment
// ======= ====       =======
//
// 1.0     16/08/2017 Released version  (JM)
//
////////////////////////////////////////////////////////////////////////

#include "../headers/grid.h"
#include "vector_operations.h"

// some helper local functions
inline double3 normproj(const double3 & v1, const double3 & v2)
{
  
    
    double l1 = dot(v1,v2);
    
    double3 vc1 = cross(v1,v2);

    double l2 = length(vc1);

    vc1 = vc1/l2*atan2(l2, l1);

    return vc1;
}



__host__ Icogrid::Icogrid (bool sprd         ,  // Spring dynamics option
                           double spring_beta,  // Parameter beta for spring dynamics
                           int glevel        ,  // Horizontal resolution level
                           int vlevel        ,  // Number of vertical layers
                           int nlat          ,
                           double A          ,  // Planet radius [m]
                           double Top_altitude, // Top model's domain [m]
                           bool sponge){

    printf("\n\n Building icosahedral grid!");

//  Number of vertices
    point_num = 2 + 10*pow(2.0,2.0*glevel); // Horizontal
    nv        = vlevel                    ; // Vertical
    nvi       = vlevel + 1                ; // Interfaces between layers

//  Number of times the main rhombi are divided.
    int divide_face; // Used to split memory on the GPU
    if (glevel == 4) divide_face = 0;
    else if (glevel == 5) divide_face = 1;
    else if (glevel == 6) divide_face = 2;
    else if (glevel == 7) divide_face = 3;
    else {
        printf("\nParameter not tested! Use the predefined divide_face values.\n");
        exit(EXIT_FAILURE);
    }
    grid_level = glevel;

//  Rhombi
    int n_region = pow(2.0,glevel)*pow(2.0,glevel) ; //
    nl_region    = pow(2.0,glevel)                 ; //
    nl_region    = nl_region/(pow(2.0,divide_face)); //
    int nl2      = nl_region*nl_region             ; //
    int nfaces   = pow(4.0,divide_face)            ; //
    int nlhalo   = nl_region+2                     ; //
    int kxl      = sqrt(nfaces*1.0)                ; //
    nr           = (point_num-2)/nl2               ; //

//  Compute standard grid.
    point_xyz  = (double*)malloc(3*point_num * sizeof(double));
    pent_ind   = (int*)malloc(12 * sizeof(int));
    sphere_ico (point_xyz  ,
                glevel     ,
                n_region   ,
                nl_region  ,
                nl2        ,
                kxl        ,
                nfaces     ,
                pent_ind   ,
                divide_face,
                point_num  );

//  Finds the closest neighbors of each point.
    point_local  = (int*)malloc(6*point_num * sizeof(int));
    neighbors_indx (point_local,
                    point_xyz  ,
                    pent_ind   ,
                    point_num  );

//  Reorder neighbors in the clockwise order.
    reorder_neighbors_indx (point_local,
                            point_xyz  ,
                            pent_ind   ,
                            point_num  );

//  Generate halos.
    nh = 10*nfaces*4*nlhalo;
    halo = (int*)malloc(10*nfaces*4*nlhalo * sizeof(int));
    generate_halos(halo       ,
                   point_local,
                   n_region   ,
                   divide_face);

//  Reorder neighbors consistent with the rhombi.
    reorder_neighbors_indx_rhombi (point_local,
                                   halo       ,
                                   pent_ind   ,
                                   nl_region  ,
                                   nl2        ,
                                   nr         ,
                                   nlhalo     ,
                                   point_num  );

//  Finds the closest neighbors at the pole.
    neighbors_indx_pl(point_local,
                      point_xyz  ,
                      pent_ind   ,
                      point_num  );

//  Reorder neighbors in the clockwise order at the pole.
    reorder_neighbors_indx_pl(point_local,
                              point_xyz  ,
                              pent_ind   ,
                              point_num  );

//  Produce rhombus' maps.
    maps = (int*)malloc( (nl_region+2)*(nl_region+2)*nr * sizeof(int));
    produce_maps(maps     ,
                 halo     ,
                 nr       ,
                 nl_region);

//  Smooths the grid applying the spring dynamic method.
    if (sprd){
        // Applies spring dynamics.
        spring_dynamics(point_local,
                        pent_ind   ,
                        glevel     ,
                        spring_beta,
                        point_xyz  ,
                        point_num  );

        //  Finds the q points.
        point_xyzq = (double*)malloc(6*3*point_num * sizeof(double));
        find_qpoints (point_local ,
                      point_xyzq  ,
                      point_xyz   ,
                      pent_ind    ,
                      point_num   );

        // Fixes the position of the centroids.
        relocate_centres(point_local ,
                         point_xyzq  ,
                         point_xyz   ,
                         pent_ind    ,
                         point_num   );

    }
    else{
//      Finds the q points.
        point_xyzq = (double*)malloc(6*3*point_num * sizeof(double));
        find_qpoints (point_local ,
                      point_xyzq  ,
                      point_xyz   ,
                      pent_ind    ,
                      point_num   );
    }

//  Recompute cartesians points for the new radius.
    correct_xyz_points ( A          ,
                         point_xyzq ,
                         point_xyz  ,
                         pent_ind   ,
                         point_num  );

//  Radial vectors with norm 1.
    func_r   = (double*)malloc(3*point_num * sizeof(double));
    compute_func(func_r   ,
                 point_xyz,
                 point_num);

//  Compute control areas.
    areas  = (double*)malloc(6*3*point_num * sizeof(double));
    areasTr= (double*)malloc(6 * point_num * sizeof(double));
    areasT = (double*)malloc(point_num * sizeof(double));
    control_areas( areasT     ,
                   areasTr    ,
                   areas      ,
                   point_local,
                   point_xyzq ,
                   point_xyz  ,
                   pent_ind   ,
                   point_num  );

//  Computes control vectors.
    nvec    = (double*)malloc(6 * 3 * point_num * sizeof(double));
    nvecoa  = (double*)malloc(6 * 3 * point_num * sizeof(double));
    nvecti  = (double*)malloc(6 * 3 * point_num * sizeof(double));
    nvecte  = (double*)malloc(6 * 3 * point_num * sizeof(double));
    control_vec(nvec       ,
                nvecoa     ,
                nvecti     ,
                nvecte     ,
                areasT     ,
                point_xyz  ,
                point_xyzq ,
                point_local,
                pent_ind   ,
                point_num  );

//  Set the Altitudes
    Altitude  = (double*)malloc(nv * sizeof(double));
    Altitudeh = (double*)malloc(nvi * sizeof(double));
    set_altitudes(Altitude    ,
                  Altitudeh   ,
                  Top_altitude,
                  nv         );

//  Converting to spherical coordinates.
    lonlat  = (double*)malloc(2*point_num * sizeof(double));
    cart2sphe ( lonlat    ,
                point_xyz ,
                point_num );

//  Computes the divergence operator.
    div = (double*)malloc(7 * 3 * point_num * sizeof(double));
    div_operator(areasT    ,
                 areas     ,
                 div       ,
                 nvec      ,
                 pent_ind  ,
                 point_num);

//  Computes the gradiente operator.
    grad = (double*)malloc(7 * 3 * point_num * sizeof(double));
    gra_operator(areasT   ,
                 areas    ,
                 grad     ,
                 nvec     ,
                 pent_ind ,
                 point_num);

//  Computes zonal mean for sponge layer operations
    if (sponge == true) {
      zonal_mean_tab = (int*)malloc(2 * point_num * sizeof(int));
      zonal_mean_tab_f(zonal_mean_tab    ,
                     lonlat            ,
                     nlat              ,
                     point_num         );
    }

    printf(" GRID DONE!\n\n");
}


void Icogrid::sphere_ico (double *xyz        ,
                          int     glevel     ,
                          int     n_region   ,
                          int     nl_region  ,
                          int     nl2        ,
                          int     kxl        ,
                          int     nfaces     ,
                          int    *pent_ind   ,
                          int     divide_face,
                          int     num        ){

//
//  Description:
//
//  Returns standard icosahedral grid points on a sphere with radius 1.
//
//  Input: - num         - Number of vertices.
//         - glevel      - Number of recursive iterations to increase horizontal resolution.
//         - pent_ind    - Pentagons' indexes.
//         - divide_face - Number of times the main rhombus are divided.
//
//  Output: - xyz       - the vertices coordinates (Cartesian; radius 1).
//          - nl_region - nl_region^2 is the number of points in the faces.
//          - kxl       - kxl^2 is the number of small rhombi inside the main rhombi.

//  Local variables
    int sizei = pow(2.0, glevel) ;
    int count_points             ;
    int sizeip1 = sizei+1        ;
    int sizei2  = sizeip1*sizeip1;

//  Temporary rhombus indexes.
    int *rhombi;
    rhombi = new int[10*sizei2]();
    int *rhomb;
    rhomb = new int[10*sizei*sizei]();

//  Temporary main vertices.
    double3 *xyzi;
    xyzi = new double3[(10*sizei2+2)]();

    double w = 2.0 * acos(1.0/(2.0*sin(M_PI/5.0)));

//  First icosahedron coordinates (pentagons)
//  Poles
//  North
    xyzi[0] = make_double3( 0.0, 0.0, 1.0);
// South
    xyzi[1] = make_double3(0.0, 0.0, -1.0);
//  Other points of the icosahedron.
//  3
    xyzi[2] = make_double3(cos(-M_PI/5.0)*cos(M_PI/2.0-w),
                           sin(-M_PI/5.0)*cos(M_PI/2.0-w),
                           sin(M_PI/2.0-w));
    rhombi[0] = 2;
//  4
    xyzi[3] = make_double3(cos(M_PI/5.0)*cos(M_PI/2.0-w),
                           sin(M_PI/5.0)*cos(M_PI/2.0-w),
                           sin(M_PI/2.0-w)),
    rhombi[1] = 3;
//  5
    xyzi[4] = make_double3(cos(3.0*M_PI/5.0)*cos(M_PI/2.0-w),
                           sin(3.0*M_PI/5.0)*cos(M_PI/2.0-w),
                           sin(M_PI/2-w));
    rhombi[2] = 4;
//  6
    xyzi[5] = make_double3(cos(M_PI)*cos(M_PI/2.0-w),
                           sin(M_PI)*cos(M_PI/2.0-w),
                           sin(M_PI/2.0-w));
    rhombi[3] = 5;
//  7
    xyzi[6] = make_double3(cos(-(3.0/5.0)*M_PI)*cos(M_PI/2.0-w),
                           sin(-(3.0/5.0)*M_PI)*cos(M_PI/2.0-w),
                           sin(M_PI/2.0-w));
    rhombi[4] = 6;
//  8
    xyzi[7] = make_double3(cos(0.0)*cos(w-M_PI/2.0),
                          sin(0.0)*cos(w-M_PI/2.0),
                          sin(w-M_PI/2.0));
    rhombi[5] = 7;
//  9
    xyzi[8] = make_double3(cos(2.0*M_PI/5.0)*cos(w-M_PI/2.0),
                           sin(2.0*M_PI/5.0)*cos(w-M_PI/2.0),
                           sin(w-M_PI/2.0));
    rhombi[6] = 8;
//  10
    xyzi[9] = make_double3(cos(4.0*M_PI/5.0)*cos(w-M_PI/2.0),
                           sin(4.0*M_PI/5.0)*cos(w-M_PI/2.0),
                           sin(w-M_PI/2.0));
    rhombi[7] = 9;
//  11
    xyzi[10] = make_double3(cos(-4.0*M_PI/5.0)*cos(w-M_PI/2.0),
                            sin(-4.0*M_PI/5.0)*cos(w-M_PI/2.0),
                            sin(w-M_PI/2.0));
    rhombi[8] = 10;
//  12
    xyzi[11] = make_double3(cos(-2.0*M_PI/5.0)*cos(w-M_PI/2.0),
                            sin(-2.0*M_PI/5.0)*cos(w-M_PI/2.0),
                            sin(w-M_PI/2.0));
    rhombi[9] = 11;
//
//  Standard grid points.
//
    int rhombi_points [10][4] = {
            {  2,  0,  7,  3 },
            {  3,  0,  8,  4 },
            {  4,  0,  9,  5 },
            {  5,  0, 10,  6 },
            {  6,  0, 11,  2 },
            {  7,  3,  1,  8 },
            {  8,  4,  1,  9 },
            {  9,  5,  1, 10 },
            { 10,  6,  1, 11 },
            { 11,  2,  1,  7 }
        };
    
    for (int faces = 0; faces < 10; faces++){
            rhombi[faces]                               = rhombi_points[faces][0];
            rhombi[sizei*sizeip1*10 + faces]            = rhombi_points[faces][1];
            rhombi[sizei*10 + faces]                    = rhombi_points[faces][2];
            rhombi[sizei*sizeip1*10 + sizei*10 + faces] = rhombi_points[faces][3];
    }

    // Recursive method
    count_points = 12;
    for (int faces = 0; faces < 10; faces++){
        for (int grd = 0; grd < glevel; grd++){
            int ind_jump = pow(2.0,glevel-grd-1);
            int ind      = ind_jump;
            for (int it = 0; it < pow(2.0,grd); it++){
                int indx = ind;
                int indy = 0  ;
                rhombi[indx*sizeip1*10 + indy*10 + faces] = count_points;
                int indx1 = ind+ind_jump;
                int indy1 = 0  ;
                int indx2 = ind-ind_jump;
                int indy2 = 0  ;

                int r_idx_1 = rhombi[indx1*sizeip1*10 + indy1*10 + faces];
                int r_idx_2 = rhombi[indx2*sizeip1*10 + indy2*10 + faces];
                                      
                xyzi[count_points] = normalize( (xyzi[r_idx_1] + xyzi[r_idx_2])*0.5);
                count_points += 1;

                indx = 0;
                indy = ind  ;
                rhombi[indx*sizeip1*10 + indy*10 + faces] = count_points;
                indx1 = 0;
                indy1 = ind+ind_jump;
                indx2 = 0;
                indy2 = ind-ind_jump;

                r_idx_1 = rhombi[indx1*sizeip1*10 + indy1*10 + faces];
                r_idx_2 = rhombi[indx2*sizeip1*10 + indy2*10 + faces];
                                      
                xyzi[count_points] = normalize( (xyzi[r_idx_1] + xyzi[r_idx_2])*0.5);
                
                count_points += 1;

                indx = ind;
                indy = ind;
                rhombi[indx*sizeip1*10 + indy*10 + faces] = count_points;
                indx1 = ind+ind_jump;
                indy1 = ind+ind_jump;
                indx2 = ind-ind_jump;
                indy2 = ind-ind_jump;
                r_idx_1 = rhombi[indx1*sizeip1*10 + indy1*10 + faces];
                r_idx_2 = rhombi[indx2*sizeip1*10 + indy2*10 + faces];
                                      
                xyzi[count_points] = normalize( (xyzi[r_idx_1] + xyzi[r_idx_2])*0.5);

                count_points += 1;

                indx = ind;
                indy = sizei;
                rhombi[indx*sizeip1*10 + indy*10 + faces] = count_points;
                indx1 = ind+ind_jump;
                indy1 = sizei;
                indx2 = ind-ind_jump;
                indy2 = sizei;
                r_idx_1 = rhombi[indx1*sizeip1*10 + indy1*10 + faces];
                r_idx_2 = rhombi[indx2*sizeip1*10 + indy2*10 + faces];
                                      
                xyzi[count_points] = normalize( (xyzi[r_idx_1] + xyzi[r_idx_2])*0.5);

                count_points += 1;

                indx = sizei;
                indy = ind;
                rhombi[indx*sizeip1*10 + indy*10 + faces] = count_points;
                indx1 = sizei;
                indy1 = ind+ind_jump;
                indx2 = sizei;
                indy2 = ind-ind_jump;
                
                r_idx_1 = rhombi[indx1*sizeip1*10 + indy1*10 + faces];
                r_idx_2 = rhombi[indx2*sizeip1*10 + indy2*10 + faces];
                                      
                xyzi[count_points] = normalize( (xyzi[r_idx_1] + xyzi[r_idx_2])*0.5);

                count_points += 1;

                ind += 2*ind_jump;
            }
        }
        for (int grd = 0; grd < glevel; grd++){
            if(grd>0){
                // Upper triangle
                int ind_jump = pow(2.0,glevel-grd-1);
                int indj     = ind_jump;
                for (int j = 0; j < pow(2.0,grd)-1; j++){
                    int indy = ind_jump;
                    int indx = 0       ;
                    for (int i = 0; i < j+1; i++){
                        int rindy = indy;
                        int rindx = sizei-indj+indx;

                        rhombi[rindx*sizeip1*10 + rindy*10 + faces] = count_points;
                        int indx1 = rindx-ind_jump;
                        int indy1 = rindy-ind_jump;
                        int indx2 = rindx+ind_jump;
                        int indy2 = rindy+ind_jump;

                        int r_idx_1 = rhombi[indx1*sizeip1*10 + indy1*10 + faces];
                        int r_idx_2 = rhombi[indx2*sizeip1*10 + indy2*10 + faces];
                        
                        xyzi[count_points] = normalize((xyzi[r_idx_1] + xyzi[r_idx_2])*0.5);
                        
                        count_points += 1;

                        rhombi[(rindx-ind_jump)*sizeip1*10 + rindy*10 + faces] = count_points;
                        indx1 = rindx-ind_jump;
                        indy1 = rindy-ind_jump;
                        indx2 = rindx-ind_jump;
                        indy2 = rindy+ind_jump;
                        r_idx_1 = rhombi[indx1*sizeip1*10 + indy1*10 + faces];
                        r_idx_2 = rhombi[indx2*sizeip1*10 + indy2*10 + faces];
                        
                        xyzi[count_points] = normalize((xyzi[r_idx_1] + xyzi[r_idx_2])*0.5);

                        count_points += 1;

                        rhombi[rindx*sizeip1*10 + (rindy+ind_jump)*10 + faces] = count_points;
                        indx1 = rindx+ind_jump;
                        indy1 = rindy+ind_jump;
                        indx2 = rindx-ind_jump;
                        indy2 = rindy+ind_jump;
                        r_idx_1 = rhombi[indx1*sizeip1*10 + indy1*10 + faces];
                        r_idx_2 = rhombi[indx2*sizeip1*10 + indy2*10 + faces];
                        
                        xyzi[count_points] = normalize((xyzi[r_idx_1] + xyzi[r_idx_2])*0.5);

                        count_points += 1;

                        indx += 2*ind_jump;
                        indy += 2*ind_jump;
                    }
                    indj += 2*ind_jump;
                }
            }
            if(grd>0){
                // Lower triangle
                int ind_jump = pow(2.0,glevel-grd-1);
                int indj      = ind_jump;
                for (int j = 0; j < pow(2.0,grd)-1; j++){
                    int indy = 0;
                    int indx = ind_jump;
                    for (int i = 0; i < j+1; i++){

                        int rindx = indx;
                        int rindy = sizei-indj+indy;

                        rhombi[rindx*sizeip1*10 + rindy*10 + faces] = count_points;
                        int indx1 = rindx-ind_jump;
                        int indy1 = rindy-ind_jump;
                        int indx2 = rindx+ind_jump;
                        int indy2 = rindy+ind_jump;

                        int r_idx_1 = rhombi[indx1*sizeip1*10 + indy1*10 + faces];
                        int r_idx_2 = rhombi[indx2*sizeip1*10 + indy2*10 + faces];
                        
                        xyzi[count_points] = normalize((xyzi[r_idx_1] + xyzi[r_idx_2])*0.5);
                        count_points += 1;

                        rhombi[rindx*sizeip1*10 + (rindy-ind_jump)*10 + faces] = count_points;
                        indx1 = rindx-ind_jump;
                        indy1 = rindy-ind_jump;
                        indx2 = rindx+ind_jump;
                        indy2 = rindy-ind_jump;

                        r_idx_1 = rhombi[indx1*sizeip1*10 + indy1*10 + faces];
                        r_idx_2 = rhombi[indx2*sizeip1*10 + indy2*10 + faces];
                        
                        xyzi[count_points] = normalize((xyzi[r_idx_1] + xyzi[r_idx_2])*0.5);
                       
                        count_points += 1;

                        rhombi[(rindx+ind_jump)*sizeip1*10 + rindy*10 + faces] = count_points;
                        indx1 = rindx+ind_jump;
                        indy1 = rindy+ind_jump;
                        indx2 = rindx+ind_jump;
                        indy2 = rindy-ind_jump;
                        r_idx_1 = rhombi[indx1*sizeip1*10 + indy1*10 + faces];
                        r_idx_2 = rhombi[indx2*sizeip1*10 + indy2*10 + faces];
                        
                        xyzi[count_points] = normalize((xyzi[r_idx_1] + xyzi[r_idx_2])*0.5);
                       
                        count_points += 1;

                        indx += 2*ind_jump;
                        indy += 2*ind_jump;
                    }
                indj += 2*ind_jump;
                }
            }
        }
    }

    int nli_region  = sqrt(n_region*1.0);

    for (int fc = 0; fc < 10; fc++)    for (int j = 0; j < sizei; j++) for (int i = 0; i < sizei; i++)
        rhomb[i*sizei*10 + j*10 + fc] = rhombi[i*sizeip1*10 + j*10 + fc];

    for (int fc = 0; fc < 10; fc++)    for (int kx = 0; kx < kxl; kx++) for (int ky = 0; ky < kxl; ky++)
        for (int i = 0; i < nl_region; i++)    for (int j = 0; j < nl_region; j++){
             int idx1 = fc*nfaces*nl2 + ky*kxl*nl2 + kx*nl2 + j*nl_region + i;
             int idx2 = rhomb[((ky*nl_region + j)*nli_region + kx*nl_region + i)*10 + fc];
                
             xyz[idx1*3 + 0] = xyzi[idx2].x;
             xyz[idx1*3 + 1] = xyzi[idx2].y;
             xyz[idx1*3 + 2] = xyzi[idx2].z;
            }
            
    
    

    //North
    xyz[(num-2)*3 + 0] =0.0;
    xyz[(num-2)*3 + 1] =0.0;
    xyz[(num-2)*3 + 2] =1.0;

    //South
    xyz[(num-1)*3 + 0] =0.0 ;
    xyz[(num-1)*3 + 1] =0.0 ;
    xyz[(num-1)*3 + 2] =-1.0;

    //Pentagons' indexes
    for (int faces = 0; faces < 10; faces++) pent_ind[faces] = faces*nfaces*nl2;
    pent_ind[10] = num-2;
    pent_ind[11] = num-1;

    delete [] rhombi;
    delete [] rhomb ;
    delete [] xyzi  ;
}

void Icogrid::neighbors_indx (int *point_local,
                              double *xyz     ,
                              int *pent_ind   ,
                              int point_num   ){

//
//  Description:
//
//  Finds the first neighbors of each vertex.
//
//  Input: xyz       - Vertices coordinates.
//         point_num - Number of vertices.
//         pent_ind  - Pentagons indexes.
//
//  Output: point_local- First neighbors indexes of each vertex.
//

//  Local variables.
    int position;
    double small;
    double *distance;
    distance = new double[point_num]();

    for (int i = 0; i < 6*point_num; i++) point_local[i] = -1;

//  Find neighbors.
    for (int i = 0; i < point_num; i++){
        for (int j = 0; j < point_num; j++){
                    distance[j] = sqrt(pow(xyz[0 + j*3] - xyz[0 + i*3],2) +
                                       pow(xyz[1 + j*3] - xyz[1 + i*3],2) +
                                       pow(xyz[2 + j*3] - xyz[2 + i*3],2));
        }
        for(int k = 0; k < 6; k++){  // Hexagons.
            small    =  2;
            position = -1;
            for (int j = 0; j < point_num; j++){
                if(k == 0){
                    if(small >= distance[j] && distance[j] > 0){
                        small    = distance[j];
                        position = j;
                    }
                }
                else{
                    if(small >= distance[j] && distance[j] > 0){
                        if(point_local[i*6 + 0] != j && point_local[i*6 + 1] != j &&
                           point_local[i*6 + 2] != j && point_local[i*6 + 3] != j &&
                           point_local[i*6 + 4] != j && point_local[i*6 + 5] != j){
                           small = distance[j];
                           position = j;
                        }
                    }
                }
            }
            point_local[i*6 + k] = position;
        }
        for(int k = 0; k < 12; k++){
            if(i==pent_ind[k]){ // Pentagons.
                point_local[i*6 + 5] = -1;
            }
        }
    }

    delete [] distance;
}

void Icogrid::neighbors_indx_pl(int *point_local,
                                double *xyz     ,
                                int *pent_ind   ,
                                int point_num   ){

//
//  Description:
//
//  Finds the first neighbors of each pole.
//
//  Input: xyz       - Vertices coordinates.
//         point_num - Number of vertices.
//         pent_ind  - Pentagons indexes.
//
//  Output: point_local- First neighbors indexes of each vertex.
//

//  Local variables.
    int position;
    double small;
    double *distance;
    distance = new double[point_num]();

    //  Find neighbors.
    for (int i = point_num-2; i < point_num; i++){
        for (int j = 0; j < point_num; j++){
            distance[j] = sqrt(pow(xyz[0 + j * 3] - xyz[0 + i * 3], 2) +
                pow(xyz[1 + j * 3] - xyz[1 + i * 3], 2) +
                pow(xyz[2 + j * 3] - xyz[2 + i * 3], 2));
        }
        for (int k = 0; k < 5; k++){
            small = 2;
            position = -1;
            for (int j = 0; j < point_num; j++){
                if (k == 0){
                    if (small >= distance[j] && distance[j] > 0){
                        small = distance[j];
                        position = j;
                    }
                }
                else{
                    if (small >= distance[j] && distance[j] > 0){
                        if (point_local[i * 6 + 0] != j && point_local[i * 6 + 1] != j &&
                            point_local[i * 6 + 2] != j && point_local[i * 6 + 3] != j &&
                            point_local[i * 6 + 4] != j && point_local[i * 6 + 5] != j){
                            small = distance[j];
                            position = j;
                        }
                    }
                }
            }
            point_local[i * 6 + k] = position;
        }
    }

    delete[] distance;
}

void Icogrid::reorder_neighbors_indx (int *point_local,
                                      double *xyz     ,
                                      int *pent_ind   ,
                                      int point_num   ){

//
//  Description:
//
//  Re-order neighbors indexes in the clockwise direction.
//
//  Input: xyz        - Vertices coordinates.
//         point_num  - Number of vertices.
//         point_local- First neighbors indexes of each vertex.
//         pent_ind   - Pentagons indexes.
//
//  Output:POINT_LOCAL- First neighbors indexes of each vertex (clockwise direction).
//

//  Local variables
    double x1, y1, z1;
    double dx21, dy21, dz21;
    double dx31, dy31, dz31;
    double n21, n31;
    double dircx, dircy, dircz, dircn;
    double temp, tempp;
    int geo;

//  Temporary arrays.
    double *ang         ;
    int    *point_local_new;
    ang   = new double[6]();
    point_local_new = new int[point_num*6];

    for (int i = 0; i < point_num; i++)    for (int j = 0; j < 6; j++) point_local_new[i*6 + j] = point_local[i*6 + j];

    // Reorder indexes.
    for (int i = 0; i < point_num; i++){
        geo = 6; // Hexagons.
        for (int k = 0; k < 12; k++) if(i == pent_ind[k]) geo = 5; // Pentagons.
        if(geo == 5){
            for (int j = 0; j < 6; j++) ang[j] = 0;
            x1   = 0;
            y1   = 0;
            z1   = 0;
            for (int j = 0; j < 5; j++){ // Central coplanar point.
                x1 = x1 + xyz[point_local[i*6 + j]*3 + 0]/5;
                y1 = y1 + xyz[point_local[i*6 + j]*3 + 1]/5;
                z1 = z1 + xyz[point_local[i*6 + j]*3 + 2]/5;
            }
            dx21 = xyz[point_local[i*6 + 0]*3 + 0] - x1;
            dy21 = xyz[point_local[i*6 + 0]*3 + 1] - y1;
            dz21 = xyz[point_local[i*6 + 0]*3 + 2] - z1;
            for (int j = 1; j < 5; j++){
                dx31 = xyz[point_local[i*6 + j]*3 + 0] - x1;
                dy31 = xyz[point_local[i*6 + j]*3 + 1] - y1;
                dz31 = xyz[point_local[i*6 + j]*3 + 2] - z1;
                n21  = sqrt(dx21*dx21 + dy21*dy21 + dz21*dz21);
                n31  = sqrt(dx31*dx31 + dy31*dy31 + dz31*dz31);
                dircx= dy21*dz31 - dy31*dz21;
                dircy= dz21*dx31 - dz31*dx21;
                dircz= dx21*dy31 - dx31*dy21;
                dircn= dircx*x1  + dircy*y1 + dircz*z1;
                if(dircn <= 0){
                    ang[j] = acos((dx21*dx31 + dy21*dy31 + dz21*dz31)/(n21*n31));
                }
                else{
                    ang[j] = 2.0*M_PI - acos((dx21*dx31 + dy21*dy31 + dz21*dz31)/(n21*n31));
                }
            }
            for (int j = 0; j < 5; j++){
                for (int k = 4; k >= j; k--){
                    if(ang[j] > ang[k]){
                        temp                     = ang[k];
                        tempp                    = point_local_new[i*6 + k];
                        ang[k]                   = ang[j];
                        point_local_new[i*6 + k] = point_local_new[i*6 + j];
                        ang[j]                   = temp;
                        point_local_new[i*6 + j] = tempp;
                    }
                }
            }
        }
        else{
            for (int j = 0; j < 6; j++) ang[j] = 0;
            x1 = 0;
            y1 = 0;
            z1 = 0;
            for (int j = 0; j < 6; j++){ // Central coplanar point.
                x1 = x1 + xyz[point_local[i*6 + j]*3 + 0]/6;
                y1 = y1 + xyz[point_local[i*6 + j]*3 + 1]/6;
                z1 = z1 + xyz[point_local[i*6 + j]*3 + 2]/6;
            }
            dx21 = xyz[point_local[i*6 + 0]*3 + 0] - x1;
            dy21 = xyz[point_local[i*6 + 0]*3 + 1] - y1;
            dz21 = xyz[point_local[i*6 + 0]*3 + 2] - z1;
            for (int j = 1; j < 6; j++){
                dx31 = xyz[point_local[i*6 + j]*3 + 0] - x1;
                dy31 = xyz[point_local[i*6 + j]*3 + 1] - y1;
                dz31 = xyz[point_local[i*6 + j]*3 + 2] - z1;
                n21  = sqrt(dx21*dx21 + dy21*dy21 + dz21*dz21);
                n31  = sqrt(dx31*dx31 + dy31*dy31 + dz31*dz31);
                dircx= dy21*dz31 - dy31*dz21;
                dircy= dz21*dx31 - dz31*dx21;
                dircz= dx21*dy31 - dx31*dy21;
                dircn= dircx*x1  + dircy*y1 + dircz*z1;
                if(dircn <= 0)  ang[j] = acos((dx21*dx31 + dy21*dy31 + dz21*dz31)/(n21*n31));
                else ang[j] = 2.0*M_PI - acos((dx21*dx31 + dy21*dy31 + dz21*dz31)/(n21*n31));
            }
            for (int j = 0; j < 6; j++){
                for (int k = 5; k >= j; k--){
                    if(ang[j] > ang[k]){
                        temp                 = ang[k];
                        tempp                = point_local_new[i*6 + k];
                        ang[k]               = ang[j];
                        point_local_new[i*6 + k] = point_local_new[i*6 + j];
                        ang[j]               = temp;
                        point_local_new[i*6 + j] = tempp;
                    }
                }
            }
        }
    }

    for (int i = 0; i < point_num; i++)    for (int j = 0; j < 6; j++) point_local[i*6 + j] =point_local_new[i*6 + j];

    delete [] ang         ;
    delete [] point_local_new;
}

void Icogrid::reorder_neighbors_indx_pl(int *point_local,
                                        double *xyz     ,
                                        int *pent_ind   ,
                                        int point_num   ){

//
//  Description:
//
//  Re-order neighbors indexes in the clockwise direction (POLES).
//
//  Input: xyz        - Vertices coordinates.
//         point_num  - Number of vertices.
//         point_local- First neighbors indexes of each vertex.
//         pent_ind   - Pentagons indexes.
//
//  Output:POINT_LOCAL- First neighbors indexes of each vertex (clockwise direction).
//

//  Local variables
    double x1, y1, z1;
    double dx21, dy21, dz21;
    double dx31, dy31, dz31;
    double n21, n31;
    double dircx, dircy, dircz, dircn;
    double temp, tempp;

    //  Temporary arrays.
    double *ang;
    int    *point_local_new;
    ang = new double[6]();
    point_local_new = new int[point_num * 6];

    for (int i = 0; i < point_num; i++)    for (int j = 0; j < 6; j++) point_local_new[i * 6 + j] = point_local[i * 6 + j];

    // Reorder indexes.
    for (int i = point_num-2; i < point_num; i++){
        for (int j = 0; j < 6; j++) ang[j] = 0;
        x1 = 0;
        y1 = 0;
        z1 = 0;
        for (int j = 0; j < 5; j++){ // Central coplanar point.
            x1 = x1 + xyz[point_local[i * 6 + j] * 3 + 0] / 5;
            y1 = y1 + xyz[point_local[i * 6 + j] * 3 + 1] / 5;
            z1 = z1 + xyz[point_local[i * 6 + j] * 3 + 2] / 5;
        }
        dx21 = xyz[point_local[i * 6 + 0] * 3 + 0] - x1;
        dy21 = xyz[point_local[i * 6 + 0] * 3 + 1] - y1;
        dz21 = xyz[point_local[i * 6 + 0] * 3 + 2] - z1;
        for (int j = 1; j < 5; j++){
            dx31 = xyz[point_local[i * 6 + j] * 3 + 0] - x1;
            dy31 = xyz[point_local[i * 6 + j] * 3 + 1] - y1;
            dz31 = xyz[point_local[i * 6 + j] * 3 + 2] - z1;
            n21 = sqrt(dx21*dx21 + dy21*dy21 + dz21*dz21);
            n31 = sqrt(dx31*dx31 + dy31*dy31 + dz31*dz31);
            dircx = dy21*dz31 - dy31*dz21;
            dircy = dz21*dx31 - dz31*dx21;
            dircz = dx21*dy31 - dx31*dy21;
            dircn = dircx*x1 + dircy*y1 + dircz*z1;
            if (dircn <= 0){
                ang[j] = acos((dx21*dx31 + dy21*dy31 + dz21*dz31) / (n21*n31));
            }
            else{
                ang[j] = 2.0*M_PI - acos((dx21*dx31 + dy21*dy31 + dz21*dz31) / (n21*n31));
            }
        }
        for (int j = 0; j < 5; j++){
            for (int k = 4; k >= j; k--){
                if (ang[j] > ang[k]){
                    temp = ang[k];
                    tempp = point_local_new[i * 6 + k];
                    ang[k] = ang[j];
                    point_local_new[i * 6 + k] = point_local_new[i * 6 + j];
                    ang[j] = temp;
                    point_local_new[i * 6 + j] = tempp;
                }
            }
        }
    }

    for (int i = point_num-2; i < point_num; i++)for (int j = 0; j < 6; j++) point_local[i * 6 + j] = point_local_new[i * 6 + j];

    delete[] ang;
    delete[] point_local_new;
}

void Icogrid::generate_halos(int *halo       ,
                             int *point_local,
                             int  n_region   ,
                             int  divide_face){

//
//  Description:
//
//  Generate the halos around the rhombi.
//
//  Input: - point_local - First neighbors indexes of each vertex.
//         - n_region    - Number of points in a rhombus.
//         - divide_face - Number of times the main rhombus are divided.
//
//  Output: - halo - Indexes of the points that form the halos.
//

//  Local variables.
    int nl_region = sqrt(n_region*1.0)        ;
    nl_region = nl_region/pow(2.0,divide_face);
    int nl2 = nl_region*nl_region             ;
    int nfaces = pow(4.0, divide_face)        ;
    int kxl = sqrt(nfaces*1.0)                ;
    int nh = 4*(nl_region+2)                  ;

    /*
                  //\\
       Side 4    //%%\\  Side 3
                //%%%%\\
                \\%%%%//
       Side 1    \\%%//  Side 2
                  \\//
    */

    int index_value1, index_value2;
    int v1, v2;

    int out_ind;
    int ind_nei = 0;

    for (int i = 0; i < 10*nfaces*4*(nl_region+2); i++)halo[i] = -1;

    // Generate halos.
    for (int fc = 0; fc < 10; fc++){
        for (int kx = 0; kx < kxl; kx++){
            for (int ky = 0; ky < kxl; ky++){
                if(kx == 0 && ky == 0){ // First point is a pentagon
                    // Side 1
                    for (int i = 0; i < nl_region-1; i++){
                        int j = 0;
                        index_value1 = fc*nfaces*nl2 + ky*kxl*nl2 + kx*nl2 + j*nl_region + i;
                        index_value2 = fc*nfaces*nl2 + ky*kxl*nl2 + kx*nl2 + j*nl_region + i + 1;
                        if(kx == 0 && ky == 0 && i == 0 && j == 0){
                            for (int k = 0; k < 5; k++){
                                for (int kk = 0; kk < 6; kk++){
                                    v1 = point_local[index_value1*6 + k];
                                    v2 = point_local[index_value2*6 + kk];
                                    if(v1==v2){
                                        out_ind = 0;
                                        for (int kkk = 0; kkk < nl2; kkk++)    if(v1==fc*nfaces*nl2 + ky*kxl*nl2 + kx*nl2 + kkk) out_ind = 1;
                                        if(out_ind == 0)halo[fc*kxl*kxl*nh + ky*kxl*nh + kx*nh + i] = v1;
                                    }
                                }
                            }
                        }
                        else{
                            for (int k = 0; k < 6; k++){
                               for (int kk = 0; kk < 6; kk++){
                                   v1 = point_local[index_value1*6 + k];
                                   v2 = point_local[index_value2*6 + kk];
                                   if(v1==v2){
                                       out_ind = 0;
                                       for (int kkk = 0; kkk < nl2; kkk++)    if(v1==fc*nfaces*nl2 + ky*kxl*nl2 + kx*nl2 + kkk) out_ind = 1;
                                       if(out_ind == 0)halo[fc*kxl*kxl*nh + ky*kxl*nh + kx*nh + i] = v1;
                                   }
                               }
                            }
                        }
                    }
                    int i = nl_region -1;
                    int j = 0           ;
                    index_value1 = fc*nfaces*nl2 + ky*kxl*nl2 + kx*nl2 + j*nl_region + i;
                    for (int k = 0; k < 6; k++){
                        v1 = point_local[index_value1*6 + k];
                        if(v1== halo[fc*kxl*kxl*nh + ky*kxl*nh + kx*nh + i - 1])ind_nei = k;
                    }
                    if(ind_nei > 0) halo[fc*kxl*kxl*nh + ky*kxl*nh + kx*nh + i] = point_local[index_value1*6 + ind_nei-1];
                    else            halo[fc*kxl*kxl*nh + ky*kxl*nh + kx*nh + i] = point_local[index_value1*6 + 5];

                    //Side 2
                    i = nl_region -1;
                    j = 0           ;
                    index_value1 = fc*nfaces*nl2 + ky*kxl*nl2 + kx*nl2 + j*nl_region + i;
                    for (int k = 0; k < 6; k++){
                        v1 = point_local[index_value1*6 + k];
                        if(v1== halo[fc*kxl*kxl*nh + ky*kxl*nh + kx*nh + nl_region - 1])ind_nei = k;
                    }
                    if(ind_nei > 0) halo[fc*kxl*kxl*nh + ky*kxl*nh + kx*nh + nl_region +2] = point_local[index_value1*6 + ind_nei-1];
                    else            halo[fc*kxl*kxl*nh + ky*kxl*nh + kx*nh + nl_region +2] = point_local[index_value1*6 + 5];

                    i = nl_region-1;
                    for (int j = 0; j < nl_region-1; j++){
                         index_value1 = fc*nfaces*nl2 + ky*kxl*nl2 + kx*nl2 + j*nl_region + i;
                         index_value2 = fc*nfaces*nl2 + ky*kxl*nl2 + kx*nl2 + (j+1)*nl_region + i;
                         for (int k = 0; k < 6; k++){
                             for (int kk = 0; kk < 6; kk++){
                                 v1 = point_local[index_value1*6 + k];
                                 v2 = point_local[index_value2*6 + kk];
                                 if(v1==v2){
                                     int out_ind = 0;
                                     for (int kkk = 0; kkk < nl2; kkk++)    if(v1==fc*nfaces*nl2 + ky*kxl*nl2 + kx*nl2 + kkk) out_ind = 1;
                                     if(out_ind == 0)halo[fc*kxl*kxl*nh + ky*kxl*nh + kx*nh + (nl_region+2) + j+1] = v1;
                                 }
                             }
                         }
                    }

                    // Side 3
                    i = nl_region -1;
                    j = nl_region -1;
                    index_value1 = fc*nfaces*nl2 + ky*kxl*nl2 + kx*nl2 + j*nl_region + i;
                    for (int k = 0; k < 6; k++){
                        v1 = point_local[index_value1*6 + k];
                        if(v1 == halo[fc*kxl*kxl*nh + ky*kxl*nh + kx*nh + nl_region + 2 + nl_region - 1])ind_nei = k;
                    }
                    if(ind_nei > 0) halo[fc*kxl*kxl*nh + ky*kxl*nh + kx*nh + 2*(nl_region +2)] = point_local[index_value1*6 + ind_nei-1];
                    else            halo[fc*kxl*kxl*nh + ky*kxl*nh + kx*nh + 2*(nl_region +2)] = point_local[index_value1*6 + 5];

                    for (int i = 0; i < nl_region-1; i++){
                         j = nl_region-1;
                         index_value1 = fc*nfaces*nl2 + ky*kxl*nl2 + kx*nl2 + j*nl_region + nl_region - i-1;
                         index_value2 = fc*nfaces*nl2 + ky*kxl*nl2 + kx*nl2 + j*nl_region + nl_region - i-2;
                         for (int k = 0; k < 6; k++){
                             for (int kk = 0; kk < 6; kk++){
                                 v1 = point_local[index_value1*6 + k];
                                 v2 = point_local[index_value2*6 + kk];

                                 if(v1==v2){
                                     int out_ind = 0;
                                     for (int kkk = 0; kkk < nl2; kkk++)    if(v1==fc*nfaces*nl2 + ky*kxl*nl2 + kx*nl2 + kkk) out_ind = 1;
                                     if(out_ind == 0)halo[fc*kxl*kxl*nh + ky*kxl*nh + kx*nh + 2*(nl_region+2) + i+1] = v1;
                                 }
                             }
                         }
                    }

                    i = 0;
                    j = nl_region -1;
                    index_value1 = fc*nfaces*nl2 + ky*kxl*nl2 + kx*nl2 + j*nl_region + i;
                    for (int k = 0; k < 6; k++){
                        v1 = point_local[index_value1*6 + k];
                        if(v1 == halo[fc*kxl*kxl*nh + ky*kxl*nh + kx*nh + 2*(nl_region + 2) + nl_region - 1])ind_nei = k;
                    }
                    if(ind_nei > 0) halo[fc*kxl*kxl*nh + ky*kxl*nh + kx*nh + 2*(nl_region +2) + nl_region] = point_local[index_value1*6 + ind_nei-1];
                    else            halo[fc*kxl*kxl*nh + ky*kxl*nh + kx*nh + 2*(nl_region +2) + nl_region] = point_local[index_value1*6 + 5];

                    // Side 4
                    i = 0;
                    j = nl_region -1;
                    index_value1 = fc*nfaces*nl2 + ky*kxl*nl2 + kx*nl2 + j*nl_region + i;
                    for (int k = 0; k < 6; k++){
                        v1 = point_local[index_value1*6 + k];
                        if(v1 == halo[fc*kxl*kxl*nh + ky*kxl*nh + kx*nh + 2*(nl_region + 2) + nl_region ])ind_nei = k;
                    }
                    if(ind_nei > 0) halo[fc*kxl*kxl*nh + ky*kxl*nh + kx*nh + 3*(nl_region +2)] = point_local[index_value1*6 + ind_nei-1];
                    else            halo[fc*kxl*kxl*nh + ky*kxl*nh + kx*nh + 3*(nl_region +2)] = point_local[index_value1*6 + 5];

                    i = 0;
                    for (int j = 0; j < nl_region-1; j++){
                         index_value1 = fc*nfaces*nl2 + ky*kxl*nl2 + kx*nl2 + (nl_region - j-1)*nl_region;
                         index_value2 = fc*nfaces*nl2 + ky*kxl*nl2 + kx*nl2 + (nl_region - j-2)*nl_region;
                         for (int k = 0; k < 6; k++){
                             for (int kk = 0; kk < 6; kk++){
                                 v1 = point_local[index_value1*6 + k];
                                 v2 = point_local[index_value2*6 + kk];
                                 if(v1==v2){
                                     out_ind = 0;
                                     for (int kkk = 0; kkk < nl2; kkk++)    if(v1==fc*nfaces*nl2 + ky*kxl*nl2 + kx*nl2 + kkk) out_ind = 1;
                                     if(out_ind == 0)halo[fc*kxl*kxl*nh + ky*kxl*nh + kx*nh + 3*(nl_region+2) + j+1] = v1;
                                 }
                             }
                         }
                    }
                }
                else{//Hexagon
                    for (int i = 0; i < nl_region-1; i++){
                        int j = 0;
                        index_value1 = fc*nfaces*nl2 + ky*kxl*nl2 + kx*nl2 + j*nl_region + i;
                        index_value2 = fc*nfaces*nl2 + ky*kxl*nl2 + kx*nl2 + j*nl_region + i + 1;
                        if(kx == 0 && ky == 0 && i == 0 && j == 0){
                            for (int k = 0; k < 6; k++){
                                for (int kk = 0; kk < 6; kk++){
                                    v1 = point_local[index_value1*6 + k];
                                    v2 = point_local[index_value2*6 + kk];
                                    if(v1==v2){
                                        int out_ind = 0;
                                        for (int kkk = 0; kkk < nl2; kkk++)    if(v1==fc*nfaces*nl2 + ky*kxl*nl2 + kx*nl2 + kkk) out_ind = 1;
                                        if(out_ind == 0)halo[fc*kxl*kxl*nh + ky*kxl*nh + kx*nh + i] = v1;
                                    }
                                }
                            }
                        }
                        else{
                            for (int k = 0; k < 6; k++){
                               for (int kk = 0; kk < 6; kk++){
                                   v1 = point_local[index_value1*6 + k];
                                   v2 = point_local[index_value2*6 + kk];
                                   if(v1==v2){
                                       out_ind = 0;
                                       for (int kkk = 0; kkk < nl2; kkk++)    if(v1==fc*nfaces*nl2 + ky*kxl*nl2 + kx*nl2 + kkk) out_ind = 1;
                                       if(out_ind == 0)halo[fc*kxl*kxl*nh + ky*kxl*nh + kx*nh + i] = v1;
                                   }
                               }
                            }
                        }
                    }
                    int i = nl_region -1;
                    int j = 0           ;
                    index_value1 = fc*nfaces*nl2 + ky*kxl*nl2 + kx*nl2 + j*nl_region + i;
                    for (int k = 0; k < 6; k++){
                        v1 = point_local[index_value1*6 + k];
                        if(v1== halo[fc*kxl*kxl*nh + ky*kxl*nh + kx*nh + i - 1])ind_nei = k;
                    }
                    if(ind_nei > 0) halo[fc*kxl*kxl*nh + ky*kxl*nh + kx*nh + i] = point_local[index_value1*6 + ind_nei-1];
                    else            halo[fc*kxl*kxl*nh + ky*kxl*nh + kx*nh + i] = point_local[index_value1*6 + 5];

                    //Side 2
                    i = nl_region -1;
                    j = 0           ;
                    index_value1 = fc*nfaces*nl2 + ky*kxl*nl2 + kx*nl2 + j*nl_region + i;
                    for (int k = 0; k < 6; k++){
                        v1 = point_local[index_value1*6 + k];
                        if(v1== halo[fc*kxl*kxl*nh + ky*kxl*nh + kx*nh + nl_region - 1])ind_nei = k;
                    }
                    if(ind_nei > 0) halo[fc*kxl*kxl*nh + ky*kxl*nh + kx*nh + nl_region +2] = point_local[index_value1*6 + ind_nei-1];
                    else            halo[fc*kxl*kxl*nh + ky*kxl*nh + kx*nh + nl_region +2] = point_local[index_value1*6 + 5];

                    i = nl_region-1;
                    for (int j = 0; j < nl_region-1; j++){
                         index_value1 = fc*nfaces*nl2 + ky*kxl*nl2 + kx*nl2 + j*nl_region + i;
                         index_value2 = fc*nfaces*nl2 + ky*kxl*nl2 + kx*nl2 + (j+1)*nl_region + i;
                         for (int k = 0; k < 6; k++){
                             for (int kk = 0; kk < 6; kk++){
                                 v1 = point_local[index_value1*6 + k];
                                 v2 = point_local[index_value2*6 + kk];
                                 if(v1==v2){
                                     out_ind = 0;
                                     for (int kkk = 0; kkk < nl2; kkk++)    if(v1==fc*nfaces*nl2 + ky*kxl*nl2 + kx*nl2 + kkk) out_ind = 1;
                                     if(out_ind == 0)halo[fc*kxl*kxl*nh + ky*kxl*nh + kx*nh + (nl_region+2) + j+1] = v1;
                                 }
                             }
                         }
                    }

                    // Side 3

                    i = nl_region -1;
                    j = nl_region -1;
                    index_value1 = fc*nfaces*nl2 + ky*kxl*nl2 + kx*nl2 + j*nl_region + i;
                    for (int k = 0; k < 6; k++){
                        v1 = point_local[index_value1*6 + k];
                        if(v1 == halo[fc*kxl*kxl*nh + ky*kxl*nh + kx*nh + nl_region + 2 + nl_region - 1])ind_nei = k;
                    }
                    if(ind_nei > 0) halo[fc*kxl*kxl*nh + ky*kxl*nh + kx*nh + 2*(nl_region +2)] = point_local[index_value1*6 + ind_nei-1];
                    else            halo[fc*kxl*kxl*nh + ky*kxl*nh + kx*nh + 2*(nl_region +2)] = point_local[index_value1*6 + 5];

                    for (int i = 0; i < nl_region-1; i++){
                         j = nl_region-1;
                         index_value1 = fc*nfaces*nl2 + ky*kxl*nl2 + kx*nl2 + j*nl_region + nl_region - i-1;
                         index_value2 = fc*nfaces*nl2 + ky*kxl*nl2 + kx*nl2 + j*nl_region + nl_region - i-2;
                         for (int k = 0; k < 6; k++){
                             for (int kk = 0; kk < 6; kk++){
                                 v1 = point_local[index_value1*6 + k];
                                 v2 = point_local[index_value2*6 + kk];
                                 if(v1==v2){
                                     out_ind = 0;
                                     for (int kkk = 0; kkk < nl2; kkk++) if(v1==fc*nfaces*nl2 + ky*kxl*nl2 + kx*nl2 + kkk) out_ind = 1;
                                     if(out_ind == 0) halo[fc*kxl*kxl*nh + ky*kxl*nh + kx*nh + 2*(nl_region+2) + i+1] = v1;
                                 }
                             }
                         }
                    }

                    i = 0;
                    j = nl_region -1;
                    index_value1 = fc*nfaces*nl2 + ky*kxl*nl2 + kx*nl2 + j*nl_region + i;
                    for (int k = 0; k < 6; k++){
                        v1 = point_local[index_value1*6 + k];
                        if(v1 == halo[fc*kxl*kxl*nh + ky*kxl*nh + kx*nh + 2*(nl_region + 2) + nl_region - 1]) ind_nei = k;
                    }
                    if(ind_nei > 0) halo[fc*kxl*kxl*nh + ky*kxl*nh + kx*nh + 2*(nl_region +2) +nl_region] = point_local[index_value1*6 + ind_nei-1];
                    else            halo[fc*kxl*kxl*nh + ky*kxl*nh + kx*nh + 2*(nl_region +2) +nl_region] = point_local[index_value1*6 + 5];

                    // Side 4
                    i = 0;
                    j = nl_region -1;
                    index_value1 = fc*nfaces*nl2 + ky*kxl*nl2 + kx*nl2 + j*nl_region + i;
                    for (int k = 0; k < 6; k++){
                        v1 = point_local[index_value1*6 + k];
                        if(v1 == halo[fc*kxl*kxl*nh + ky*kxl*nh + kx*nh + 2*(nl_region + 2) + nl_region ]) ind_nei = k;
                    }
                    if(ind_nei > 0) halo[fc*kxl*kxl*nh + ky*kxl*nh + kx*nh + 3*(nl_region +2)] = point_local[index_value1*6 + ind_nei-1];
                    else            halo[fc*kxl*kxl*nh + ky*kxl*nh + kx*nh + 3*(nl_region +2)] = point_local[index_value1*6 + 5];

                    i = 0;
                    for (int j = 0; j < nl_region-1; j++){
                         index_value1 = fc*nfaces*nl2 + ky*kxl*nl2 + kx*nl2 + (nl_region - j-1)*nl_region;
                         index_value2 = fc*nfaces*nl2 + ky*kxl*nl2 + kx*nl2 + (nl_region - j-2)*nl_region;
                         for (int k = 0; k < 6; k++){
                             for (int kk = 0; kk < 6; kk++){
                                 v1 = point_local[index_value1*6 + k];
                                 v2 = point_local[index_value2*6 + kk];
                                 if(v1==v2){
                                     out_ind = 0;
                                     for (int kkk = 0; kkk < nl2; kkk++) if(v1==fc*nfaces*nl2 + ky*kxl*nl2 + kx*nl2 + kkk) out_ind = 1;
                                     if(out_ind == 0)halo[fc*kxl*kxl*nh + ky*kxl*nh + kx*nh + 3*(nl_region+2) + j+1] = v1;
                                 }
                             }
                         }
                    }

                    i = 0;
                    j = 0;
                    index_value1 = fc*nfaces*nl2 + ky*kxl*nl2 + kx*nl2 + j*nl_region + i;
                    for (int k = 0; k < 6; k++){
                        v1 = point_local[index_value1*6 + k];
                        if(v1 == halo[fc*kxl*kxl*nh + ky*kxl*nh + kx*nh + 3*(nl_region + 2) + nl_region -1]) ind_nei = k;
                    }
                    if(ind_nei > 0) halo[fc*kxl*kxl*nh + ky*kxl*nh + kx*nh + 3*(nl_region +2)+nl_region] = point_local[index_value1*6 + ind_nei-1];
                    else            halo[fc*kxl*kxl*nh + ky*kxl*nh + kx*nh + 3*(nl_region +2)+nl_region] = point_local[index_value1*6 + 5];

                }
            }
        }
    }
}

void Icogrid::reorder_neighbors_indx_rhombi (int *point_local,
                                             int *halo       ,
                                             int *pent_ind   ,
                                             int  nl_region  ,
                                             int  nl2        ,
                                             int  ni         ,
                                             int  nlhalo     ,
                                             int  point_num  ){

//
//  Description:
//
//
//  Input:  - halo        - Indexes of the points that form the halos.
//          - pent_ind    - Pentagon's indexes
//          - nl_region   - Number of points at each side of the rhombi.
//          - nl2         - Number of points in a rhombi.
//          - point_num   - Number of points.
//          - nlhalo      - Length of the halo side.
//
//  Output: - point_local - First neighbors.
//

//  Local variables
    int geo, np;
    int nh2 = nlhalo*4;
    for (int i = 0; i < point_num; i++)    for (int j = 0; j < 6; j++) point_local[i*6 + j] = -1;

    for (int i = 0; i < ni; i++){

        // Corner
        geo = 6; // Hexagons.
        np = nl2*i;
        for (int l = 0; l <12; l++) if (np == pent_ind[l]) geo = 5; // Pentagons.

        // Find neibours.
        for (int j = 0; j < nl_region; j++){
            for (int k = 0; k < nl_region; k++){
                np  = nl2*i + j*nl_region + k;
                // Right
                if(j==0){
                    point_local[np*6 + 0] = nl2*i + (j+1)*nl_region + k;
                    point_local[np * 6 + 3] = halo[i*nh2 + k];
                }
                else if(j==nl_region-1){
                    point_local[np*6 + 0] = halo[i*nh2 + 2*nlhalo+nl_region-k];
                    point_local[np * 6 + 3] = nl2*i + (j-1)*nl_region + k;
                }
                else{
                    point_local[np*6 + 0] = nl2*i + (j+1)*nl_region + k;
                    point_local[np*6 + 3] = nl2*i + (j-1)*nl_region + k;
                }
                // Up
                if(k == 0){
                    point_local[np*6 + 2] = nl2*i + j*nl_region + k+1;
                    if (geo == 6) point_local[np * 6 + 5] = halo[i*nh2 + 3 * nlhalo + nl_region - j - 1];
                    else{
                        if (j == 0) point_local[np * 6 + 5] = -1;
                        else  point_local[np * 6 + 5] = halo[i*nh2 + 3 * nlhalo + nl_region - j - 1];
                    }
                }
                else if(k == nl_region-1){
                    point_local[np*6 + 2] = halo[i*nh2 + nlhalo+j];
                    point_local[np * 6 + 5] = nl2*i + j*nl_region + k - 1;
                }
                else{
                    point_local[np*6 + 2] = nl2*i + j*nl_region + k+1;
                    point_local[np * 6 + 5] = nl2*i + j*nl_region + k - 1;
                }
                // Front
                if(j==0){
                    if(k == 0){
                        if (geo == 6)point_local[np * 6 + 4] = halo[i*nh2 + 3 * nlhalo + nl_region];
                        else point_local[np * 6 + 4] = halo[i*nh2 + 3 * nlhalo + nl_region-1];
                        point_local[np*6 + 1] = nl2*i + nl_region + k+1;
                    }
                    else if(k == nl_region -1){
                        point_local[np*6 + 4] = halo[i*nh2 + nl_region-2];
                        point_local[np*6 + 1] = halo[i*nh2 + nl_region+3];
                    }
                    else{
                        point_local[np*6 + 4] = halo[i*nh2 + k-1];
                        point_local[np*6 + 1] = nl2*i + (j+1)*nl_region + k+1;
                    }
                }
                else if(j == nl_region-1){
                    if(k == 0){
                        point_local[np*6 + 4] = halo[i*nh2 + 3*nlhalo+1];
                        point_local[np * 6 + 1] = halo[i*nh2 + 2 * nlhalo + nl_region - 1];
                    }
                    else if(k == nl_region -1){
                        point_local[np*6 + 4] = nl2*i + (j-1)*nl_region + k-1;
                        point_local[np*6 + 1] = halo[i*nh2 + 2*nlhalo];
                    }
                    else{
                        point_local[np*6 + 4] = nl2*i + (j-1)*nl_region + k-1;
                        point_local[np*6 + 1] = halo[i*nh2 + 2*nlhalo-k+nl_region-1];
                    }
                }
                else{
                    if(k == 0){
                        point_local[np*6 + 4] = halo[i*nh2 + 3*nlhalo-j+nl_region];
                        point_local[np*6 + 1] = nl2*i + (j+1)*nl_region + k+1;
                    }
                    else if(k == nl_region -1){
                        point_local[np*6 + 4] = nl2*i + (j-1)*nl_region + k-1;
                        point_local[np*6 + 1] = halo[i*nh2 + nl_region+2+j+1];
                    }
                    else{
                        point_local[np*6 + 4] = nl2*i + (j-1)*nl_region + k-1;
                        point_local[np*6 + 1] = nl2*i + (j+1)*nl_region + k+1;
                    }
                }
            }
        }
    }
}



void Icogrid::produce_maps(int *maps     ,
                           int *halo     ,
                           int  nr       ,
                           int nl_region ){

//
//  Description:
//
//  Produces regions that contain the halos.
//
//  Input:  - halo         - halos of the main rhombi.
//          - nr           - number of the rhombi.
//          - nl_region    - length of the rhombus side.
//
//  Output: - maps - indexes of the maps that contain the main rhombi and halos.
//

    // Local variables
    int nl2 = nl_region + 2;
    int nl22 = nl2*nl2;
    int nh2 = 4 * (nl_region + 2);
    int nl_region2 = nl_region*nl_region;

    for (int l = 0; l < nr; l++){
        for (int i = 0; i < nl_region; i++){
            for (int j = 0; j < nl_region; j++){
                maps[l*nl22 + (j+1)*nl2 + i+1] = l*nl_region2 + j*nl_region + i;
            }
        }
        //Side 1
        for (int i = 0; i < nl_region+1; i++){
            maps[l*nl22 + 0 * (nl_region + 2) + i + 1] = halo[l*nh2 + i];
        }
        //Side 2
        for (int i = 0; i < nl_region + 2; i++){
            maps[l*nl22 + i * (nl_region + 2) + 0] = halo[l*nh2 + 3 * nl2 + nl_region - i];
        }
        //Side 3
        for (int i = 0; i < nl_region + 1; i++){
            maps[l*nl22 + (nl_region + 1) * (nl_region + 2) + i+1] = halo[l*nh2 + 2 * nl2 + nl_region - i];
        }
        //Side 4
        for (int i = 0; i < nl_region; i++){
            maps[l*nl22 + (i + 1) * (nl_region + 2) + nl_region + 1] = halo[l*nh2 + nl2 + i];
        }
    }
}

void Icogrid::spring_dynamics(int *point_local    ,
                              int *pent_ind       ,
                              int  glevel         ,
                              double spring_beta  ,
                              double *xyz         ,
                              int point_num       ){

//
//  Description:
//
//  Method developed in Tomita et al. 2002.
//
//  Input:  - xyz         - Standard icosahedral points;
//          - point_local - First neighbors.
//          - pent_ind    - Pentagon's indexes
//          - spring_beta - Natural sring constant.
//          - point_num   - Number of points.
//
//  Output: - xyz - New vertex's positions
//

    //  Local variables
    double Wx, Wy, Wz;
    double Vx, Vy, Vz;
    double velx_val, vely_val, velz_val;
    double Fxnet, Fynet, Fznet;
    double X, Y, Z;
    double l, d, H, max_v;

    //Local arrays.
    double *xyzi;
    double *velx, *vely, *velz;
    double *Px, *Py, *Pz;
    double *Fx, *Fy, *Fz;
    double *Px_Nei, *Py_Nei, *Pz_Nei;
    double *Fx_Nei, *Fy_Nei, *Fz_Nei;
    double *Fnet;

    // Central point
    xyzi = new double[12*3]();
    velx = new double[point_num]();
    vely = new double[point_num]();
    velz = new double[point_num]();
    Fnet = new double[point_num]();
    Fx   = new double[point_num]();
    Fy   = new double[point_num]();
    Fz   = new double[point_num]();
    Px   = new double[point_num]();
    Py   = new double[point_num]();
    Pz   = new double[point_num]();

    // First neighbors.
    Fx_Nei   = new double[point_num*6]();
    Fy_Nei   = new double[point_num*6]();
    Fz_Nei   = new double[point_num*6]();
    Px_Nei   = new double[point_num*6]();
    Py_Nei   = new double[point_num*6]();
    Pz_Nei   = new double[point_num*6]();

    // Routine parameters.
    double lba       = 2.0*M_PI/(10.0*pow(2.0,glevel-1));
    double dbar      = spring_beta*lba;
    double drag      = 1;                                   // Friction coefficient.
    double tstep     = 2.0E-2;                              // time step.
    double cond      = 1.0E-4;                              // Criteria for convergence.

    int lim = 100000; // Limit of iterations.

    printf("\n\n Running spring dynamics.\n\n");

    for (int i = 0; i < point_num; i++){
        velx[i] = 0.0;
        vely[i] = 0.0;
        velz[i] = 0.0;
    }

    for (int i = 0; i < 12; i++){
        xyzi[i * 3 + 0] = xyz[pent_ind[i] * 3 + 0];
        xyzi[i * 3 + 1] = xyz[pent_ind[i] * 3 + 1];
        xyzi[i * 3 + 2] = xyz[pent_ind[i] * 3 + 2];
        point_local[pent_ind[i]*6 + 5] = point_local[pent_ind[i]*6];
    }

    // Solving spring dynamics.
    for (int it = 0; it < lim; it++){

        for (int i = 0; i < point_num; i++){
            Px[i] = xyz[i*3 + 0];
            Py[i] = xyz[i*3 + 1];
            Pz[i] = xyz[i*3 + 2];
            for (int j = 0; j < 6; j++){
                Px_Nei[i*6 + j] = xyz[point_local[i*6 + j]*3 + 0];
                Py_Nei[i*6 + j] = xyz[point_local[i*6 + j]*3 + 1];
                Pz_Nei[i*6 + j] = xyz[point_local[i*6 + j]*3 + 2];
            }
        }
        for (int i = 0; i < point_num; i++){
            for (int j = 0; j < 6; j++){

                Wx = Py[i]*Pz_Nei[i*6 + j] - Pz[i]*Py_Nei[i*6 + j];
                Wy = Pz[i]*Px_Nei[i*6 + j] - Px[i]*Pz_Nei[i*6 + j];
                Wz = Px[i]*Py_Nei[i*6 + j] - Py[i]*Px_Nei[i*6 + j];

                Vx = Wy*Pz[i] - Wz*Py[i];
                Vy = Wz*Px[i] - Wx*Pz[i];
                Vz = Wx*Py[i] - Wy*Px[i];

                //norm
                l = sqrt(Vx*Vx + Vy*Vy + Vz*Vz);

                d = acos(Px[i]*Px_Nei[i*6 + j] + Py[i]*Py_Nei[i*6 + j] + Pz[i]*Pz_Nei[i*6 + j]);

                Fx_Nei[i*6 + j] = (d - dbar) * Vx;
                Fy_Nei[i*6 + j] = (d - dbar) * Vy;
                Fz_Nei[i*6 + j] = (d - dbar) * Vz;

                Fx_Nei[i*6 + j] = Fx_Nei[i*6 + j]/l;
                Fy_Nei[i*6 + j] = Fy_Nei[i*6 + j]/l;
                Fz_Nei[i*6 + j] = Fz_Nei[i*6 + j]/l;
            }
        }

        for (int i = 0; i < 12; i++){

            Fx_Nei[pent_ind[i]*6 + 5] = 0.0;
            Fy_Nei[pent_ind[i]*6 + 5] = 0.0;
            Fz_Nei[pent_ind[i]*6 + 5] = 0.0;

        }

        for (int i = 0; i < point_num; i++){

                // Calculate net forces in each point.
                Fxnet = Fx_Nei[i*6 + 0] + Fx_Nei[i*6 + 1] + Fx_Nei[i*6 + 2] + Fx_Nei[i*6 + 3] + Fx_Nei[i*6 + 4] + Fx_Nei[i*6 + 5];
                Fynet = Fy_Nei[i*6 + 0] + Fy_Nei[i*6 + 1] + Fy_Nei[i*6 + 2] + Fy_Nei[i*6 + 3] + Fy_Nei[i*6 + 4] + Fy_Nei[i*6 + 5];
                Fznet = Fz_Nei[i*6 + 0] + Fz_Nei[i*6 + 1] + Fz_Nei[i*6 + 2] + Fz_Nei[i*6 + 3] + Fz_Nei[i*6 + 4] + Fz_Nei[i*6 + 5];

                // Used to check convergence.
                Fnet[i] = sqrt(Fxnet * Fxnet + Fynet * Fynet + Fznet * Fznet) / lba;

                // Drag move.
                Fx[i] = Fxnet - drag * velx[i];
                Fy[i] = Fynet - drag * vely[i];
                Fz[i] = Fznet - drag * velz[i];
        }

        for (int i = 0; i < point_num; i++){

            // Update points
            X = xyz[i*3 + 0] + velx[i]*tstep;
            Y = xyz[i*3 + 1] + vely[i]*tstep;
            Z = xyz[i*3 + 2] + velz[i]*tstep;

            l = sqrt(X*X + Y*Y + Z*Z);

            xyz[i*3 + 0] = X/l;
            xyz[i*3 + 1] = Y/l;
            xyz[i*3 + 2] = Z/l;

        }

        for (int i = 0; i < point_num; i++){

            // Update vel
            velx_val = velx[i] + Fx[i]*tstep;
            vely_val = vely[i] + Fy[i]*tstep;
            velz_val = velz[i] + Fz[i]*tstep;

            H = xyz[i*3 + 0]*velx_val + xyz[i*3 + 1]*vely_val + xyz[i*3 + 2]*velz_val;

            // Remove radial component (if any).
            velx[i] = velx_val - H*xyz[i*3 + 0];
            vely[i] = vely_val - H*xyz[i*3 + 1];
            velz[i] = velz_val - H*xyz[i*3 + 2];

        }

        // Fix Petagon's position.
        for (int i = 0; i < 12; i++){

            velx[pent_ind[i]]      = 0.0;
            vely[pent_ind[i]]      = 0.0;
            velz[pent_ind[i]]      = 0.0;
            Fnet[pent_ind[i]]      = 0.0;
            xyz[pent_ind[i] * 3 + 0] = xyzi[i * 3 + 0];
            xyz[pent_ind[i] * 3 + 1] = xyzi[i * 3 + 1];
            xyz[pent_ind[i] * 3 + 2] = xyzi[i * 3 + 2];

        }

        // Check convergence.
        max_v = -1E-10;
        for (int i = 0; i < point_num; i++)    if( Fnet[i] > max_v ) max_v = Fnet[i];
        if (max_v < cond) break;

    }

    printf(" Done!\n\n");

    delete [] xyzi;
    delete [] velx;
    delete [] vely;
    delete [] velz;
    delete [] Fx;
    delete [] Fy;
    delete [] Fz;
    delete [] Px;
    delete [] Py;
    delete [] Pz;
    delete [] Fnet;
    delete [] Fx_Nei;
    delete [] Fy_Nei;
    delete [] Fz_Nei;
    delete [] Px_Nei;
    delete [] Py_Nei;
    delete [] Pz_Nei;

}

void Icogrid::find_qpoints (int    *point_local,
                            double *xyzq       ,
                            double *xyz        ,
                            int    *pent_ind   ,
                            int point_num      ){

//
//  Description:
//
//  Finds the coordinates of the q-points of each vertex.
//
//  Input: - pointlocal - First neighbors indexes of each vertex.
//         - xyz        - Vertices coordinates.
//         - point_num  - Number of vertices.
//         - pen_ind    - Pentagon's indexes.
//
//  Output:- xyzq       - Q-points coordinates.
//

    double3 * xyz3 = (double3*)xyz;
    double3 * xyzq3 = (double3*)xyzq;
    
    double3 vc1, vc2, vc3;

    int geo;
    
    for (int i = 0; i < point_num; i++){
        geo = 6; // Hexagons.
        for (int k = 0; k < 12; k++) if(i == pent_ind[k]) geo = 5; // Pentagons.
        if(geo == 5){
            for (int j = 0; j < 4; j++){
                vc1 = normproj( xyz3[point_local[i*6 + j]], xyz3[i] );
                vc2 = normproj( xyz3[point_local[i*6 + j+1]], xyz3[point_local[i*6 + j]] );
                vc3 = normproj( xyz3[i], xyz3[point_local[i*6 + j+1]] );
                xyzq3[i*6 + j] = normalize(vc1 + vc2 + vc3);
            }

            vc1 = normproj( xyz3[point_local[i*6 + 4]], xyz3[i] );
            vc2 = normproj( xyz3[point_local[i*6 + 0]], xyz3[point_local[i*6 + 4]] );
            vc3 = normproj( xyz3[i], xyz3[point_local[i*6 + 0]] );

            xyzq3[i*6 + 4] = normalize(vc1 + vc2 + vc3);
        }
        else{   // Hexagons
            for (int j = 0; j < 5; j++){
                vc1 = normproj( xyz3[point_local[i*6 + j]], xyz3[i] );
                vc2 = normproj( xyz3[point_local[i*6 + j+1]], xyz3[point_local[i*6 + j]] );
                vc3 = normproj( xyz3[i], xyz3[point_local[i*6 + j+1]] );
                
                xyzq3[i*6 + j] = normalize(vc1 + vc2 + vc3);
            }
            vc1 = normproj( xyz3[point_local[i*6 + 5]], xyz3[i] );
            vc2 = normproj( xyz3[point_local[i*6 + 0]], xyz3[point_local[i*6 + 5]] );
            vc3 = normproj( xyz3[i], xyz3[point_local[i*6 + 0]] );
            
            xyzq3[i*6 + 5] = normalize(vc1 + vc2 + vc3);
        }
    }
}


void Icogrid::relocate_centres(int    *point_local,
                               double *xyzq       ,
                               double *xyz        ,
                               int    *pent_ind   ,
                               int point_num      ){

//
//  Description:
//
//  Corrects the centroids positions.
//
//  Input: - point_local - First neighbors indexes.
//         - xyz         - Centers of the control volumes.
//         - xyzq        - q-points coordinates.
//         - pent_ind    - Indexes of the pentagons.
//         - point_num   - Number of vertices.
//
//  Output: - xyz - New centroids.
//
//

    // Local variables.
    int geo;

    // Local arrays-
    double3 * xyzq3 = (double3*)xyzq;
    double3 * xyz3 = (double3*)xyz;

    // local vectors
    double3 vc1, vc2, vc3, vc4, vc5, vc6;
    double3 vgc;

    for (int i = 0; i < point_num; i++){
        geo = 6; // Hexagons.
        for (int k = 0; k < 12; k++) if(i == pent_ind[k]) geo = 5; // Pentagons.
        if(geo == 5){


            vc1 = normproj(xyzq3[i*6 + 1], xyzq3[i*6 + 0]);
            vc2 = normproj(xyzq3[i*6 + 2], xyzq3[i*6 + 1]);
            vc3 = normproj(xyzq3[i*6 + 3], xyzq3[i*6 + 2]);
            vc4 = normproj(xyzq3[i*6 + 4], xyzq3[i*6 + 3]);
            vc5 = normproj(xyzq3[i*6 + 0], xyzq3[i*6 + 4]);
            
            vgc = vc1 + vc2 + vc3 + vc4 + vc5;

            xyz3[i] = normalize(vgc);
        }
        else{
            vc1 = normproj(xyzq3[i*6 + 1], xyzq3[i*6 + 0]);
            vc2 = normproj(xyzq3[i*6 + 2], xyzq3[i*6 + 1]);
            vc3 = normproj(xyzq3[i*6 + 3], xyzq3[i*6 + 2]);
            vc4 = normproj(xyzq3[i*6 + 4], xyzq3[i*6 + 3]);
            vc5 = normproj(xyzq3[i*6 + 5], xyzq3[i*6 + 4]);
            vc6 = normproj(xyzq3[i*6 + 0], xyzq3[i*6 + 5]);
            
            vgc = vc1 + vc2 + vc3 + vc4 + vc5 + vc6;

            xyz3[i] = normalize(vgc);
        }
    }
}

void Icogrid::set_altitudes(double *Altitude    ,
                            double *Altitudeh   ,
                            double  Top_altitude,
                            int nv              ){

//
//  Description:
//
//  Sets the layers and interfaces altitudes.
//
//  Input:  - nv - Number of vertical layers.
//          - top_altitude - Altitude of the top model domain.
//
//  Output: - Altitude  - Layers altitudes.
//          - Altitudeh - Interfaces altitudes.
//

    double res_vert = Top_altitude / nv;
    Altitudeh[0] = 0.0;
    for (int lev = 0; lev < nv; lev++)  Altitudeh[lev + 1] = Altitudeh[lev] + res_vert;
    for (int lev = 0; lev < nv; lev++)  Altitude[lev] = (Altitudeh[lev] + Altitudeh[lev+1])/2.0;
}


void Icogrid::cart2sphe ( double *lonlat   ,
                          double *xyz      ,
                          int point_num    ){

//
//  Description:
//
//  Transform Cartesian to spherical coordinates.
//
//  Input: xyz       - Vertices coordinates.
//         point_num - Number of vertices.
//
//  Output: - lonlat - Spherical coordinates.
//
    for (int i = 0; i < point_num; i++){
        lonlat[i*2 + 0] = atan2(xyz[i*3 + 1], xyz[i*3 + 0]);
        lonlat[i*2 + 1] = atan2(xyz[i*3 + 2], sqrt(pow(xyz[i*3 + 0],2) + pow(xyz[i*3 + 1],2)));
    }
}

void Icogrid::correct_xyz_points ( double A     ,
                                   double *xyzq ,
                                   double *xyz  ,
                                   int *pent_ind,
                                   int point_num){

//
//  Description:
//
//  Corrects the points positions taking into account the radius of the planet.
//
//  Input: - A - Planet radius.
//         - xyz         - Centers of the control volumes.
//         - xyzq        - q-points coordinates.
//         - pent_ind    - Indexes of the pentagons.
//         - point_num   - Number of vertices.
//
//
//  Output:- xyz         - Centers of the control volumes.
//
    // Local variables.
    int geo;

    double3 * xyzq3 = (double3*)xyzq;
    double3 * xyz3 = (double3*)xyz;    
    
    for (int i = 0; i < point_num; i++){
        geo = 6; // Hexagons.
        for (int k = 0; k < 12; k++) if(i == pent_ind[k]) geo = 5; // Pentagons.
        if(geo == 5){
            for (int j = 0; j < 5; j++){    // Pentagons.
                xyzq3[i*6 + j] = xyzq3[i*6 + j]*A;
            }
            xyz3[i] = xyz3[i]*A;
        }
        else{
            for (int j = 0; j < 6; j++){    // Hexagons.
                xyzq3[i*6 + j] = xyzq3[i*6 + j + 0]*A;
            }
            xyz3[i] = xyz3[i]*A;
        }
    }
}


void Icogrid::control_areas ( double *areasT  ,
                              double *areasTr ,
                              double *areasq  ,
                              int *point_local,
                              double *xyzq    ,
                              double *xyz     ,
                              int *pent_ind   ,
                              int point_num   ){

//
//  Description:
//
//  Compute control areas.
//
//  Input: - xyz         - Vertices coordinates.
//         - xyzq        - Q-points coordinates.
//         - point_local - First neighbors indexes of each vertex.
//         - pent_ind    - Pentagon' indexes.
//         - point_num   - Number of vertices.
//
//  Output: areas      - Sub-areas (alpha, beta and gamma).
//          areasT     - Control volume area.
//

    // Local variables.
    int geo;
    double t11;
    double t22;
    double t33;
    double fac11, fac21, fac12, fac22, fac13, fac23;

    double radius;

    // Local arrays
    double *a, *b, *c;
    double *v01, *v02, *v03;
    double *v11, *v12, *v13;
    double *v21, *v22, *v23;
    double *w11, *w21, *w12, *w22;
    double *w13, *w23;
    double *ang, areav;

    a   = new double[3]();
    b   = new double[3]();
    c   = new double[3]();
    v01 = new double[3]();
    v02 = new double[3]();
    v03 = new double[3]();
    v11 = new double[3]();
    v12 = new double[3]();
    v13 = new double[3]();
    v21 = new double[3]();
    v22 = new double[3]();
    v23 = new double[3]();
    w11 = new double[3]();
    w21 = new double[3]();
    w12 = new double[3]();
    w22 = new double[3]();
    w13 = new double[3]();
    w23 = new double[3]();
    ang = new double[3]();

    for (int i = 0; i < point_num*6*3; i++) areas[i] = 0.0;

    for (int i = 0; i < point_num; i++){
//
// Areas - Alpha, Beta, Gamma
//
        geo = 6; // Hexagons.
        for (int k = 0; k < 12; k++) if(i == pent_ind[k]) geo = 5; // Pentagons.
        if(geo == 5){
            for (int j = 0; j < 5; j++){    //Pentagons
                for (int k = 0; k < 3; k++){
                    if ( j < 4 ){
                        if (k == 0){
                            a[0] = xyz[point_local[i*6 + j+1]*3 + 0];
                            a[1] = xyz[point_local[i*6 + j+1]*3 + 1];
                            a[2] = xyz[point_local[i*6 + j+1]*3 + 2];
                            b[0] = xyz[point_local[i*6 + j]*3 + 0];
                            b[1] = xyz[point_local[i*6 + j]*3 + 1];
                            b[2] = xyz[point_local[i*6 + j]*3 + 2];
                            c[0] = xyzq[i*6*3 + j*3 + 0];
                            c[1] = xyzq[i*6*3 + j*3 + 1];
                            c[2] = xyzq[i*6*3 + j*3 + 2];
                        }
                        else if ( k ==1 ){
                            a[0] = xyz[point_local[i*6 + j+1]*3 + 0];
                            a[1] = xyz[point_local[i*6 + j+1]*3 + 1];
                            a[2] = xyz[point_local[i*6 + j+1]*3 + 2];
                            b[0] = xyz[i*3 + 0];
                            b[1] = xyz[i*3 + 1];
                            b[2] = xyz[i*3 + 2];
                            c[0] = xyzq[i*6*3 + j*3 + 0];
                            c[1] = xyzq[i*6*3 + j*3 + 1];
                            c[2] = xyzq[i*6*3 + j*3 + 2];
                        }
                        else{
                            a[0] = xyz[i*3 + 0];
                            a[1] = xyz[i*3 + 1];
                            a[2] = xyz[i*3 + 2];
                            b[0] = xyz[point_local[i*6 + j]*3 + 0];
                            b[1] = xyz[point_local[i*6 + j]*3 + 1];
                            b[2] = xyz[point_local[i*6 + j]*3 + 2];
                            c[0] = xyzq[i*6*3 + j*3 + 0];
                            c[1] = xyzq[i*6*3 + j*3 + 1];
                            c[2] = xyzq[i*6*3 + j*3 + 2];
                        }
                    }
                    else{
                        if (k == 0){
                            a[0] = xyz[point_local[i*6 + 0]*3 + 0];
                            a[1] = xyz[point_local[i*6 + 0]*3 + 1];
                            a[2] = xyz[point_local[i*6 + 0]*3 + 2];
                            b[0] = xyz[point_local[i*6 + 4]*3 + 0];
                            b[1] = xyz[point_local[i*6 + 4]*3 + 1];
                            b[2] = xyz[point_local[i*6 + 4]*3 + 2];
                            c[0] = xyzq[i*6*3 + 4*3 + 0];
                            c[1] = xyzq[i*6*3 + 4*3 + 1];
                            c[2] = xyzq[i*6*3 + 4*3 + 2];
                        }
                        else if ( k ==1 ){
                            a[0] = xyz[point_local[i*6 + 0]*3 + 0];
                            a[1] = xyz[point_local[i*6 + 0]*3 + 1];
                            a[2] = xyz[point_local[i*6 + 0]*3 + 2];
                            b[0] = xyz[i*3 + 0];
                            b[1] = xyz[i*3 + 1];
                            b[2] = xyz[i*3 + 2];
                            c[0] = xyzq[i*6*3 + 4*3 + 0];
                            c[1] = xyzq[i*6*3 + 4*3 + 1];
                            c[2] = xyzq[i*6*3 + 4*3 + 2];
                        }
                        else{
                            a[0] = xyz[i*3 + 0];
                            a[1] = xyz[i*3 + 1];
                            a[2] = xyz[i*3 + 2];
                            b[0] = xyz[point_local[i*6 + 4]*3 + 0];
                            b[1] = xyz[point_local[i*6 + 4]*3 + 1];
                            b[2] = xyz[point_local[i*6 + 4]*3 + 2];
                            c[0] = xyzq[i*6*3 + 4*3 + 0];
                            c[1] = xyzq[i*6*3 + 4*3 + 1];
                            c[2] = xyzq[i*6*3 + 4*3 + 2];
                        }
                    }
                    for (int v = 0; v < 3; v++){
                        v01[v] = a[v];
                        v02[v] = b[v];
                        v03[v] = c[v];

                        v11[v] = b[v] - a[v];
                        v12[v] = a[v] - b[v];
                        v13[v] = a[v] - c[v];

                        v21[v] = c[v] - a[v];
                        v22[v] = c[v] - b[v];
                        v23[v] = b[v] - c[v];
                    }

                    t11 = 1/(v01[0]*v01[0] + v01[1]*v01[1] + v01[2]*v01[2]);

                    fac11       = (v01[0]*v11[0] + v01[1]*v11[1] + v01[2]*v11[2])*t11;
                    fac21       = (v01[0]*v21[0] + v01[1]*v21[1] + v01[2]*v21[2])*t11;

                    w11[0] = v11[0] - fac11*v01[0];
                    w11[1] = v11[1] - fac11*v01[1];
                    w11[2] = v11[2] - fac11*v01[2];

                    w21[0] = v21[0] - fac21*v01[0];
                    w21[1] = v21[1] - fac21*v01[1];
                    w21[2] = v21[2] - fac21*v01[2];

                    ang[0] = (w11[0]*w21[0] + w11[1]*w21[1] + w11[2]*w21[2])/
                        (sqrt(w11[0]*w11[0] + w11[1]*w11[1] + w11[2]*w11[2])*
                         sqrt(w21[0]*w21[0] + w21[1]*w21[1] + w21[2]*w21[2]));

                    if (ang[0] > 1) ang[0] =  1;
                    if (ang[0] < -1)ang[0] = -1;

                    ang[0] = acos(ang[0]);

                    t22 = 1/(v02[0]*v02[0] + v02[1]*v02[1] + v02[2]*v02[2]);

                    fac12       = (v02[0]*v12[0] + v02[1]*v12[1] + v02[2]*v12[2])*t22;
                    fac22       = (v02[0]*v22[0] + v02[1]*v22[1] + v02[2]*v22[2])*t22;

                    w12[0] = v12[0] - fac12*v02[0];
                    w12[1] = v12[1] - fac12*v02[1];
                    w12[2] = v12[2] - fac12*v02[2];

                    w22[0] = v22[0] - fac22*v02[0];
                    w22[1] = v22[1] - fac22*v02[1];
                    w22[2] = v22[2] - fac22*v02[2];

                    ang[1] = (w12[0]*w22[0] + w12[1]*w22[1] + w12[2]*w22[2])/
                        (sqrt(w12[0]*w12[0] + w12[1]*w12[1] + w12[2]*w12[2])*
                         sqrt(w22[0]*w22[0] + w22[1]*w22[1] + w22[2]*w22[2]));

                    if ( ang[1] > 1) ang[1] =  1;
                    if ( ang[1] < -1)ang[1] = -1;

                    ang[1] = acos(ang[1]);

                    t33 = 1/(v03[0]*v03[0] + v03[1]*v03[1] + v03[2]*v03[2]);

                    fac13       = (v03[0]*v13[0] + v03[1]*v13[1] + v03[2]*v13[2])*t33;
                    fac23       = (v03[0]*v23[0] + v03[1]*v23[1] + v03[2]*v23[2])*t33;

                    w13[0] = v13[0] - fac13*v03[0];
                    w13[1] = v13[1] - fac13*v03[1];
                    w13[2] = v13[2] - fac13*v03[2];

                    w23[0] = v23[0] - fac23*v03[0];
                    w23[1] = v23[1] - fac23*v03[1];
                    w23[2] = v23[2] - fac23*v03[2];

                    ang[2] = (w13[0]*w23[0] + w13[1]*w23[1] + w13[2]*w23[2])/
                        (sqrt(w13[0]*w13[0] + w13[1]*w13[1] + w13[2]*w13[2])*
                         sqrt(w23[0]*w23[0] + w23[1]*w23[1] + w23[2]*w23[2]));

                    if (ang[2] > 1)    ang[2] =  1;
                    if (ang[2] < -1)ang[2] = -1;

                    ang[2] = acos(ang[2]);

                    radius = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]) ;
                    areasq[i*6*3 + j*3 + k] = (ang[0] + ang[1] + ang[2] - M_PI)*pow(radius,2);
                }
                areasTr[i * 6 + j] = areasq[i * 6 * 3 + j * 3 + 0] + areasq[i * 6 * 3 + j * 3 + 1] + areasq[i * 6 * 3 + j * 3 + 2];
            }
            areasTr[i * 6 + 5] = 0.0;
        }
        else{
            for (int j = 0; j < 6; j++){    //Hexagons
                for (int k = 0; k < 3; k++){
                    if ( j < 5 ){
                        if (k == 0){
                            a[0] = xyz[point_local[i*6 + j+1]*3 + 0];
                            a[1] = xyz[point_local[i*6 + j+1]*3 + 1];
                            a[2] = xyz[point_local[i*6 + j+1]*3 + 2];
                            b[0] = xyz[point_local[i*6 + j]*3 + 0];
                            b[1] = xyz[point_local[i*6 + j]*3 + 1];
                            b[2] = xyz[point_local[i*6 + j]*3 + 2];
                            c[0] = xyzq[i*6*3 + j*3 + 0];
                            c[1] = xyzq[i*6*3 + j*3 + 1];
                            c[2] = xyzq[i*6*3 + j*3 + 2];
                        }
                        else if ( k ==1 ){
                            a[0] = xyz[point_local[i*6 + j+1]*3 + 0];
                            a[1] = xyz[point_local[i*6 + j+1]*3 + 1];
                            a[2] = xyz[point_local[i*6 + j+1]*3 + 2];
                            b[0] = xyz[i*3 + 0];
                            b[1] = xyz[i*3 + 1];
                            b[2] = xyz[i*3 + 2];
                            c[0] = xyzq[i*6*3 + j*3 + 0];
                            c[1] = xyzq[i*6*3 + j*3 + 1];
                            c[2] = xyzq[i*6*3 + j*3 + 2];
                        }
                        else{
                            a[0] = xyz[i*3 + 0];
                            a[1] = xyz[i*3 + 1];
                            a[2] = xyz[i*3 + 2];
                            b[0] = xyz[point_local[i*6 + j]*3 + 0];
                            b[1] = xyz[point_local[i*6 + j]*3 + 1];
                            b[2] = xyz[point_local[i*6 + j]*3 + 2];
                            c[0] = xyzq[i*6*3 + j*3 + 0];
                            c[1] = xyzq[i*6*3 + j*3 + 1];
                            c[2] = xyzq[i*6*3 + j*3 + 2];
                        }
                    }
                    else{
                        if (k == 0){
                            a[0] = xyz[point_local[i*6 + 0]*3 + 0];
                            a[1] = xyz[point_local[i*6 + 0]*3 + 1];
                            a[2] = xyz[point_local[i*6 + 0]*3 + 2];
                            b[0] = xyz[point_local[i*6 + 5]*3 + 0];
                            b[1] = xyz[point_local[i*6 + 5]*3 + 1];
                            b[2] = xyz[point_local[i*6 + 5]*3 + 2];
                            c[0] = xyzq[i*6*3 + 5*3 + 0];
                            c[1] = xyzq[i*6*3 + 5*3 + 1];
                            c[2] = xyzq[i*6*3 + 5*3 + 2];
                        }
                        else if ( k ==1 ){
                            a[0] = xyz[point_local[i*6 + 0]*3 + 0];
                            a[1] = xyz[point_local[i*6 + 0]*3 + 1];
                            a[2] = xyz[point_local[i*6 + 0]*3 + 2];
                            b[0] = xyz[i*3 + 0];
                            b[1] = xyz[i*3 + 1];
                            b[2] = xyz[i*3 + 2];
                            c[0] = xyzq[i*6*3 + 5*3 + 0];
                            c[1] = xyzq[i*6*3 + 5*3 + 1];
                            c[2] = xyzq[i*6*3 + 5*3 + 2];
                        }
                        else{
                            a[0] = xyz[i*3 + 0];
                            a[1] = xyz[i*3 + 1];
                            a[2] = xyz[i*3 + 2];
                            b[0] = xyz[point_local[i*6 + 5]*3 + 0];
                            b[1] = xyz[point_local[i*6 + 5]*3 + 1];
                            b[2] = xyz[point_local[i*6 + 5]*3 + 2];
                            c[0] = xyzq[i*6*3 + 5*3 + 0];
                            c[1] = xyzq[i*6*3 + 5*3 + 1];
                            c[2] = xyzq[i*6*3 + 5*3 + 2];
                        }
                    }

                    for (int v = 0; v < 3; v++){
                        v01[v] = a[v];
                        v02[v] = b[v];
                        v03[v] = c[v];

                        v11[v] = b[v] - a[v];
                        v12[v] = a[v] - b[v];
                        v13[v] = a[v] - c[v];

                        v21[v] = c[v] - a[v];
                        v22[v] = c[v] - b[v];
                        v23[v] = b[v] - c[v];
                    }

                    t11 = 1/(v01[0]*v01[0] + v01[1]*v01[1] + v01[2]*v01[2]);

                    fac11       = (v01[0]*v11[0] + v01[1]*v11[1] + v01[2]*v11[2])*t11;
                    fac21       = (v01[0]*v21[0] + v01[1]*v21[1] + v01[2]*v21[2])*t11;

                    w11[0] = v11[0] - fac11*v01[0];
                    w11[1] = v11[1] - fac11*v01[1];
                    w11[2] = v11[2] - fac11*v01[2];

                    w21[0] = v21[0] - fac21*v01[0];
                    w21[1] = v21[1] - fac21*v01[1];
                    w21[2] = v21[2] - fac21*v01[2];

                    ang[0] = (w11[0]*w21[0] + w11[1]*w21[1] + w11[2]*w21[2])/
                        (sqrt(w11[0]*w11[0] + w11[1]*w11[1] + w11[2]*w11[2])*
                         sqrt(w21[0]*w21[0] + w21[1]*w21[1] + w21[2]*w21[2]));

                    if (ang[0] > 1)    ang[0] =  1;
                    if (ang[0] < -1)ang[0] = -1;

                    ang[0] = acos(ang[0]);

                    t22 = 1/(v02[0]*v02[0] + v02[1]*v02[1] + v02[2]*v02[2]);

                    fac12       = (v02[0]*v12[0] + v02[1]*v12[1] + v02[2]*v12[2])*t22;
                    fac22       = (v02[0]*v22[0] + v02[1]*v22[1] + v02[2]*v22[2])*t22;

                    w12[0] = v12[0] - fac12*v02[0];
                    w12[1] = v12[1] - fac12*v02[1];
                    w12[2] = v12[2] - fac12*v02[2];

                    w22[0] = v22[0] - fac22*v02[0];
                    w22[1] = v22[1] - fac22*v02[1];
                    w22[2] = v22[2] - fac22*v02[2];

                    ang[1] = (w12[0]*w22[0] + w12[1]*w22[1] + w12[2]*w22[2])/
                        (sqrt(w12[0]*w12[0] + w12[1]*w12[1] + w12[2]*w12[2])*
                         sqrt(w22[0]*w22[0] + w22[1]*w22[1] + w22[2]*w22[2]));

                    if ( ang[1] > 1) ang[1] =  1;
                    if ( ang[1] < -1)ang[1] = -1;

                    ang[1] = acos(ang[1]);

                    t33 = 1/(v03[0]*v03[0] + v03[1]*v03[1] + v03[2]*v03[2]);

                    fac13       = (v03[0]*v13[0] + v03[1]*v13[1] + v03[2]*v13[2])*t33;
                    fac23       = (v03[0]*v23[0] + v03[1]*v23[1] + v03[2]*v23[2])*t33;

                    w13[0] = v13[0] - fac13*v03[0];
                    w13[1] = v13[1] - fac13*v03[1];
                    w13[2] = v13[2] - fac13*v03[2];

                    w23[0] = v23[0] - fac23*v03[0];
                    w23[1] = v23[1] - fac23*v03[1];
                    w23[2] = v23[2] - fac23*v03[2];

                    ang[2] = (w13[0]*w23[0] + w13[1]*w23[1] + w13[2]*w23[2])/
                        (sqrt(w13[0]*w13[0] + w13[1]*w13[1] + w13[2]*w13[2])*
                         sqrt(w23[0]*w23[0] + w23[1]*w23[1] + w23[2]*w23[2]));

                    if (ang[2] > 1)    ang[2] =  1;
                    if (ang[2] < -1)ang[2] = -1;

                    ang[2] = acos(ang[2]);

                    radius = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]) ;
                    areasq[i*6*3 + j*3 + k] =
                        (ang[0] + ang[1] + ang[2] - M_PI)*pow(radius,2);
                }
                areasTr[i * 6 + j] = areasq[i * 6 * 3 + j * 3 + 0] + areasq[i * 6 * 3 + j * 3 + 1] + areasq[i * 6 * 3 + j * 3 + 2];
            }
        }
    }

    for (int i = 0; i < point_num; i++){
//
//      Control Volumes
//
        areasT[i] = 0.0;

        geo = 6; // Hexagons.
        for (int k = 0; k < 12; k++) if(i == pent_ind[k]) geo = 5; // Pentagons.
        for (int j = 0; j < 6; j++){
            if (geo == 5 && j == 4){
                a[0] = xyz[i*3 + 0];
                a[1] = xyz[i*3 + 1];
                a[2] = xyz[i*3 + 2];
                b[0] = xyzq[i*3*6 + 4*3 + 0];
                b[1] = xyzq[i*3*6 + 4*3 + 1];
                b[2] = xyzq[i*3*6 + 4*3 + 2];
                c[0] = xyzq[i*3*6 + 0*3 + 0];
                c[1] = xyzq[i*3*6 + 0*3 + 1];
                c[2] = xyzq[i*3*6 + 0*3 + 2];

            }
            else if (geo == 5 && j == 5){
                break;
            }
            else if (geo == 6 && j == 5){
                a[0] = xyz[i*3 + 0];
                a[1] = xyz[i*3 + 1];
                a[2] = xyz[i*3 + 2];
                b[0] = xyzq[i*3*6 + 5*3 + 0];
                b[1] = xyzq[i*3*6 + 5*3 + 1];
                b[2] = xyzq[i*3*6 + 5*3 + 2];
                c[0] = xyzq[i*3*6 + 0*3 + 0];
                c[1] = xyzq[i*3*6 + 0*3 + 1];
                c[2] = xyzq[i*3*6 + 0*3 + 2];
            }
            else{
                a[0] = xyz[i*3 + 0];
                a[1] = xyz[i*3 + 1];
                a[2] = xyz[i*3 + 2];
                b[0] = xyzq[i*3*6 + j*3 + 0];
                b[1] = xyzq[i*3*6 + j*3 + 1];
                b[2] = xyzq[i*3*6 + j*3 + 2];
                c[0] = xyzq[i*3*6 + (j+1)*3 + 0];
                c[1] = xyzq[i*3*6 + (j+1)*3 + 1];
                c[2] = xyzq[i*3*6 + (j+1)*3 + 2];
            }
            for (int v = 0; v < 3; v++){
                v01[v] = a[v];
                v02[v] = b[v];
                v03[v] = c[v];

                v11[v] = b[v] - a[v];
                v12[v] = a[v] - b[v];
                v13[v] = a[v] - c[v];

                v21[v] = c[v] - a[v];
                v22[v] = c[v] - b[v];
                v23[v] = b[v] - c[v];
            }

            t11 = 1/(v01[0]*v01[0] + v01[1]*v01[1] + v01[2]*v01[2]);

            fac11       = (v01[0]*v11[0] + v01[1]*v11[1] + v01[2]*v11[2])*t11;
            fac21       = (v01[0]*v21[0] + v01[1]*v21[1] + v01[2]*v21[2])*t11;

            w11[0] = v11[0] - fac11*v01[0];
            w11[1] = v11[1] - fac11*v01[1];
            w11[2] = v11[2] - fac11*v01[2];

            w21[0] = v21[0] - fac21*v01[0];
            w21[1] = v21[1] - fac21*v01[1];
            w21[2] = v21[2] - fac21*v01[2];

            ang[0] = (w11[0]*w21[0] + w11[1]*w21[1] + w11[2]*w21[2])/
                    (sqrt(w11[0]*w11[0] + w11[1]*w11[1] + w11[2]*w11[2])*
                     sqrt(w21[0]*w21[0] + w21[1]*w21[1] + w21[2]*w21[2]));

            if (ang[0] > 1)    ang[0] = 1 ;
            if (ang[0] < -1)ang[0] = -1;

            ang[0] = acos(ang[0]);

            t22 = 1/(v02[0]*v02[0] + v02[1]*v02[1] + v02[2]*v02[2]);
            fac12 = (v02[0]*v12[0] + v02[1]*v12[1] + v02[2]*v12[2])*t22;
            fac22 = (v02[0]*v22[0] + v02[1]*v22[1] + v02[2]*v22[2])*t22;

            w12[0] = v12[0] - fac12*v02[0];
            w12[1] = v12[1] - fac12*v02[1];
            w12[2] = v12[2] - fac12*v02[2];

            w22[0] = v22[0] - fac22*v02[0];
            w22[1] = v22[1] - fac22*v02[1];
            w22[2] = v22[2] - fac22*v02[2];

            ang[1] = (w12[0]*w22[0] + w12[1]*w22[1] + w12[2]*w22[2])/
                    (sqrt(w12[0]*w12[0] + w12[1]*w12[1] + w12[2]*w12[2])*
                     sqrt(w22[0]*w22[0] + w22[1]*w22[1] + w22[2]*w22[2]));

            if (ang[1] > 1)    ang[1] = 1 ;
            if (ang[1] <-1)    ang[1] = -1;

            ang[1] = acos(ang[1]);

            t33 = 1/(v03[0]*v03[0] + v03[1]*v03[1] + v03[2]*v03[2]);
            fac13 = (v03[0]*v13[0] + v03[1]*v13[1] + v03[2]*v13[2])*t33;
            fac23 = (v03[0]*v23[0] + v03[1]*v23[1] + v03[2]*v23[2])*t33;

            w13[0] = v13[0] - fac13*v03[0];
            w13[1] = v13[1] - fac13*v03[1];
            w13[2] = v13[2] - fac13*v03[2];

            w23[0] = v23[0] - fac23*v03[0];
            w23[1] = v23[1] - fac23*v03[1];
            w23[2] = v23[2] - fac23*v03[2];

            ang[2] = (w13[0]*w23[0] + w13[1]*w23[1] + w13[2]*w23[2])/
                    (sqrt(w13[0]*w13[0] + w13[1]*w13[1] + w13[2]*w13[2])*
                     sqrt(w23[0]*w23[0] + w23[1]*w23[1] + w23[2]*w23[2]));

            if (ang[2] > 1)    ang[2] = 1 ;
            if (ang[2] <-1)    ang[2] = -1;

            ang[2] = acos(ang[2]);

            radius = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]) ;
            areav  = (ang[0] + ang[1] + ang[2] - M_PI)*pow(radius,2);

            areasT[i] += areav;
        }
    }

    delete [] a;
    delete [] b;
    delete [] c;
    delete [] v01;
    delete [] v02;
    delete [] v03;
    delete [] v11;
    delete [] v12;
    delete [] v13;
    delete [] v21;
    delete [] v22;
    delete [] v23;
    delete [] w13;
    delete [] w23;
    delete [] ang;
}


void Icogrid::control_vec(double *nvec    ,
                          double *nevcoa  ,
                          double *nvecti  ,
                          double *nvecte  ,
                          double *areasT  ,
                          double *xyz     ,
                          double *xyzq    ,
                          int *point_local,
                          int *pent_ind   ,
                          int point_num   ){

//
//  Description:
//
//  Compute normal and tangential vectors.
//
//  Input: - xyz         - Vertices coordinates.
//         - xyzq        - Q-points coordinates.
//         - point_local - First neighbors indexes of each vertex.
//         - pent_ind    - Pentagon' indexes.
//         - point_num   - Number of vertices.
//
//  Output: - nvec   - Vectors normal to the edges of the control volume.
//          - nvecti - Vectors normal to the side edges of the triangles.
//          - nvecte - Vectors normal to the outward edges of the triangles.
//

    // Local variables.
    double vec_l, l;
    double fac_nv;

    int geo;

    // Local arrays.
    double *v1, *v2;
    double *nv;

    v1   = new double[3]();
    v2   = new double[3]();
    nv   = new double[3]();

    for (int i = 0; i < point_num; i++){

        geo = 6; // Hexagons.
        for (int k = 0; k < 12; k++) if(i == pent_ind[k]) geo = 5; // Pentagons.

        if(geo == 5){

            for (int j = 0; j < 5; j++){

                if (j < 4){
                    v2[0] = xyzq[i*6*3 + (j+1)*3 + 0];
                    v2[1] = xyzq[i*6*3 + (j+1)*3 + 1];
                    v2[2] = xyzq[i*6*3 + (j+1)*3 + 2];
                    v1[0] = xyzq[i*6*3 + j*3 + 0];
                    v1[1] = xyzq[i*6*3 + j*3 + 1];
                    v1[2] = xyzq[i*6*3 + j*3 + 2];
                }
                else{
                    v2[0] = xyzq[i*6*3 + 0*3 + 0];
                    v2[1] = xyzq[i*6*3 + 0*3 + 1];
                    v2[2] = xyzq[i*6*3 + 0*3 + 2];
                    v1[0] = xyzq[i*6*3 + 4*3 + 0];
                    v1[1] = xyzq[i*6*3 + 4*3 + 1];
                    v1[2] = xyzq[i*6*3 + 4*3 + 2];
                }

                vec_l = sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);

                l = acos((v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2])/(pow(vec_l,2.0)))*vec_l;

                nv[0] = v1[1]*v2[2] - v1[2]*v2[1];
                nv[1] = v1[2]*v2[0] - v1[0]*v2[2];
                nv[2] = v1[0]*v2[1] - v1[1]*v2[0];

                fac_nv = l/sqrt(nv[0]*nv[0] + nv[1]*nv[1] + nv[2]*nv[2]);

                nvec[i*6*3 + j*3 + 0] = nv[0]*fac_nv;
                nvec[i*6*3 + j*3 + 1] = nv[1]*fac_nv;
                nvec[i*6*3 + j*3 + 2] = nv[2]*fac_nv;

                nvecoa[i * 6 * 3 + j * 3 + 0] = nvec[i * 6 * 3 + j * 3 + 0] / areasT[i];
                nvecoa[i * 6 * 3 + j * 3 + 1] = nvec[i * 6 * 3 + j * 3 + 1] / areasT[i];
                nvecoa[i * 6 * 3 + j * 3 + 2] = nvec[i * 6 * 3 + j * 3 + 2] / areasT[i];
            }

            nvec[i * 6 * 3 + 5 * 3 + 0] = 0.0;
            nvec[i * 6 * 3 + 5 * 3 + 1] = 0.0;
            nvec[i * 6 * 3 + 5 * 3 + 2] = 0.0;

            nvecoa[i * 6 * 3 + 5 * 3 + 0] = 0.0;
            nvecoa[i * 6 * 3 + 5 * 3 + 1] = 0.0;
            nvecoa[i * 6 * 3 + 5 * 3 + 2] = 0.0;

        }
        else{

            for (int j = 0; j < 6; j++){

                if (j < 5){
                    v2[0] = xyzq[i*6*3 + (j+1)*3 + 0];
                    v2[1] = xyzq[i*6*3 + (j+1)*3 + 1];
                    v2[2] = xyzq[i*6*3 + (j+1)*3 + 2];
                    v1[0] = xyzq[i*6*3 + j*3 + 0];
                    v1[1] = xyzq[i*6*3 + j*3 + 1];
                    v1[2] = xyzq[i*6*3 + j*3 + 2];
                }
                else{
                    v2[0] = xyzq[i*6*3 + 0*3 + 0];
                    v2[1] = xyzq[i*6*3 + 0*3 + 1];
                    v2[2] = xyzq[i*6*3 + 0*3 + 2];
                    v1[0] = xyzq[i*6*3 + 5*3 + 0];
                    v1[1] = xyzq[i*6*3 + 5*3 + 1];
                    v1[2] = xyzq[i*6*3 + 5*3 + 2];
                }

                vec_l = sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);

                l = acos((v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2])/(pow(vec_l,2.0)))*vec_l;

                nv[0] = v1[1]*v2[2] - v1[2]*v2[1];
                nv[1] = v1[2]*v2[0] - v1[0]*v2[2];
                nv[2] = v1[0]*v2[1] - v1[1]*v2[0];

                fac_nv = l/sqrt(nv[0]*nv[0] + nv[1]*nv[1] + nv[2]*nv[2]);

                nvec[i*6*3 + j*3 + 0] = nv[0]*fac_nv;
                nvec[i*6*3 + j*3 + 1] = nv[1]*fac_nv;
                nvec[i*6*3 + j*3 + 2] = nv[2]*fac_nv;

                nvecoa[i * 6 * 3 + j * 3 + 0] = nvec[i * 6 * 3 + j * 3 + 0] / areasT[i];
                nvecoa[i * 6 * 3 + j * 3 + 1] = nvec[i * 6 * 3 + j * 3 + 1] / areasT[i];
                nvecoa[i * 6 * 3 + j * 3 + 2] = nvec[i * 6 * 3 + j * 3 + 2] / areasT[i];
            }
        }

// nvecti & nvecte! Triangles.
        if(geo == 5){

            for (int j = 0; j < 5; j++){

                v1[0] = xyz[point_local[i*6 + j]*3 + 0];
                v1[1] = xyz[point_local[i*6 + j]*3 + 1];
                v1[2] = xyz[point_local[i*6 + j]*3 + 2];
                v2[0] = xyz[i*3 + 0];
                v2[1] = xyz[i*3 + 1];
                v2[2] = xyz[i*3 + 2];

                vec_l = sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);

                l = acos((v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2])/(pow(vec_l,2.0)))*vec_l;

                nv[0] = v1[1]*v2[2] - v1[2]*v2[1];
                nv[1] = v1[2]*v2[0] - v1[0]*v2[2];
                nv[2] = v1[0]*v2[1] - v1[1]*v2[0];

                fac_nv = l/sqrt(nv[0]*nv[0] + nv[1]*nv[1] + nv[2]*nv[2]);

                nvecti[i*6*3 + j*3 + 0] = nv[0]*fac_nv;
                nvecti[i*6*3 + j*3 + 1] = nv[1]*fac_nv;
                nvecti[i*6*3 + j*3 + 2] = nv[2]*fac_nv;
            }
            nvecti[i * 6 * 3 + 5 * 3 + 0] = 0.0;
            nvecti[i * 6 * 3 + 5 * 3 + 1] = 0.0;
            nvecti[i * 6 * 3 + 5 * 3 + 2] = 0.0;
        }
        else{
            for (int j = 0; j < 6; j++){
                v1[0] = xyz[point_local[i*6 + j]*3 + 0];
                v1[1] = xyz[point_local[i*6 + j]*3 + 1];
                v1[2] = xyz[point_local[i*6 + j]*3 + 2];
                v2[0] = xyz[i*3 + 0];
                v2[1] = xyz[i*3 + 1];
                v2[2] = xyz[i*3 + 2];

                vec_l = sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);

                l = acos((v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2])/(pow(vec_l,2.0)))*vec_l;

                nv[0] = v1[1]*v2[2] - v1[2]*v2[1];
                nv[1] = v1[2]*v2[0] - v1[0]*v2[2];
                nv[2] = v1[0]*v2[1] - v1[1]*v2[0];

                fac_nv = l/sqrt(nv[0]*nv[0] + nv[1]*nv[1] + nv[2]*nv[2]);

                nvecti[i*6*3 + j*3 + 0] = nv[0]*fac_nv;
                nvecti[i*6*3 + j*3 + 1] = nv[1]*fac_nv;
                nvecti[i*6*3 + j*3 + 2] = nv[2]*fac_nv;
            }
        }

        if(geo == 5){

            for (int j = 0; j < 5; j++){

                if (j < 4){
                    v1[0] = xyz[point_local[i*6 + j]*3 + 0];
                    v1[1] = xyz[point_local[i*6 + j]*3 + 1];
                    v1[2] = xyz[point_local[i*6 + j]*3 + 2];
                    v2[0] = xyz[point_local[i*6 + j+1]*3 + 0];
                    v2[1] = xyz[point_local[i*6 + j+1]*3 + 1];
                    v2[2] = xyz[point_local[i*6 + j+1]*3 + 2];
                }
                else{
                    v1[0] = xyz[point_local[i*6 + 4]*3 + 0];
                    v1[1] = xyz[point_local[i*6 + 4]*3 + 1];
                    v1[2] = xyz[point_local[i*6 + 4]*3 + 2];
                    v2[0] = xyz[point_local[i*6 + 0]*3 + 0];
                    v2[1] = xyz[point_local[i*6 + 0]*3 + 1];
                    v2[2] = xyz[point_local[i*6 + 0]*3 + 2];
                }

                vec_l = sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);

                l = acos((v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2])/(pow(vec_l,2.0)))*vec_l;

                nv[0] = v1[1]*v2[2] - v1[2]*v2[1];
                nv[1] = v1[2]*v2[0] - v1[0]*v2[2];
                nv[2] = v1[0]*v2[1] - v1[1]*v2[0];

                fac_nv = l/sqrt(nv[0]*nv[0] + nv[1]*nv[1] + nv[2]*nv[2]);

                nvecte[i*6*3 + j*3 + 0] = nv[0]*fac_nv;
                nvecte[i*6*3 + j*3 + 1] = nv[1]*fac_nv;
                nvecte[i*6*3 + j*3 + 2] = nv[2]*fac_nv;
            }
            nvecte[i * 6 * 3 + 5 * 3 + 0] = 0.0;
            nvecte[i * 6 * 3 + 5 * 3 + 1] = 0.0;
            nvecte[i * 6 * 3 + 5 * 3 + 2] = 0.0;
        }
        else{
            for (int j = 0; j < 6; j++){

                if (j < 5){
                    v1[0] = xyz[point_local[i*6 + j]*3 + 0];
                    v1[1] = xyz[point_local[i*6 + j]*3 + 1];
                    v1[2] = xyz[point_local[i*6 + j]*3 + 2];
                    v2[0] = xyz[point_local[i*6 + j+1]*3 + 0];
                    v2[1] = xyz[point_local[i*6 + j+1]*3 + 1];
                    v2[2] = xyz[point_local[i*6 + j+1]*3 + 2];
                }
                else{
                    v1[0] = xyz[point_local[i*6 + 5]*3 + 0];
                    v1[1] = xyz[point_local[i*6 + 5]*3 + 1];
                    v1[2] = xyz[point_local[i*6 + 5]*3 + 2];
                    v2[0] = xyz[point_local[i*6 + 0]*3 + 0];
                    v2[1] = xyz[point_local[i*6 + 0]*3 + 1];
                    v2[2] = xyz[point_local[i*6 + 0]*3 + 2];
                }

                vec_l = sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);

                l = acos((v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2])/(pow(vec_l,2.0)))*vec_l;

                nv[0] = v1[1]*v2[2] - v1[2]*v2[1];
                nv[1] = v1[2]*v2[0] - v1[0]*v2[2];
                nv[2] = v1[0]*v2[1] - v1[1]*v2[0];

                fac_nv = l/sqrt(nv[0]*nv[0] + nv[1]*nv[1] + nv[2]*nv[2]);

                nvecte[i*6*3 + j*3 + 0] = nv[0]*fac_nv;
                nvecte[i*6*3 + j*3 + 1] = nv[1]*fac_nv;
                nvecte[i*6*3 + j*3 + 2] = nv[2]*fac_nv;
            }
        }
    }

    delete [] v1;
    delete [] v2;
    delete [] nv;
}


void Icogrid::compute_func(double *func_r,
                           double *xyz   ,
                           int point_num ){

//
//  Description:
//
//  Creates radial vectors with norm 1 at the center of each control volume.
//
//  Input:  - xyz       - Position of the centroids.
//          - point_num - Number of vertices.
//
//  Output: - func_r    - vector position with norm 1.
//

    // Local variables:
    double norm;

    for (int i = 0; i < point_num; i++){

        norm = sqrt(pow(xyz[i*3 + 0],2.0) + pow(xyz[i*3 + 1],2.0) + pow(xyz[i*3 + 2],2.0));

        func_r[i*3 + 0] = xyz[i*3 + 0]/norm;
        func_r[i*3 + 1] = xyz[i*3 + 1]/norm;
        func_r[i*3 + 2] = xyz[i*3 + 2]/norm;

    }
}

void Icogrid::div_operator(double *areasT,
                           double *areas ,
                           double *div   ,
                           double *nvec  ,
                           int *pent_ind ,
                           int point_num){

//
//  Description:
//
//  Computes the divergence operator from equations 15 and 16 of M. Satoh et al. 2008.
//
//  Input: - areas     - Sub-areas (alpha, beta and gamma).
//         - areasT    - Control volume area.
//         - nvec      - Vectors normal to the edges of the control volume.
//         - pent_ind  - Pentagon' indexes.
//         - point_num - Number of vertices.
//
//  Output: - div - divergence operator.
//
    // Local variables:
    int geo;
    double area1,
           area2,
           area3,
           area4,
           area5,
           area6;

    double *areasq2;

    areasq2 = new double[6 * 3 * point_num]();

    for (int i = 0; i < point_num; i++){

        geo = 6; // Hexagons.
        for (int k = 0; k < 12; k++) if (i == pent_ind[k]) geo = 5; // Pentagons.

        area1 = areas[i * 6 * 3 + 0 * 3 + 0] +
                areas[i * 6 * 3 + 0 * 3 + 1] +
                areas[i * 6 * 3 + 0 * 3 + 2];
        area2 = areas[i * 6 * 3 + 1 * 3 + 0] +
                areas[i * 6 * 3 + 1 * 3 + 1] +
                areas[i * 6 * 3 + 1 * 3 + 2];
        area3 = areas[i * 6 * 3 + 2 * 3 + 0] +
                areas[i * 6 * 3 + 2 * 3 + 1] +
                areas[i * 6 * 3 + 2 * 3 + 2];
        area4 = areas[i * 6 * 3 + 3 * 3 + 0] +
                areas[i * 6 * 3 + 3 * 3 + 1] +
                areas[i * 6 * 3 + 3 * 3 + 2];
        area5 = areas[i * 6 * 3 + 4 * 3 + 0] +
                areas[i * 6 * 3 + 4 * 3 + 1] +
                areas[i * 6 * 3 + 4 * 3 + 2];

        if (geo == 5) area6 = 0;
        else area6 = areas[i * 6 * 3 + 5 * 3 + 0] +
                     areas[i * 6 * 3 + 5 * 3 + 1] +
                     areas[i * 6 * 3 + 5 * 3 + 2];

        for (int k = 0; k < 3; k++){
            areasq2[i * 6 * 3 + 0 * 3 + k] = areas[i * 6 * 3 + 0 * 3 + k] / area1;
            areasq2[i * 6 * 3 + 1 * 3 + k] = areas[i * 6 * 3 + 1 * 3 + k] / area2;
            areasq2[i * 6 * 3 + 2 * 3 + k] = areas[i * 6 * 3 + 2 * 3 + k] / area3;
            areasq2[i * 6 * 3 + 3 * 3 + k] = areas[i * 6 * 3 + 3 * 3 + k] / area4;
            areasq2[i * 6 * 3 + 4 * 3 + k] = areas[i * 6 * 3 + 4 * 3 + k] / area5;
            if (geo == 5)areasq2[i * 6 * 3 + 5 * 3 + k] = 0;
            else areasq2[i * 6 * 3 + 5 * 3 + k] =    areas[i * 6 * 3 + 5 * 3 + k] / area6;
        }
    }

    for (int i = 0; i < point_num; i++){

        geo = 6; // Hexagons.
        for (int k = 0; k < 12; k++) if (i == pent_ind[k]) geo = 5; // Pentagons.

        for (int k = 0; k < 3; k++){
            if (geo == 5){    //Pentagons.
                div[i * 7 * 3 + 0 * 3 + k] = nvec[i * 6 * 3 + 0 * 3 + k] * (areasq2[i * 6 * 3 + 0 * 3 + 0] + areasq2[i * 6 * 3 + 1 * 3 + 0]) +
                                             nvec[i * 6 * 3 + 1 * 3 + k] * (areasq2[i * 6 * 3 + 1 * 3 + 0] + areasq2[i * 6 * 3 + 2 * 3 + 0]) +
                                             nvec[i * 6 * 3 + 2 * 3 + k] * (areasq2[i * 6 * 3 + 2 * 3 + 0] + areasq2[i * 6 * 3 + 3 * 3 + 0]) +
                                             nvec[i * 6 * 3 + 3 * 3 + k] * (areasq2[i * 6 * 3 + 3 * 3 + 0] + areasq2[i * 6 * 3 + 4 * 3 + 0]) +
                                             nvec[i * 6 * 3 + 4 * 3 + k] * (areasq2[i * 6 * 3 + 4 * 3 + 0] + areasq2[i * 6 * 3 + 0 * 3 + 0]);

                div[i * 7 * 3 + 1 * 3 + k] = nvec[i * 6 * 3 + 4 * 3 + k] * (areasq2[i * 6 * 3 + 0 * 3 + 1] + areasq2[i * 6 * 3 + 4 * 3 + 2]) +
                                             nvec[i * 6 * 3 + 3 * 3 + k] * areasq2[i * 6 * 3 + 4 * 3 + 2] +
                                             nvec[i * 6 * 3 + 0 * 3 + k] * areasq2[i * 6 * 3 + 0 * 3 + 1];

                div[i * 7 * 3 + 2 * 3 + k] = nvec[i * 6 * 3 + 0 * 3 + k] * (areasq2[i * 6 * 3 + 1 * 3 + 1] + areasq2[i * 6 * 3 + 0 * 3 + 2]) +
                                             nvec[i * 6 * 3 + 4 * 3 + k] * areasq2[i * 6 * 3 + 0 * 3 + 2] +
                                             nvec[i * 6 * 3 + 1 * 3 + k] * areasq2[i * 6 * 3 + 1 * 3 + 1];

                div[i * 7 * 3 + 3 * 3 + k] = nvec[i * 6 * 3 + 1 * 3 + k] * (areasq2[i * 6 * 3 + 2 * 3 + 1] + areasq2[i * 6 * 3 + 1 * 3 + 2]) +
                                             nvec[i * 6 * 3 + 0 * 3 + k] * areasq2[i * 6 * 3 + 1 * 3 + 2] +
                                             nvec[i * 6 * 3 + 2 * 3 + k] * areasq2[i * 6 * 3 + 2 * 3 + 1];

                div[i * 7 * 3 + 4 * 3 + k] = nvec[i * 6 * 3 + 2 * 3 + k] * (areasq2[i * 6 * 3 + 3 * 3 + 1] + areasq2[i * 6 * 3 + 2 * 3 + 2]) +
                                             nvec[i * 6 * 3 + 1 * 3 + k] * areasq2[i * 6 * 3 + 2 * 3 + 2] +
                                             nvec[i * 6 * 3 + 3 * 3 + k] * areasq2[i * 6 * 3 + 3 * 3 + 1];

                div[i * 7 * 3 + 5 * 3 + k] = nvec[i * 6 * 3 + 3 * 3 + k] * (areasq2[i * 6 * 3 + 4 * 3 + 1] + areasq2[i * 6 * 3 + 3 * 3 + 2]) +
                                             nvec[i * 6 * 3 + 2 * 3 + k] * areasq2[i * 6 * 3 + 3 * 3 + 2] +
                                             nvec[i * 6 * 3 + 4 * 3 + k] * areasq2[i * 6 * 3 + 4 * 3 + 1];

                div[i * 7 * 3 + 6 * 3 + k] = 0.0;
            }
            else{
                div[i * 7 * 3 + 0 * 3 + k] = nvec[i * 6 * 3 + 0 * 3 + k] * (areasq2[i * 6 * 3 + 0 * 3 + 0] + areasq2[i * 6 * 3 + 1 * 3 + 0]) +
                                             nvec[i * 6 * 3 + 1 * 3 + k] * (areasq2[i * 6 * 3 + 1 * 3 + 0] + areasq2[i * 6 * 3 + 2 * 3 + 0]) +
                                             nvec[i * 6 * 3 + 2 * 3 + k] * (areasq2[i * 6 * 3 + 2 * 3 + 0] + areasq2[i * 6 * 3 + 3 * 3 + 0]) +
                                             nvec[i * 6 * 3 + 3 * 3 + k] * (areasq2[i * 6 * 3 + 3 * 3 + 0] + areasq2[i * 6 * 3 + 4 * 3 + 0]) +
                                             nvec[i * 6 * 3 + 4 * 3 + k] * (areasq2[i * 6 * 3 + 4 * 3 + 0] + areasq2[i * 6 * 3 + 5 * 3 + 0]) +
                                             nvec[i * 6 * 3 + 5 * 3 + k] * (areasq2[i * 6 * 3 + 5 * 3 + 0] + areasq2[i * 6 * 3 + 0 * 3 + 0]);

                div[i * 7 * 3 + 1 * 3 + k] = nvec[i * 6 * 3 + 5 * 3 + k] * (areasq2[i * 6 * 3 + 0 * 3 + 1] + areasq2[i * 6 * 3 + 5 * 3 + 2]) +
                                             nvec[i * 6 * 3 + 4 * 3 + k] * areasq2[i * 6 * 3 + 5 * 3 + 2] +
                                             nvec[i * 6 * 3 + 0 * 3 + k] * areasq2[i * 6 * 3 + 0 * 3 + 1];

                div[i * 7 * 3 + 2 * 3 + k] = nvec[i * 6 * 3 + 0 * 3 + k] * (areasq2[i * 6 * 3 + 1 * 3 + 1] + areasq2[i * 6 * 3 + 0 * 3 + 2]) +
                                             nvec[i * 6 * 3 + 5 * 3 + k] * areasq2[i * 6 * 3 + 0 * 3 + 2] +
                                             nvec[i * 6 * 3 + 1 * 3 + k] * areasq2[i * 6 * 3 + 1 * 3 + 1];

                div[i * 7 * 3 + 3 * 3 + k] = nvec[i * 6 * 3 + 1 * 3 + k] * (areasq2[i * 6 * 3 + 2 * 3 + 1] + areasq2[i * 6 * 3 + 1 * 3 + 2]) +
                                             nvec[i * 6 * 3 + 0 * 3 + k] * areasq2[i * 6 * 3 + 1 * 3 + 2] +
                                             nvec[i * 6 * 3 + 2 * 3 + k] * areasq2[i * 6 * 3 + 2 * 3 + 1];

                div[i * 7 * 3 + 4 * 3 + k] = nvec[i * 6 * 3 + 2 * 3 + k] * (areasq2[i * 6 * 3 + 3 * 3 + 1] + areasq2[i * 6 * 3 + 2 * 3 + 2]) +
                                             nvec[i * 6 * 3 + 1 * 3 + k] * areasq2[i * 6 * 3 + 2 * 3 + 2] +
                                             nvec[i * 6 * 3 + 3 * 3 + k] * areasq2[i * 6 * 3 + 3 * 3 + 1];

                div[i * 7 * 3 + 5 * 3 + k] = nvec[i * 6 * 3 + 3 * 3 + k] * (areasq2[i * 6 * 3 + 4 * 3 + 1] + areasq2[i * 6 * 3 + 3 * 3 + 2]) +
                                             nvec[i * 6 * 3 + 2 * 3 + k] * areasq2[i * 6 * 3 + 3 * 3 + 2] +
                                             nvec[i * 6 * 3 + 4 * 3 + k] * areasq2[i * 6 * 3 + 4 * 3 + 1];

                div[i * 7 * 3 + 6 * 3 + k] = nvec[i * 6 * 3 + 4 * 3 + k] * (areasq2[i * 6 * 3 + 5 * 3 + 1] + areasq2[i * 6 * 3 + 4 * 3 + 2]) +
                                             nvec[i * 6 * 3 + 3 * 3 + k] * areasq2[i * 6 * 3 + 4 * 3 + 2] +
                                             nvec[i * 6 * 3 + 5 * 3 + k] * areasq2[i * 6 * 3 + 5 * 3 + 1];
            }
        }
        for (int j = 0; j < 7; j++)for (int k = 0; k < 3; k++)    div[i * 7 * 3 + j * 3 + k] = div[i * 7 * 3 + j * 3 + k] / (2.0*areasT[i]);
    }
    delete[] areasq2;
}

void Icogrid::gra_operator(double *areasT,
                           double *areas ,
                           double *grad  ,
                           double *nvec  ,
                           int *pent_ind ,
                           int point_num ){
//
//  Description:
//
//  Computes the gradient operator from equation 17 of M. Satoh et al. 2008.
//
//  Input: - areas     - Sub-areas (alpha, beta and gamma).
//         - areasT    - Control volume area.
//         - nvec      - Vectors normal to the edges of the control volume.
//         - pent_ind  - Pentagon' indexes.
//         - point_num - Number of vertices.
//
//  Output: - grad - gradient operator.
//
    // Local variables:
    int geo;
    double area1,
           area2,
           area3,
           area4,
           area5,
           area6;

    double *areasq2;

    areasq2 = new double[6 * 3 * point_num]();

    for (int i = 0; i < point_num; i++){

        geo = 6; // Hexagons.
        for (int k = 0; k < 12; k++) if (i == pent_ind[k]) geo = 5; // Pentagons.

        area1 = areas[i * 6 * 3 + 0 * 3 + 0] +
                areas[i * 6 * 3 + 0 * 3 + 1] +
                areas[i * 6 * 3 + 0 * 3 + 2];
        area2 = areas[i * 6 * 3 + 1 * 3 + 0] +
                areas[i * 6 * 3 + 1 * 3 + 1] +
                areas[i * 6 * 3 + 1 * 3 + 2];
        area3 = areas[i * 6 * 3 + 2 * 3 + 0] +
                areas[i * 6 * 3 + 2 * 3 + 1] +
                areas[i * 6 * 3 + 2 * 3 + 2];
        area4 = areas[i * 6 * 3 + 3 * 3 + 0] +
                areas[i * 6 * 3 + 3 * 3 + 1] +
                areas[i * 6 * 3 + 3 * 3 + 2];
        area5 = areas[i * 6 * 3 + 4 * 3 + 0] +
                areas[i * 6 * 3 + 4 * 3 + 1] +
                areas[i * 6 * 3 + 4 * 3 + 2];

        if (geo == 5) area6 = 0;
        else area6 = areas[i * 6 * 3 + 5 * 3 + 0] +
                     areas[i * 6 * 3 + 5 * 3 + 1] +
                     areas[i * 6 * 3 + 5 * 3 + 2];

        for (int k = 0; k < 3; k++){
            areasq2[i * 6 * 3 + 0 * 3 + k] = areas[i * 6 * 3 + 0 * 3 + k] / area1;
            areasq2[i * 6 * 3 + 1 * 3 + k] = areas[i * 6 * 3 + 1 * 3 + k] / area2;
            areasq2[i * 6 * 3 + 2 * 3 + k] = areas[i * 6 * 3 + 2 * 3 + k] / area3;
            areasq2[i * 6 * 3 + 3 * 3 + k] = areas[i * 6 * 3 + 3 * 3 + k] / area4;
            areasq2[i * 6 * 3 + 4 * 3 + k] = areas[i * 6 * 3 + 4 * 3 + k] / area5;
            if (geo == 5) areasq2[i * 6 * 3 + 5 * 3 + k] = 0;
            else areasq2[i * 6 * 3 + 5 * 3 + k] = areas[i * 6 * 3 + 5 * 3 + k] / area6;
        }
    }

    for (int i = 0; i < point_num; i++){

        geo = 6; // Hexagons.

        for (int k = 0; k < 12; k++) if (i == pent_ind[k]) geo = 5; // Pentagons.
        for (int k = 0; k < 3; k++){
            if (geo == 5){    //Pentagons.
                grad[i * 7 * 3 + 0 * 3 + k] = nvec[i * 6 * 3 + 0 * 3 + k] * (areasq2[i * 6 * 3 + 0 * 3 + 0] + areasq2[i * 6 * 3 + 1 * 3 + 0]) +
                                              nvec[i * 6 * 3 + 1 * 3 + k] * (areasq2[i * 6 * 3 + 1 * 3 + 0] + areasq2[i * 6 * 3 + 2 * 3 + 0]) +
                                              nvec[i * 6 * 3 + 2 * 3 + k] * (areasq2[i * 6 * 3 + 2 * 3 + 0] + areasq2[i * 6 * 3 + 3 * 3 + 0]) +
                                              nvec[i * 6 * 3 + 3 * 3 + k] * (areasq2[i * 6 * 3 + 3 * 3 + 0] + areasq2[i * 6 * 3 + 4 * 3 + 0]) +
                                              nvec[i * 6 * 3 + 4 * 3 + k] * (areasq2[i * 6 * 3 + 4 * 3 + 0] + areasq2[i * 6 * 3 + 0 * 3 + 0]);
                grad[i * 7 * 3 + 0 * 3 + k] = grad[i * 7 * 3 + 0 * 3 + k] -
                                              2 * (nvec[i * 6 * 3 + 0 * 3 + k] + nvec[i * 6 * 3 + 1 * 3 + k] + nvec[i * 6 * 3 + 2 * 3 + k] +
                                              nvec[i * 6 * 3 + 3 * 3 + k] + nvec[i * 6 * 3 + 4 * 3 + k]);

                grad[i * 7 * 3 + 1 * 3 + k] = nvec[i * 6 * 3 + 4 * 3 + k] * (areasq2[i * 6 * 3 + 0 * 3 + 1] + areasq2[i * 6 * 3 + 4 * 3 + 2]) +
                                              nvec[i * 6 * 3 + 3 * 3 + k] * areasq2[i * 6 * 3 + 4 * 3 + 2] +
                                              nvec[i * 6 * 3 + 0 * 3 + k] * areasq2[i * 6 * 3 + 0 * 3 + 1];

                grad[i * 7 * 3 + 2 * 3 + k] = nvec[i * 6 * 3 + 0 * 3 + k] * (areasq2[i * 6 * 3 + 1 * 3 + 1] + areasq2[i * 6 * 3 + 0 * 3 + 2]) +
                                              nvec[i * 6 * 3 + 4 * 3 + k] * areasq2[i * 6 * 3 + 0 * 3 + 2] +
                                              nvec[i * 6 * 3 + 1 * 3 + k] * areasq2[i * 6 * 3 + 1 * 3 + 1];

                grad[i * 7 * 3 + 3 * 3 + k] = nvec[i * 6 * 3 + 1 * 3 + k] * (areasq2[i * 6 * 3 + 2 * 3 + 1] + areasq2[i * 6 * 3 + 1 * 3 + 2]) +
                                              nvec[i * 6 * 3 + 0 * 3 + k] * areasq2[i * 6 * 3 + 1 * 3 + 2] +
                                              nvec[i * 6 * 3 + 2 * 3 + k] * areasq2[i * 6 * 3 + 2 * 3 + 1];

                grad[i * 7 * 3 + 4 * 3 + k] = nvec[i * 6 * 3 + 2 * 3 + k] * (areasq2[i * 6 * 3 + 3 * 3 + 1] + areasq2[i * 6 * 3 + 2 * 3 + 2]) +
                                              nvec[i * 6 * 3 + 1 * 3 + k] * areasq2[i * 6 * 3 + 2 * 3 + 2] +
                                              nvec[i * 6 * 3 + 3 * 3 + k] * areasq2[i * 6 * 3 + 3 * 3 + 1];

                grad[i * 7 * 3 + 5 * 3 + k] = nvec[i * 6 * 3 + 3 * 3 + k] * (areasq2[i * 6 * 3 + 4 * 3 + 1] + areasq2[i * 6 * 3 + 3 * 3 + 2]) +
                                              nvec[i * 6 * 3 + 2 * 3 + k] * areasq2[i * 6 * 3 + 3 * 3 + 2] +
                                              nvec[i * 6 * 3 + 4 * 3 + k] * areasq2[i * 6 * 3 + 4 * 3 + 1];

                grad[i * 7 * 3 + 6 * 3 + k] = 0.0;
            }
            else{
                grad[i * 7 * 3 + 0 * 3 + k] = nvec[i * 6 * 3 + 0 * 3 + k] * (areasq2[i * 6 * 3 + 0 * 3 + 0] + areasq2[i * 6 * 3 + 1 * 3 + 0]) +
                                              nvec[i * 6 * 3 + 1 * 3 + k] * (areasq2[i * 6 * 3 + 1 * 3 + 0] + areasq2[i * 6 * 3 + 2 * 3 + 0]) +
                                              nvec[i * 6 * 3 + 2 * 3 + k] * (areasq2[i * 6 * 3 + 2 * 3 + 0] + areasq2[i * 6 * 3 + 3 * 3 + 0]) +
                                              nvec[i * 6 * 3 + 3 * 3 + k] * (areasq2[i * 6 * 3 + 3 * 3 + 0] + areasq2[i * 6 * 3 + 4 * 3 + 0]) +
                                              nvec[i * 6 * 3 + 4 * 3 + k] * (areasq2[i * 6 * 3 + 4 * 3 + 0] + areasq2[i * 6 * 3 + 5 * 3 + 0]) +
                                              nvec[i * 6 * 3 + 5 * 3 + k] * (areasq2[i * 6 * 3 + 5 * 3 + 0] + areasq2[i * 6 * 3 + 0 * 3 + 0]);
                grad[i * 7 * 3 + 0 * 3 + k] = grad[i * 7 * 3 + 0 * 3 + k] -
                                              2 * (nvec[i * 6 * 3 + 0 * 3 + k] + nvec[i * 6 * 3 + 1 * 3 + k] + nvec[i * 6 * 3 + 2 * 3 + k] +
                                              nvec[i * 6 * 3 + 3 * 3 + k] + nvec[i * 6 * 3 + 4 * 3 + k] + nvec[i * 6 * 3 + 5 * 3 + k]);

                grad[i * 7 * 3 + 1 * 3 + k] = nvec[i * 6 * 3 + 5 * 3 + k] * (areasq2[i * 6 * 3 + 0 * 3 + 1] + areasq2[i * 6 * 3 + 5 * 3 + 2]) +
                                              nvec[i * 6 * 3 + 4 * 3 + k] * areasq2[i * 6 * 3 + 5 * 3 + 2] +
                                              nvec[i * 6 * 3 + 0 * 3 + k] * areasq2[i * 6 * 3 + 0 * 3 + 1];

                grad[i * 7 * 3 + 2 * 3 + k] = nvec[i * 6 * 3 + 0 * 3 + k] * (areasq2[i * 6 * 3 + 1 * 3 + 1] + areasq2[i * 6 * 3 + 0 * 3 + 2]) +
                                              nvec[i * 6 * 3 + 5 * 3 + k] * areasq2[i * 6 * 3 + 0 * 3 + 2] +
                                              nvec[i * 6 * 3 + 1 * 3 + k] * areasq2[i * 6 * 3 + 1 * 3 + 1];

                grad[i * 7 * 3 + 3 * 3 + k] = nvec[i * 6 * 3 + 1 * 3 + k] * (areasq2[i * 6 * 3 + 2 * 3 + 1] + areasq2[i * 6 * 3 + 1 * 3 + 2]) +
                                              nvec[i * 6 * 3 + 0 * 3 + k] * areasq2[i * 6 * 3 + 1 * 3 + 2] +
                                              nvec[i * 6 * 3 + 2 * 3 + k] * areasq2[i * 6 * 3 + 2 * 3 + 1];

                grad[i * 7 * 3 + 4 * 3 + k] = nvec[i * 6 * 3 + 2 * 3 + k] * (areasq2[i * 6 * 3 + 3 * 3 + 1] + areasq2[i * 6 * 3 + 2 * 3 + 2]) +
                                              nvec[i * 6 * 3 + 1 * 3 + k] * areasq2[i * 6 * 3 + 2 * 3 + 2] +
                                              nvec[i * 6 * 3 + 3 * 3 + k] * areasq2[i * 6 * 3 + 3 * 3 + 1];

                grad[i * 7 * 3 + 5 * 3 + k] = nvec[i * 6 * 3 + 3 * 3 + k] * (areasq2[i * 6 * 3 + 4 * 3 + 1] + areasq2[i * 6 * 3 + 3 * 3 + 2]) +
                                              nvec[i * 6 * 3 + 2 * 3 + k] * areasq2[i * 6 * 3 + 3 * 3 + 2] +
                                              nvec[i * 6 * 3 + 4 * 3 + k] * areasq2[i * 6 * 3 + 4 * 3 + 1];

                grad[i * 7 * 3 + 6 * 3 + k] = nvec[i * 6 * 3 + 4 * 3 + k] * (areasq2[i * 6 * 3 + 5 * 3 + 1] + areasq2[i * 6 * 3 + 4 * 3 + 2]) +
                                              nvec[i * 6 * 3 + 3 * 3 + k] * areasq2[i * 6 * 3 + 4 * 3 + 2] +
                                              nvec[i * 6 * 3 + 5 * 3 + k] * areasq2[i * 6 * 3 + 5 * 3 + 1];
            }
        }
        for (int j = 0; j < 7; j++)    for (int k = 0; k < 3; k++)    grad[i * 7 * 3 + j * 3 + k] = grad[i * 7 * 3 + j * 3 + k] / (2.0*areasT[i]);
    }
    delete[] areasq2;
}

void Icogrid::zonal_mean_tab_f(int    *zonal_mean_tab,
							   double *lonlat        ,
							   int nlat              ,
						       int point_num         ){

//
//  Description:
//
//  .
//
//  Input: - .
//         - point_num - Number of vertices.
//
//  Output: - .
//
	// Local variables:

	double des_lat = M_PI/nlat;
	int *count_num;
	double *lat_array;

	count_num   = new int[nlat]();
	lat_array   = new double[nlat+1]();

	for (int i = 0; i < point_num; i++)	for (int k = 0; k < 2; k++) zonal_mean_tab[i*2 + k]	= 0;

	for (int j = 0; j < nlat; j++)	count_num[j] = 0;
	lat_array[0] = -M_PI/2.0;
	for (int j = 1; j < nlat+1; j++) lat_array[j] = lat_array[j-1] + des_lat;
	for (int i = 0; i < point_num; i++){
		for (int j = 0; j < nlat; j++){
			if(lonlat[i*2 + 1] >= lat_array[j] && lonlat[i*2 + 1] < lat_array[j+1]){
				zonal_mean_tab[i*2] = j;
				count_num[j] = count_num[j] + 1;
				break;
			}
		}
		if(lonlat[i*2 + 1] >= lat_array[nlat]){
			zonal_mean_tab[i*2] = nlat-1;
			count_num[nlat-1] = count_num[nlat-1] + 1;
		}
	}
	for (int i = 0; i < point_num; i++){
		int ind = zonal_mean_tab[i*2];
		zonal_mean_tab[i*2 + 1] = count_num[ind];
	}

	delete[] count_num;
	delete[] lat_array;
}


//END OF GRID.CU
