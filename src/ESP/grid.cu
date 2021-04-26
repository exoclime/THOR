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

#include "grid.h"
#include "log_writer.h"
#include "vector_operations.h"

// some helper local functions
inline double3 normproj(const double3 &v1, const double3 &v2) {


    double l1 = dot(v1, v2);

    double3 vc1 = cross(v1, v2);

    double l2 = length(vc1);

    vc1 = vc1 / l2 * atan2(l2, l1);

    return vc1;
}


__host__ Icogrid::Icogrid(bool   sprd,        // Spring dynamics option
                          double spring_beta, // Parameter beta for spring dynamics
                          int    glevel,      // Horizontal resolution level
                          int    vlevel,      // Number of vertical layers
                          int    nlat,
                          double A,            // Planet radius [m]
                          double Top_altitude, // Top model's domain [m]
                          bool   sponge,
                          int *  max_count,
                          bool   vert_refined,
                          double lowest_layer_thickness,
                          double transition_altitude) {

    log::printf("\n\n Building icosahedral grid!");

    //  Number of vertices
    point_num = 2 + 10 * pow(2.0, 2.0 * glevel); // Horizontal
    nv        = vlevel;                          // Vertical
    nvi       = vlevel + 1;                      // Interfaces between layers

    //  Number of times the main rhombi are divided.
    int divide_face; // Used to split memory on the GPU
    if (glevel == 4)
        divide_face = 0;
    else if (glevel == 5)
        divide_face = 1;
    else if (glevel == 6)
        divide_face = 2;
    else if (glevel == 7)
        divide_face = 3;
    else {
        log::printf("\nParameter not tested! Use the predefined divide_face values.\n");
        exit(EXIT_FAILURE);
    }
    grid_level = glevel;

    //  Rhombi
    int n_region = pow(2.0, glevel) * pow(2.0, glevel); //
    nl_region    = pow(2.0, glevel);                    //
    nl_region    = nl_region / (pow(2.0, divide_face)); //
    int nl2      = nl_region * nl_region;               //
    int nfaces   = pow(4.0, divide_face);               //
    int nlhalo   = nl_region + 2;                       //
    int kxl      = sqrt(nfaces * 1.0);                  //
    nr           = (point_num - 2) / nl2;               //

    //  Compute standard grid.
    point_xyz = (double *)malloc(3 * point_num * sizeof(double));
    pent_ind  = (int *)malloc(12 * sizeof(int));
    sphere_ico(
        point_xyz, glevel, n_region, nl_region, nl2, kxl, nfaces, pent_ind, divide_face, point_num);

    //  Finds the closest neighbors of each point.
    point_local = (int *)malloc(6 * point_num * sizeof(int));
    neighbors_indx(point_local, point_xyz, pent_ind, point_num);

    //  Reorder neighbors in the clockwise order.
    reorder_neighbors_indx(point_local, point_xyz, pent_ind, point_num);

    //  Generate halos.
    nh   = 10 * nfaces * 4 * nlhalo;
    halo = (int *)malloc(10 * nfaces * 4 * nlhalo * sizeof(int));
    generate_halos(halo, point_local, n_region, divide_face);

    //  Reorder neighbors consistent with the rhombi.
    reorder_neighbors_indx_rhombi(
        point_local, halo, pent_ind, nl_region, nl2, nr, nlhalo, point_num);

    //  Finds the closest neighbors at the pole.
    neighbors_indx_pl(point_local, point_xyz, pent_ind, point_num);

    //  Reorder neighbors in the clockwise order at the pole.
    reorder_neighbors_indx_pl(point_local, point_xyz, pent_ind, point_num);

    //  Produce rhombus' maps.
    maps = (int *)malloc((nl_region + 2) * (nl_region + 2) * nr * sizeof(int));
    produce_maps(maps, halo, nr, nl_region);

    //  Smooths the grid applying the spring dynamic method.
    if (sprd) {
        // Applies spring dynamics.
        spring_dynamics(point_local, pent_ind, glevel, spring_beta, point_xyz, point_num);

        //  Finds the q points.
        point_xyzq = (double *)malloc(6 * 3 * point_num * sizeof(double));

        find_qpoints(point_local, point_xyzq, point_xyz, pent_ind, point_num);

        // Fixes the position of the centroids.
        relocate_centres(point_local, point_xyzq, point_xyz, pent_ind, point_num);
    }
    else {
        //      Finds the q points.
        point_xyzq = (double *)malloc(6 * 3 * point_num * sizeof(double));
        find_qpoints(point_local, point_xyzq, point_xyz, pent_ind, point_num);
    }

    //  Recompute cartesians points for the new radius.
    correct_xyz_points(A, point_xyzq, point_xyz, pent_ind, point_num);

    //  Radial vectors with norm 1.
    func_r = (double *)malloc(3 * point_num * sizeof(double));
    compute_func(func_r, point_xyz, point_num);

    //  Compute control areas.
    areas   = (double *)malloc(6 * 3 * point_num * sizeof(double));
    areasTr = (double *)malloc(6 * point_num * sizeof(double));
    areasT  = (double *)malloc(point_num * sizeof(double));
    control_areas(
        areasT, areasTr, areas, point_local, point_xyzq, point_xyz, pent_ind, point_num, A);

    //  Computes control vectors.
    nvec   = (double *)malloc(6 * 3 * point_num * sizeof(double));
    nvecoa = (double *)malloc(6 * 3 * point_num * sizeof(double));
    nvecti = (double *)malloc(6 * 3 * point_num * sizeof(double));
    nvecte = (double *)malloc(6 * 3 * point_num * sizeof(double));
    mvec   = (double *)malloc(6 * 3 * point_num * sizeof(double));
    control_vec(nvec,
                nvecoa,
                nvecti,
                nvecte,
                areasT,
                point_xyz,
                point_xyzq,
                point_local,
                pent_ind,
                point_num,
                mvec);

    //  Set the Altitudes
    Altitude  = (double *)malloc(nv * sizeof(double));
    Altitudeh = (double *)malloc(nvi * sizeof(double));
    if (vert_refined) {
        // set_altitudes_refined(Altitude, Altitudeh, Top_altitude, nv, n_bl_layers);
        set_altitudes_softplus(
            Altitude, Altitudeh, Top_altitude, lowest_layer_thickness, transition_altitude, nv);
    }
    else {
        set_altitudes_uniform(Altitude, Altitudeh, Top_altitude, nv);
    }

    //  Converting to spherical coordinates.
    lonlat = (double *)malloc(2 * point_num * sizeof(double));
    cart2sphe(lonlat, point_xyz, point_num);

    //  Computes the divergence operator.
    div = (double *)malloc(7 * 3 * point_num * sizeof(double));
    div_operator(areasT, areas, div, nvec, pent_ind, point_num);

    //  Computes the gradient operator.
    grad = (double *)malloc(7 * 3 * point_num * sizeof(double));
    gra_operator(areasT, areas, grad, nvec, pent_ind, point_num);

    //  Computes the vertical curl operator.
    curlz = (double *)malloc(7 * 3 * point_num * sizeof(double));
    curlz_operator(areasT, areas, curlz, mvec, pent_ind, point_num);

    //  Computes zonal mean for sponge layer operations
    if (sponge == true) {
        zonal_mean_tab = (int *)malloc(3 * point_num * sizeof(int));
        zonal_mean_tab_f(zonal_mean_tab, lonlat, nlat, point_num, max_count);
    }

    log::printf(" GRID DONE!\n\n");
}


void Icogrid::sphere_ico(double *xyz_,
                         int     glevel,
                         int     n_region,
                         int     nl_region,
                         int     nl2,
                         int     kxl,
                         int     nfaces,
                         int *   pent_ind,
                         int     divide_face,
                         int     num) {

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

    double3 *xyz = (double3 *)xyz_;

    //  Local variables
    int sizei = pow(2.0, glevel);
    int count_points;
    int sizeip1 = sizei + 1;
    int sizei2  = sizeip1 * sizeip1;

    //  Temporary rhombus indexes.
    int *rhombi;
    rhombi = new int[10 * sizei2]();
    int *rhomb;
    rhomb = new int[10 * sizei * sizei]();

    //  Temporary main vertices.
    double3 *xyzi = new double3[(10 * sizei2 + 2)]();

    double w = 2.0 * acos(1.0 / (2.0 * sin(M_PI / 5.0)));

    //  First icosahedron coordinates (pentagons)
    //  Poles
    //  North
    xyzi[0] = make_double3(0.0, 0.0, 1.0);
    // South
    xyzi[1] = make_double3(0.0, 0.0, -1.0);
    //  Other points of the icosahedron.
    //  3
    xyzi[2]   = make_double3(cos(-M_PI / 5.0) * cos(M_PI / 2.0 - w),
                           sin(-M_PI / 5.0) * cos(M_PI / 2.0 - w),
                           sin(M_PI / 2.0 - w));
    rhombi[0] = 2;
    //  4
    xyzi[3]   = make_double3(cos(M_PI / 5.0) * cos(M_PI / 2.0 - w),
                           sin(M_PI / 5.0) * cos(M_PI / 2.0 - w),
                           sin(M_PI / 2.0 - w)),
    rhombi[1] = 3;
    //  5
    xyzi[4]   = make_double3(cos(3.0 * M_PI / 5.0) * cos(M_PI / 2.0 - w),
                           sin(3.0 * M_PI / 5.0) * cos(M_PI / 2.0 - w),
                           sin(M_PI / 2 - w));
    rhombi[2] = 4;
    //  6
    xyzi[5] = make_double3(
        cos(M_PI) * cos(M_PI / 2.0 - w), sin(M_PI) * cos(M_PI / 2.0 - w), sin(M_PI / 2.0 - w));
    rhombi[3] = 5;
    //  7
    xyzi[6]   = make_double3(cos(-(3.0 / 5.0) * M_PI) * cos(M_PI / 2.0 - w),
                           sin(-(3.0 / 5.0) * M_PI) * cos(M_PI / 2.0 - w),
                           sin(M_PI / 2.0 - w));
    rhombi[4] = 6;
    //  8
    xyzi[7] = make_double3(
        cos(0.0) * cos(w - M_PI / 2.0), sin(0.0) * cos(w - M_PI / 2.0), sin(w - M_PI / 2.0));
    rhombi[5] = 7;
    //  9
    xyzi[8]   = make_double3(cos(2.0 * M_PI / 5.0) * cos(w - M_PI / 2.0),
                           sin(2.0 * M_PI / 5.0) * cos(w - M_PI / 2.0),
                           sin(w - M_PI / 2.0));
    rhombi[6] = 8;
    //  10
    xyzi[9]   = make_double3(cos(4.0 * M_PI / 5.0) * cos(w - M_PI / 2.0),
                           sin(4.0 * M_PI / 5.0) * cos(w - M_PI / 2.0),
                           sin(w - M_PI / 2.0));
    rhombi[7] = 9;
    //  11
    xyzi[10]  = make_double3(cos(-4.0 * M_PI / 5.0) * cos(w - M_PI / 2.0),
                            sin(-4.0 * M_PI / 5.0) * cos(w - M_PI / 2.0),
                            sin(w - M_PI / 2.0));
    rhombi[8] = 10;
    //  12
    xyzi[11]  = make_double3(cos(-2.0 * M_PI / 5.0) * cos(w - M_PI / 2.0),
                            sin(-2.0 * M_PI / 5.0) * cos(w - M_PI / 2.0),
                            sin(w - M_PI / 2.0));
    rhombi[9] = 11;
    //
    //  Standard grid points.
    //
    int rhombi_points[10][4] = {{2, 0, 7, 3},
                                {3, 0, 8, 4},
                                {4, 0, 9, 5},
                                {5, 0, 10, 6},
                                {6, 0, 11, 2},
                                {7, 3, 1, 8},
                                {8, 4, 1, 9},
                                {9, 5, 1, 10},
                                {10, 6, 1, 11},
                                {11, 2, 1, 7}};

    for (int faces = 0; faces < 10; faces++) {
        rhombi[faces]                                     = rhombi_points[faces][0];
        rhombi[sizei * sizeip1 * 10 + faces]              = rhombi_points[faces][1];
        rhombi[sizei * 10 + faces]                        = rhombi_points[faces][2];
        rhombi[sizei * sizeip1 * 10 + sizei * 10 + faces] = rhombi_points[faces][3];
    }

    // Recursive method
    count_points = 12;
    for (int faces = 0; faces < 10; faces++) {
        for (int grd = 0; grd < glevel; grd++) {
            int ind_jump = pow(2.0, glevel - grd - 1);
            int ind      = ind_jump;
            for (int it = 0; it < pow(2.0, grd); it++) {
                int indx                                        = ind;
                int indy                                        = 0;
                rhombi[indx * sizeip1 * 10 + indy * 10 + faces] = count_points;
                int indx1                                       = ind + ind_jump;
                int indy1                                       = 0;
                int indx2                                       = ind - ind_jump;
                int indy2                                       = 0;

                int r_idx_1 = rhombi[indx1 * sizeip1 * 10 + indy1 * 10 + faces];
                int r_idx_2 = rhombi[indx2 * sizeip1 * 10 + indy2 * 10 + faces];

                xyzi[count_points] = normalize((xyzi[r_idx_1] + xyzi[r_idx_2]) * 0.5);
                count_points += 1;

                indx                                            = 0;
                indy                                            = ind;
                rhombi[indx * sizeip1 * 10 + indy * 10 + faces] = count_points;
                indx1                                           = 0;
                indy1                                           = ind + ind_jump;
                indx2                                           = 0;
                indy2                                           = ind - ind_jump;

                r_idx_1 = rhombi[indx1 * sizeip1 * 10 + indy1 * 10 + faces];
                r_idx_2 = rhombi[indx2 * sizeip1 * 10 + indy2 * 10 + faces];

                xyzi[count_points] = normalize((xyzi[r_idx_1] + xyzi[r_idx_2]) * 0.5);

                count_points += 1;

                indx                                            = ind;
                indy                                            = ind;
                rhombi[indx * sizeip1 * 10 + indy * 10 + faces] = count_points;
                indx1                                           = ind + ind_jump;
                indy1                                           = ind + ind_jump;
                indx2                                           = ind - ind_jump;
                indy2                                           = ind - ind_jump;
                r_idx_1 = rhombi[indx1 * sizeip1 * 10 + indy1 * 10 + faces];
                r_idx_2 = rhombi[indx2 * sizeip1 * 10 + indy2 * 10 + faces];

                xyzi[count_points] = normalize((xyzi[r_idx_1] + xyzi[r_idx_2]) * 0.5);

                count_points += 1;

                indx                                            = ind;
                indy                                            = sizei;
                rhombi[indx * sizeip1 * 10 + indy * 10 + faces] = count_points;
                indx1                                           = ind + ind_jump;
                indy1                                           = sizei;
                indx2                                           = ind - ind_jump;
                indy2                                           = sizei;
                r_idx_1 = rhombi[indx1 * sizeip1 * 10 + indy1 * 10 + faces];
                r_idx_2 = rhombi[indx2 * sizeip1 * 10 + indy2 * 10 + faces];

                xyzi[count_points] = normalize((xyzi[r_idx_1] + xyzi[r_idx_2]) * 0.5);

                count_points += 1;

                indx                                            = sizei;
                indy                                            = ind;
                rhombi[indx * sizeip1 * 10 + indy * 10 + faces] = count_points;
                indx1                                           = sizei;
                indy1                                           = ind + ind_jump;
                indx2                                           = sizei;
                indy2                                           = ind - ind_jump;

                r_idx_1 = rhombi[indx1 * sizeip1 * 10 + indy1 * 10 + faces];
                r_idx_2 = rhombi[indx2 * sizeip1 * 10 + indy2 * 10 + faces];

                xyzi[count_points] = normalize((xyzi[r_idx_1] + xyzi[r_idx_2]) * 0.5);

                count_points += 1;

                ind += 2 * ind_jump;
            }
        }
        for (int grd = 0; grd < glevel; grd++) {
            if (grd > 0) {
                // Upper triangle
                int ind_jump = pow(2.0, glevel - grd - 1);
                int indj     = ind_jump;
                for (int j = 0; j < pow(2.0, grd) - 1; j++) {
                    int indy = ind_jump;
                    int indx = 0;
                    for (int i = 0; i < j + 1; i++) {
                        int rindy = indy;
                        int rindx = sizei - indj + indx;

                        rhombi[rindx * sizeip1 * 10 + rindy * 10 + faces] = count_points;
                        int indx1                                         = rindx - ind_jump;
                        int indy1                                         = rindy - ind_jump;
                        int indx2                                         = rindx + ind_jump;
                        int indy2                                         = rindy + ind_jump;

                        int r_idx_1 = rhombi[indx1 * sizeip1 * 10 + indy1 * 10 + faces];
                        int r_idx_2 = rhombi[indx2 * sizeip1 * 10 + indy2 * 10 + faces];

                        xyzi[count_points] = normalize((xyzi[r_idx_1] + xyzi[r_idx_2]) * 0.5);

                        count_points += 1;

                        rhombi[(rindx - ind_jump) * sizeip1 * 10 + rindy * 10 + faces] =
                            count_points;
                        indx1   = rindx - ind_jump;
                        indy1   = rindy - ind_jump;
                        indx2   = rindx - ind_jump;
                        indy2   = rindy + ind_jump;
                        r_idx_1 = rhombi[indx1 * sizeip1 * 10 + indy1 * 10 + faces];
                        r_idx_2 = rhombi[indx2 * sizeip1 * 10 + indy2 * 10 + faces];

                        xyzi[count_points] = normalize((xyzi[r_idx_1] + xyzi[r_idx_2]) * 0.5);

                        count_points += 1;

                        rhombi[rindx * sizeip1 * 10 + (rindy + ind_jump) * 10 + faces] =
                            count_points;
                        indx1   = rindx + ind_jump;
                        indy1   = rindy + ind_jump;
                        indx2   = rindx - ind_jump;
                        indy2   = rindy + ind_jump;
                        r_idx_1 = rhombi[indx1 * sizeip1 * 10 + indy1 * 10 + faces];
                        r_idx_2 = rhombi[indx2 * sizeip1 * 10 + indy2 * 10 + faces];

                        xyzi[count_points] = normalize((xyzi[r_idx_1] + xyzi[r_idx_2]) * 0.5);

                        count_points += 1;

                        indx += 2 * ind_jump;
                        indy += 2 * ind_jump;
                    }
                    indj += 2 * ind_jump;
                }
            }
            if (grd > 0) {
                // Lower triangle
                int ind_jump = pow(2.0, glevel - grd - 1);
                int indj     = ind_jump;
                for (int j = 0; j < pow(2.0, grd) - 1; j++) {
                    int indy = 0;
                    int indx = ind_jump;
                    for (int i = 0; i < j + 1; i++) {

                        int rindx = indx;
                        int rindy = sizei - indj + indy;

                        rhombi[rindx * sizeip1 * 10 + rindy * 10 + faces] = count_points;
                        int indx1                                         = rindx - ind_jump;
                        int indy1                                         = rindy - ind_jump;
                        int indx2                                         = rindx + ind_jump;
                        int indy2                                         = rindy + ind_jump;

                        int r_idx_1 = rhombi[indx1 * sizeip1 * 10 + indy1 * 10 + faces];
                        int r_idx_2 = rhombi[indx2 * sizeip1 * 10 + indy2 * 10 + faces];

                        xyzi[count_points] = normalize((xyzi[r_idx_1] + xyzi[r_idx_2]) * 0.5);
                        count_points += 1;

                        rhombi[rindx * sizeip1 * 10 + (rindy - ind_jump) * 10 + faces] =
                            count_points;
                        indx1 = rindx - ind_jump;
                        indy1 = rindy - ind_jump;
                        indx2 = rindx + ind_jump;
                        indy2 = rindy - ind_jump;

                        r_idx_1 = rhombi[indx1 * sizeip1 * 10 + indy1 * 10 + faces];
                        r_idx_2 = rhombi[indx2 * sizeip1 * 10 + indy2 * 10 + faces];

                        xyzi[count_points] = normalize((xyzi[r_idx_1] + xyzi[r_idx_2]) * 0.5);

                        count_points += 1;

                        rhombi[(rindx + ind_jump) * sizeip1 * 10 + rindy * 10 + faces] =
                            count_points;
                        indx1   = rindx + ind_jump;
                        indy1   = rindy + ind_jump;
                        indx2   = rindx + ind_jump;
                        indy2   = rindy - ind_jump;
                        r_idx_1 = rhombi[indx1 * sizeip1 * 10 + indy1 * 10 + faces];
                        r_idx_2 = rhombi[indx2 * sizeip1 * 10 + indy2 * 10 + faces];

                        xyzi[count_points] = normalize((xyzi[r_idx_1] + xyzi[r_idx_2]) * 0.5);

                        count_points += 1;

                        indx += 2 * ind_jump;
                        indy += 2 * ind_jump;
                    }
                    indj += 2 * ind_jump;
                }
            }
        }
    }

    int nli_region = sqrt(n_region * 1.0);

    for (int fc = 0; fc < 10; fc++)
        for (int j = 0; j < sizei; j++)
            for (int i = 0; i < sizei; i++)
                rhomb[i * sizei * 10 + j * 10 + fc] = rhombi[i * sizeip1 * 10 + j * 10 + fc];

    for (int fc = 0; fc < 10; fc++)
        for (int kx = 0; kx < kxl; kx++)
            for (int ky = 0; ky < kxl; ky++)
                for (int i = 0; i < nl_region; i++)
                    for (int j = 0; j < nl_region; j++) {
                        int idx1 =
                            fc * nfaces * nl2 + ky * kxl * nl2 + kx * nl2 + j * nl_region + i;
                        int idx2 =
                            rhomb[((ky * nl_region + j) * nli_region + kx * nl_region + i) * 10
                                  + fc];

                        xyz[idx1] = xyzi[idx2];
                    }


    //North
    xyz[num - 2] = make_double3(0.0, 0.0, 1.0);

    //South
    xyz[num - 1] = make_double3(0.0, 0.0, -1.0);

    //Pentagons' indexes
    for (int faces = 0; faces < 10; faces++)
        pent_ind[faces] = faces * nfaces * nl2;
    pent_ind[10] = num - 2;
    pent_ind[11] = num - 1;

    delete[] rhombi;
    delete[] rhomb;
    delete[] xyzi;
}

void Icogrid::neighbors_indx(int *point_local, double *xyz, int *pent_ind, int point_num) {

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
    int     position;
    double  small;
    double *distance;
    distance = new double[point_num]();

    for (int i = 0; i < 6 * point_num; i++)
        point_local[i] = -1;

    //  Find neighbors.
    for (int i = 0; i < point_num; i++) {
        for (int j = 0; j < point_num; j++) {
            distance[j] = sqrt(pow(xyz[0 + j * 3] - xyz[0 + i * 3], 2)
                               + pow(xyz[1 + j * 3] - xyz[1 + i * 3], 2)
                               + pow(xyz[2 + j * 3] - xyz[2 + i * 3], 2));
        }
        for (int k = 0; k < 6; k++) { // Hexagons.
            small    = 2;
            position = -1;
            for (int j = 0; j < point_num; j++) {
                if (k == 0) {
                    if (small >= distance[j] && distance[j] > 0) {
                        small    = distance[j];
                        position = j;
                    }
                }
                else {
                    if (small >= distance[j] && distance[j] > 0) {
                        if (point_local[i * 6 + 0] != j && point_local[i * 6 + 1] != j
                            && point_local[i * 6 + 2] != j && point_local[i * 6 + 3] != j
                            && point_local[i * 6 + 4] != j && point_local[i * 6 + 5] != j) {
                            small    = distance[j];
                            position = j;
                        }
                    }
                }
            }
            point_local[i * 6 + k] = position;
        }
        for (int k = 0; k < 12; k++) {
            if (i == pent_ind[k]) { // Pentagons.
                point_local[i * 6 + 5] = -1;
            }
        }
    }

    delete[] distance;
}

void Icogrid::neighbors_indx_pl(int *point_local, double *xyz, int *pent_ind, int point_num) {

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
    int     position;
    double  small;
    double *distance;
    distance = new double[point_num]();

    //  Find neighbors.
    for (int i = point_num - 2; i < point_num; i++) {
        for (int j = 0; j < point_num; j++) {
            distance[j] = sqrt(pow(xyz[0 + j * 3] - xyz[0 + i * 3], 2)
                               + pow(xyz[1 + j * 3] - xyz[1 + i * 3], 2)
                               + pow(xyz[2 + j * 3] - xyz[2 + i * 3], 2));
        }
        for (int k = 0; k < 5; k++) {
            small    = 2;
            position = -1;
            for (int j = 0; j < point_num; j++) {
                if (k == 0) {
                    if (small >= distance[j] && distance[j] > 0) {
                        small    = distance[j];
                        position = j;
                    }
                }
                else {
                    if (small >= distance[j] && distance[j] > 0) {
                        if (point_local[i * 6 + 0] != j && point_local[i * 6 + 1] != j
                            && point_local[i * 6 + 2] != j && point_local[i * 6 + 3] != j
                            && point_local[i * 6 + 4] != j && point_local[i * 6 + 5] != j) {
                            small    = distance[j];
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

void Icogrid::reorder_neighbors_indx(int *point_local, double *xyz, int *pent_ind, int point_num) {

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
    int    geo;

    //  Temporary arrays.
    double *ang;
    int *   point_local_new;
    ang             = new double[6]();
    point_local_new = new int[point_num * 6];

    for (int i = 0; i < point_num; i++)
        for (int j = 0; j < 6; j++)
            point_local_new[i * 6 + j] = point_local[i * 6 + j];

    // Reorder indexes.
    for (int i = 0; i < point_num; i++) {
        geo = 6; // Hexagons.
        for (int k = 0; k < 12; k++)
            if (i == pent_ind[k])
                geo = 5; // Pentagons.
        if (geo == 5) {
            for (int j = 0; j < 6; j++)
                ang[j] = 0;
            x1 = 0;
            y1 = 0;
            z1 = 0;
            for (int j = 0; j < 5; j++) { // Central coplanar point.
                x1 = x1 + xyz[point_local[i * 6 + j] * 3 + 0] / 5;
                y1 = y1 + xyz[point_local[i * 6 + j] * 3 + 1] / 5;
                z1 = z1 + xyz[point_local[i * 6 + j] * 3 + 2] / 5;
            }
            dx21 = xyz[point_local[i * 6 + 0] * 3 + 0] - x1;
            dy21 = xyz[point_local[i * 6 + 0] * 3 + 1] - y1;
            dz21 = xyz[point_local[i * 6 + 0] * 3 + 2] - z1;
            for (int j = 1; j < 5; j++) {
                dx31  = xyz[point_local[i * 6 + j] * 3 + 0] - x1;
                dy31  = xyz[point_local[i * 6 + j] * 3 + 1] - y1;
                dz31  = xyz[point_local[i * 6 + j] * 3 + 2] - z1;
                n21   = sqrt(dx21 * dx21 + dy21 * dy21 + dz21 * dz21);
                n31   = sqrt(dx31 * dx31 + dy31 * dy31 + dz31 * dz31);
                dircx = dy21 * dz31 - dy31 * dz21;
                dircy = dz21 * dx31 - dz31 * dx21;
                dircz = dx21 * dy31 - dx31 * dy21;
                dircn = dircx * x1 + dircy * y1 + dircz * z1;
                if (dircn <= 0) {
                    ang[j] = acos((dx21 * dx31 + dy21 * dy31 + dz21 * dz31) / (n21 * n31));
                }
                else {
                    ang[j] =
                        2.0 * M_PI - acos((dx21 * dx31 + dy21 * dy31 + dz21 * dz31) / (n21 * n31));
                }
            }
            for (int j = 0; j < 5; j++) {
                for (int k = 4; k >= j; k--) {
                    if (ang[j] > ang[k]) {
                        temp                       = ang[k];
                        tempp                      = point_local_new[i * 6 + k];
                        ang[k]                     = ang[j];
                        point_local_new[i * 6 + k] = point_local_new[i * 6 + j];
                        ang[j]                     = temp;
                        point_local_new[i * 6 + j] = tempp;
                    }
                }
            }
        }
        else {
            for (int j = 0; j < 6; j++)
                ang[j] = 0;
            x1 = 0;
            y1 = 0;
            z1 = 0;
            for (int j = 0; j < 6; j++) { // Central coplanar point.
                x1 = x1 + xyz[point_local[i * 6 + j] * 3 + 0] / 6;
                y1 = y1 + xyz[point_local[i * 6 + j] * 3 + 1] / 6;
                z1 = z1 + xyz[point_local[i * 6 + j] * 3 + 2] / 6;
            }
            dx21 = xyz[point_local[i * 6 + 0] * 3 + 0] - x1;
            dy21 = xyz[point_local[i * 6 + 0] * 3 + 1] - y1;
            dz21 = xyz[point_local[i * 6 + 0] * 3 + 2] - z1;
            for (int j = 1; j < 6; j++) {
                dx31  = xyz[point_local[i * 6 + j] * 3 + 0] - x1;
                dy31  = xyz[point_local[i * 6 + j] * 3 + 1] - y1;
                dz31  = xyz[point_local[i * 6 + j] * 3 + 2] - z1;
                n21   = sqrt(dx21 * dx21 + dy21 * dy21 + dz21 * dz21);
                n31   = sqrt(dx31 * dx31 + dy31 * dy31 + dz31 * dz31);
                dircx = dy21 * dz31 - dy31 * dz21;
                dircy = dz21 * dx31 - dz31 * dx21;
                dircz = dx21 * dy31 - dx31 * dy21;
                dircn = dircx * x1 + dircy * y1 + dircz * z1;
                if (dircn <= 0)
                    ang[j] = acos((dx21 * dx31 + dy21 * dy31 + dz21 * dz31) / (n21 * n31));
                else
                    ang[j] =
                        2.0 * M_PI - acos((dx21 * dx31 + dy21 * dy31 + dz21 * dz31) / (n21 * n31));
            }
            for (int j = 0; j < 6; j++) {
                for (int k = 5; k >= j; k--) {
                    if (ang[j] > ang[k]) {
                        temp                       = ang[k];
                        tempp                      = point_local_new[i * 6 + k];
                        ang[k]                     = ang[j];
                        point_local_new[i * 6 + k] = point_local_new[i * 6 + j];
                        ang[j]                     = temp;
                        point_local_new[i * 6 + j] = tempp;
                    }
                }
            }
        }
    }

    for (int i = 0; i < point_num; i++)
        for (int j = 0; j < 6; j++)
            point_local[i * 6 + j] = point_local_new[i * 6 + j];

    delete[] ang;
    delete[] point_local_new;
}

void Icogrid::reorder_neighbors_indx_pl(int *   point_local,
                                        double *xyz,
                                        int *   pent_ind,
                                        int     point_num) {

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
    int *   point_local_new;
    ang             = new double[6]();
    point_local_new = new int[point_num * 6];

    for (int i = 0; i < point_num; i++)
        for (int j = 0; j < 6; j++)
            point_local_new[i * 6 + j] = point_local[i * 6 + j];

    // Reorder indexes.
    for (int i = point_num - 2; i < point_num; i++) {
        for (int j = 0; j < 6; j++)
            ang[j] = 0;
        x1 = 0;
        y1 = 0;
        z1 = 0;
        for (int j = 0; j < 5; j++) { // Central coplanar point.
            x1 = x1 + xyz[point_local[i * 6 + j] * 3 + 0] / 5;
            y1 = y1 + xyz[point_local[i * 6 + j] * 3 + 1] / 5;
            z1 = z1 + xyz[point_local[i * 6 + j] * 3 + 2] / 5;
        }
        dx21 = xyz[point_local[i * 6 + 0] * 3 + 0] - x1;
        dy21 = xyz[point_local[i * 6 + 0] * 3 + 1] - y1;
        dz21 = xyz[point_local[i * 6 + 0] * 3 + 2] - z1;
        for (int j = 1; j < 5; j++) {
            dx31  = xyz[point_local[i * 6 + j] * 3 + 0] - x1;
            dy31  = xyz[point_local[i * 6 + j] * 3 + 1] - y1;
            dz31  = xyz[point_local[i * 6 + j] * 3 + 2] - z1;
            n21   = sqrt(dx21 * dx21 + dy21 * dy21 + dz21 * dz21);
            n31   = sqrt(dx31 * dx31 + dy31 * dy31 + dz31 * dz31);
            dircx = dy21 * dz31 - dy31 * dz21;
            dircy = dz21 * dx31 - dz31 * dx21;
            dircz = dx21 * dy31 - dx31 * dy21;
            dircn = dircx * x1 + dircy * y1 + dircz * z1;
            if (dircn <= 0) {
                ang[j] = acos((dx21 * dx31 + dy21 * dy31 + dz21 * dz31) / (n21 * n31));
            }
            else {
                ang[j] = 2.0 * M_PI - acos((dx21 * dx31 + dy21 * dy31 + dz21 * dz31) / (n21 * n31));
            }
        }
        for (int j = 0; j < 5; j++) {
            for (int k = 4; k >= j; k--) {
                if (ang[j] > ang[k]) {
                    temp                       = ang[k];
                    tempp                      = point_local_new[i * 6 + k];
                    ang[k]                     = ang[j];
                    point_local_new[i * 6 + k] = point_local_new[i * 6 + j];
                    ang[j]                     = temp;
                    point_local_new[i * 6 + j] = tempp;
                }
            }
        }
    }

    for (int i = point_num - 2; i < point_num; i++)
        for (int j = 0; j < 6; j++)
            point_local[i * 6 + j] = point_local_new[i * 6 + j];

    delete[] ang;
    delete[] point_local_new;
}

void Icogrid::generate_halos(int *halo, int *point_local, int n_region, int divide_face) {

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
    int nl_region = sqrt(n_region * 1.0);
    nl_region     = nl_region / pow(2.0, divide_face);
    int nl2       = nl_region * nl_region;
    int nfaces    = pow(4.0, divide_face);
    int kxl       = sqrt(nfaces * 1.0);
    int nh        = 4 * (nl_region + 2);

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

    for (int i = 0; i < 10 * nfaces * 4 * (nl_region + 2); i++)
        halo[i] = -1;

    // Generate halos.
    for (int fc = 0; fc < 10; fc++) {
        for (int kx = 0; kx < kxl; kx++) {
            for (int ky = 0; ky < kxl; ky++) {
                if (kx == 0 && ky == 0) { // First point is a pentagon
                    // Side 1
                    for (int i = 0; i < nl_region - 1; i++) {
                        int j = 0;
                        index_value1 =
                            fc * nfaces * nl2 + ky * kxl * nl2 + kx * nl2 + j * nl_region + i;
                        index_value2 =
                            fc * nfaces * nl2 + ky * kxl * nl2 + kx * nl2 + j * nl_region + i + 1;
                        if (kx == 0 && ky == 0 && i == 0 && j == 0) {
                            for (int k = 0; k < 5; k++) {
                                for (int kk = 0; kk < 6; kk++) {
                                    v1 = point_local[index_value1 * 6 + k];
                                    v2 = point_local[index_value2 * 6 + kk];
                                    if (v1 == v2) {
                                        out_ind = 0;
                                        for (int kkk = 0; kkk < nl2; kkk++)
                                            if (v1
                                                == fc * nfaces * nl2 + ky * kxl * nl2 + kx * nl2
                                                       + kkk)
                                                out_ind = 1;
                                        if (out_ind == 0)
                                            halo[fc * kxl * kxl * nh + ky * kxl * nh + kx * nh
                                                 + i] = v1;
                                    }
                                }
                            }
                        }
                        else {
                            for (int k = 0; k < 6; k++) {
                                for (int kk = 0; kk < 6; kk++) {
                                    v1 = point_local[index_value1 * 6 + k];
                                    v2 = point_local[index_value2 * 6 + kk];
                                    if (v1 == v2) {
                                        out_ind = 0;
                                        for (int kkk = 0; kkk < nl2; kkk++)
                                            if (v1
                                                == fc * nfaces * nl2 + ky * kxl * nl2 + kx * nl2
                                                       + kkk)
                                                out_ind = 1;
                                        if (out_ind == 0)
                                            halo[fc * kxl * kxl * nh + ky * kxl * nh + kx * nh
                                                 + i] = v1;
                                    }
                                }
                            }
                        }
                    }
                    int i = nl_region - 1;
                    int j = 0;
                    index_value1 =
                        fc * nfaces * nl2 + ky * kxl * nl2 + kx * nl2 + j * nl_region + i;
                    for (int k = 0; k < 6; k++) {
                        v1 = point_local[index_value1 * 6 + k];
                        if (v1 == halo[fc * kxl * kxl * nh + ky * kxl * nh + kx * nh + i - 1])
                            ind_nei = k;
                    }
                    if (ind_nei > 0)
                        halo[fc * kxl * kxl * nh + ky * kxl * nh + kx * nh + i] =
                            point_local[index_value1 * 6 + ind_nei - 1];
                    else
                        halo[fc * kxl * kxl * nh + ky * kxl * nh + kx * nh + i] =
                            point_local[index_value1 * 6 + 5];

                    //Side 2
                    i = nl_region - 1;
                    j = 0;
                    index_value1 =
                        fc * nfaces * nl2 + ky * kxl * nl2 + kx * nl2 + j * nl_region + i;
                    for (int k = 0; k < 6; k++) {
                        v1 = point_local[index_value1 * 6 + k];
                        if (v1
                            == halo[fc * kxl * kxl * nh + ky * kxl * nh + kx * nh + nl_region - 1])
                            ind_nei = k;
                    }
                    if (ind_nei > 0)
                        halo[fc * kxl * kxl * nh + ky * kxl * nh + kx * nh + nl_region + 2] =
                            point_local[index_value1 * 6 + ind_nei - 1];
                    else
                        halo[fc * kxl * kxl * nh + ky * kxl * nh + kx * nh + nl_region + 2] =
                            point_local[index_value1 * 6 + 5];

                    i = nl_region - 1;
                    for (int j = 0; j < nl_region - 1; j++) {
                        index_value1 =
                            fc * nfaces * nl2 + ky * kxl * nl2 + kx * nl2 + j * nl_region + i;
                        index_value2 =
                            fc * nfaces * nl2 + ky * kxl * nl2 + kx * nl2 + (j + 1) * nl_region + i;
                        for (int k = 0; k < 6; k++) {
                            for (int kk = 0; kk < 6; kk++) {
                                v1 = point_local[index_value1 * 6 + k];
                                v2 = point_local[index_value2 * 6 + kk];
                                if (v1 == v2) {
                                    int out_ind = 0;
                                    for (int kkk = 0; kkk < nl2; kkk++)
                                        if (v1
                                            == fc * nfaces * nl2 + ky * kxl * nl2 + kx * nl2 + kkk)
                                            out_ind = 1;
                                    if (out_ind == 0)
                                        halo[fc * kxl * kxl * nh + ky * kxl * nh + kx * nh
                                             + (nl_region + 2) + j + 1] = v1;
                                }
                            }
                        }
                    }

                    // Side 3
                    i = nl_region - 1;
                    j = nl_region - 1;
                    index_value1 =
                        fc * nfaces * nl2 + ky * kxl * nl2 + kx * nl2 + j * nl_region + i;
                    for (int k = 0; k < 6; k++) {
                        v1 = point_local[index_value1 * 6 + k];
                        if (v1
                            == halo[fc * kxl * kxl * nh + ky * kxl * nh + kx * nh + nl_region + 2
                                    + nl_region - 1])
                            ind_nei = k;
                    }
                    if (ind_nei > 0)
                        halo[fc * kxl * kxl * nh + ky * kxl * nh + kx * nh + 2 * (nl_region + 2)] =
                            point_local[index_value1 * 6 + ind_nei - 1];
                    else
                        halo[fc * kxl * kxl * nh + ky * kxl * nh + kx * nh + 2 * (nl_region + 2)] =
                            point_local[index_value1 * 6 + 5];

                    for (int i = 0; i < nl_region - 1; i++) {
                        j            = nl_region - 1;
                        index_value1 = fc * nfaces * nl2 + ky * kxl * nl2 + kx * nl2 + j * nl_region
                                       + nl_region - i - 1;
                        index_value2 = fc * nfaces * nl2 + ky * kxl * nl2 + kx * nl2 + j * nl_region
                                       + nl_region - i - 2;
                        for (int k = 0; k < 6; k++) {
                            for (int kk = 0; kk < 6; kk++) {
                                v1 = point_local[index_value1 * 6 + k];
                                v2 = point_local[index_value2 * 6 + kk];

                                if (v1 == v2) {
                                    int out_ind = 0;
                                    for (int kkk = 0; kkk < nl2; kkk++)
                                        if (v1
                                            == fc * nfaces * nl2 + ky * kxl * nl2 + kx * nl2 + kkk)
                                            out_ind = 1;
                                    if (out_ind == 0)
                                        halo[fc * kxl * kxl * nh + ky * kxl * nh + kx * nh
                                             + 2 * (nl_region + 2) + i + 1] = v1;
                                }
                            }
                        }
                    }

                    i = 0;
                    j = nl_region - 1;
                    index_value1 =
                        fc * nfaces * nl2 + ky * kxl * nl2 + kx * nl2 + j * nl_region + i;
                    for (int k = 0; k < 6; k++) {
                        v1 = point_local[index_value1 * 6 + k];
                        if (v1
                            == halo[fc * kxl * kxl * nh + ky * kxl * nh + kx * nh
                                    + 2 * (nl_region + 2) + nl_region - 1])
                            ind_nei = k;
                    }
                    if (ind_nei > 0)
                        halo[fc * kxl * kxl * nh + ky * kxl * nh + kx * nh + 2 * (nl_region + 2)
                             + nl_region] = point_local[index_value1 * 6 + ind_nei - 1];
                    else
                        halo[fc * kxl * kxl * nh + ky * kxl * nh + kx * nh + 2 * (nl_region + 2)
                             + nl_region] = point_local[index_value1 * 6 + 5];

                    // Side 4
                    i = 0;
                    j = nl_region - 1;
                    index_value1 =
                        fc * nfaces * nl2 + ky * kxl * nl2 + kx * nl2 + j * nl_region + i;
                    for (int k = 0; k < 6; k++) {
                        v1 = point_local[index_value1 * 6 + k];
                        if (v1
                            == halo[fc * kxl * kxl * nh + ky * kxl * nh + kx * nh
                                    + 2 * (nl_region + 2) + nl_region])
                            ind_nei = k;
                    }
                    if (ind_nei > 0)
                        halo[fc * kxl * kxl * nh + ky * kxl * nh + kx * nh + 3 * (nl_region + 2)] =
                            point_local[index_value1 * 6 + ind_nei - 1];
                    else
                        halo[fc * kxl * kxl * nh + ky * kxl * nh + kx * nh + 3 * (nl_region + 2)] =
                            point_local[index_value1 * 6 + 5];

                    i = 0;
                    for (int j = 0; j < nl_region - 1; j++) {
                        index_value1 = fc * nfaces * nl2 + ky * kxl * nl2 + kx * nl2
                                       + (nl_region - j - 1) * nl_region;
                        index_value2 = fc * nfaces * nl2 + ky * kxl * nl2 + kx * nl2
                                       + (nl_region - j - 2) * nl_region;
                        for (int k = 0; k < 6; k++) {
                            for (int kk = 0; kk < 6; kk++) {
                                v1 = point_local[index_value1 * 6 + k];
                                v2 = point_local[index_value2 * 6 + kk];
                                if (v1 == v2) {
                                    out_ind = 0;
                                    for (int kkk = 0; kkk < nl2; kkk++)
                                        if (v1
                                            == fc * nfaces * nl2 + ky * kxl * nl2 + kx * nl2 + kkk)
                                            out_ind = 1;
                                    if (out_ind == 0)
                                        halo[fc * kxl * kxl * nh + ky * kxl * nh + kx * nh
                                             + 3 * (nl_region + 2) + j + 1] = v1;
                                }
                            }
                        }
                    }
                }
                else { //Hexagon
                    for (int i = 0; i < nl_region - 1; i++) {
                        int j = 0;
                        index_value1 =
                            fc * nfaces * nl2 + ky * kxl * nl2 + kx * nl2 + j * nl_region + i;
                        index_value2 =
                            fc * nfaces * nl2 + ky * kxl * nl2 + kx * nl2 + j * nl_region + i + 1;
                        if (kx == 0 && ky == 0 && i == 0 && j == 0) {
                            for (int k = 0; k < 6; k++) {
                                for (int kk = 0; kk < 6; kk++) {
                                    v1 = point_local[index_value1 * 6 + k];
                                    v2 = point_local[index_value2 * 6 + kk];
                                    if (v1 == v2) {
                                        int out_ind = 0;
                                        for (int kkk = 0; kkk < nl2; kkk++)
                                            if (v1
                                                == fc * nfaces * nl2 + ky * kxl * nl2 + kx * nl2
                                                       + kkk)
                                                out_ind = 1;
                                        if (out_ind == 0)
                                            halo[fc * kxl * kxl * nh + ky * kxl * nh + kx * nh
                                                 + i] = v1;
                                    }
                                }
                            }
                        }
                        else {
                            for (int k = 0; k < 6; k++) {
                                for (int kk = 0; kk < 6; kk++) {
                                    v1 = point_local[index_value1 * 6 + k];
                                    v2 = point_local[index_value2 * 6 + kk];
                                    if (v1 == v2) {
                                        out_ind = 0;
                                        for (int kkk = 0; kkk < nl2; kkk++)
                                            if (v1
                                                == fc * nfaces * nl2 + ky * kxl * nl2 + kx * nl2
                                                       + kkk)
                                                out_ind = 1;
                                        if (out_ind == 0)
                                            halo[fc * kxl * kxl * nh + ky * kxl * nh + kx * nh
                                                 + i] = v1;
                                    }
                                }
                            }
                        }
                    }
                    int i = nl_region - 1;
                    int j = 0;
                    index_value1 =
                        fc * nfaces * nl2 + ky * kxl * nl2 + kx * nl2 + j * nl_region + i;
                    for (int k = 0; k < 6; k++) {
                        v1 = point_local[index_value1 * 6 + k];
                        if (v1 == halo[fc * kxl * kxl * nh + ky * kxl * nh + kx * nh + i - 1])
                            ind_nei = k;
                    }
                    if (ind_nei > 0)
                        halo[fc * kxl * kxl * nh + ky * kxl * nh + kx * nh + i] =
                            point_local[index_value1 * 6 + ind_nei - 1];
                    else
                        halo[fc * kxl * kxl * nh + ky * kxl * nh + kx * nh + i] =
                            point_local[index_value1 * 6 + 5];

                    //Side 2
                    i = nl_region - 1;
                    j = 0;
                    index_value1 =
                        fc * nfaces * nl2 + ky * kxl * nl2 + kx * nl2 + j * nl_region + i;
                    for (int k = 0; k < 6; k++) {
                        v1 = point_local[index_value1 * 6 + k];
                        if (v1
                            == halo[fc * kxl * kxl * nh + ky * kxl * nh + kx * nh + nl_region - 1])
                            ind_nei = k;
                    }
                    if (ind_nei > 0)
                        halo[fc * kxl * kxl * nh + ky * kxl * nh + kx * nh + nl_region + 2] =
                            point_local[index_value1 * 6 + ind_nei - 1];
                    else
                        halo[fc * kxl * kxl * nh + ky * kxl * nh + kx * nh + nl_region + 2] =
                            point_local[index_value1 * 6 + 5];

                    i = nl_region - 1;
                    for (int j = 0; j < nl_region - 1; j++) {
                        index_value1 =
                            fc * nfaces * nl2 + ky * kxl * nl2 + kx * nl2 + j * nl_region + i;
                        index_value2 =
                            fc * nfaces * nl2 + ky * kxl * nl2 + kx * nl2 + (j + 1) * nl_region + i;
                        for (int k = 0; k < 6; k++) {
                            for (int kk = 0; kk < 6; kk++) {
                                v1 = point_local[index_value1 * 6 + k];
                                v2 = point_local[index_value2 * 6 + kk];
                                if (v1 == v2) {
                                    out_ind = 0;
                                    for (int kkk = 0; kkk < nl2; kkk++)
                                        if (v1
                                            == fc * nfaces * nl2 + ky * kxl * nl2 + kx * nl2 + kkk)
                                            out_ind = 1;
                                    if (out_ind == 0)
                                        halo[fc * kxl * kxl * nh + ky * kxl * nh + kx * nh
                                             + (nl_region + 2) + j + 1] = v1;
                                }
                            }
                        }
                    }

                    // Side 3

                    i = nl_region - 1;
                    j = nl_region - 1;
                    index_value1 =
                        fc * nfaces * nl2 + ky * kxl * nl2 + kx * nl2 + j * nl_region + i;
                    for (int k = 0; k < 6; k++) {
                        v1 = point_local[index_value1 * 6 + k];
                        if (v1
                            == halo[fc * kxl * kxl * nh + ky * kxl * nh + kx * nh + nl_region + 2
                                    + nl_region - 1])
                            ind_nei = k;
                    }
                    if (ind_nei > 0)
                        halo[fc * kxl * kxl * nh + ky * kxl * nh + kx * nh + 2 * (nl_region + 2)] =
                            point_local[index_value1 * 6 + ind_nei - 1];
                    else
                        halo[fc * kxl * kxl * nh + ky * kxl * nh + kx * nh + 2 * (nl_region + 2)] =
                            point_local[index_value1 * 6 + 5];

                    for (int i = 0; i < nl_region - 1; i++) {
                        j            = nl_region - 1;
                        index_value1 = fc * nfaces * nl2 + ky * kxl * nl2 + kx * nl2 + j * nl_region
                                       + nl_region - i - 1;
                        index_value2 = fc * nfaces * nl2 + ky * kxl * nl2 + kx * nl2 + j * nl_region
                                       + nl_region - i - 2;
                        for (int k = 0; k < 6; k++) {
                            for (int kk = 0; kk < 6; kk++) {
                                v1 = point_local[index_value1 * 6 + k];
                                v2 = point_local[index_value2 * 6 + kk];
                                if (v1 == v2) {
                                    out_ind = 0;
                                    for (int kkk = 0; kkk < nl2; kkk++)
                                        if (v1
                                            == fc * nfaces * nl2 + ky * kxl * nl2 + kx * nl2 + kkk)
                                            out_ind = 1;
                                    if (out_ind == 0)
                                        halo[fc * kxl * kxl * nh + ky * kxl * nh + kx * nh
                                             + 2 * (nl_region + 2) + i + 1] = v1;
                                }
                            }
                        }
                    }

                    i = 0;
                    j = nl_region - 1;
                    index_value1 =
                        fc * nfaces * nl2 + ky * kxl * nl2 + kx * nl2 + j * nl_region + i;
                    for (int k = 0; k < 6; k++) {
                        v1 = point_local[index_value1 * 6 + k];
                        if (v1
                            == halo[fc * kxl * kxl * nh + ky * kxl * nh + kx * nh
                                    + 2 * (nl_region + 2) + nl_region - 1])
                            ind_nei = k;
                    }
                    if (ind_nei > 0)
                        halo[fc * kxl * kxl * nh + ky * kxl * nh + kx * nh + 2 * (nl_region + 2)
                             + nl_region] = point_local[index_value1 * 6 + ind_nei - 1];
                    else
                        halo[fc * kxl * kxl * nh + ky * kxl * nh + kx * nh + 2 * (nl_region + 2)
                             + nl_region] = point_local[index_value1 * 6 + 5];

                    // Side 4
                    i = 0;
                    j = nl_region - 1;
                    index_value1 =
                        fc * nfaces * nl2 + ky * kxl * nl2 + kx * nl2 + j * nl_region + i;
                    for (int k = 0; k < 6; k++) {
                        v1 = point_local[index_value1 * 6 + k];
                        if (v1
                            == halo[fc * kxl * kxl * nh + ky * kxl * nh + kx * nh
                                    + 2 * (nl_region + 2) + nl_region])
                            ind_nei = k;
                    }
                    if (ind_nei > 0)
                        halo[fc * kxl * kxl * nh + ky * kxl * nh + kx * nh + 3 * (nl_region + 2)] =
                            point_local[index_value1 * 6 + ind_nei - 1];
                    else
                        halo[fc * kxl * kxl * nh + ky * kxl * nh + kx * nh + 3 * (nl_region + 2)] =
                            point_local[index_value1 * 6 + 5];

                    i = 0;
                    for (int j = 0; j < nl_region - 1; j++) {
                        index_value1 = fc * nfaces * nl2 + ky * kxl * nl2 + kx * nl2
                                       + (nl_region - j - 1) * nl_region;
                        index_value2 = fc * nfaces * nl2 + ky * kxl * nl2 + kx * nl2
                                       + (nl_region - j - 2) * nl_region;
                        for (int k = 0; k < 6; k++) {
                            for (int kk = 0; kk < 6; kk++) {
                                v1 = point_local[index_value1 * 6 + k];
                                v2 = point_local[index_value2 * 6 + kk];
                                if (v1 == v2) {
                                    out_ind = 0;
                                    for (int kkk = 0; kkk < nl2; kkk++)
                                        if (v1
                                            == fc * nfaces * nl2 + ky * kxl * nl2 + kx * nl2 + kkk)
                                            out_ind = 1;
                                    if (out_ind == 0)
                                        halo[fc * kxl * kxl * nh + ky * kxl * nh + kx * nh
                                             + 3 * (nl_region + 2) + j + 1] = v1;
                                }
                            }
                        }
                    }

                    i = 0;
                    j = 0;
                    index_value1 =
                        fc * nfaces * nl2 + ky * kxl * nl2 + kx * nl2 + j * nl_region + i;
                    for (int k = 0; k < 6; k++) {
                        v1 = point_local[index_value1 * 6 + k];
                        if (v1
                            == halo[fc * kxl * kxl * nh + ky * kxl * nh + kx * nh
                                    + 3 * (nl_region + 2) + nl_region - 1])
                            ind_nei = k;
                    }
                    if (ind_nei > 0)
                        halo[fc * kxl * kxl * nh + ky * kxl * nh + kx * nh + 3 * (nl_region + 2)
                             + nl_region] = point_local[index_value1 * 6 + ind_nei - 1];
                    else
                        halo[fc * kxl * kxl * nh + ky * kxl * nh + kx * nh + 3 * (nl_region + 2)
                             + nl_region] = point_local[index_value1 * 6 + 5];
                }
            }
        }
    }
}

void Icogrid::reorder_neighbors_indx_rhombi(int *point_local,
                                            int *halo,
                                            int *pent_ind,
                                            int  nl_region,
                                            int  nl2,
                                            int  ni,
                                            int  nlhalo,
                                            int  point_num) {

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
    int nh2 = nlhalo * 4;
    for (int i = 0; i < point_num; i++)
        for (int j = 0; j < 6; j++)
            point_local[i * 6 + j] = -1;

    for (int i = 0; i < ni; i++) {

        // Corner
        geo = 6; // Hexagons.
        np  = nl2 * i;
        for (int l = 0; l < 12; l++)
            if (np == pent_ind[l])
                geo = 5; // Pentagons.

        // Find neibours.
        for (int j = 0; j < nl_region; j++) {
            for (int k = 0; k < nl_region; k++) {
                np = nl2 * i + j * nl_region + k;
                // Right
                if (j == 0) {
                    point_local[np * 6 + 0] = nl2 * i + (j + 1) * nl_region + k;
                    point_local[np * 6 + 3] = halo[i * nh2 + k];
                }
                else if (j == nl_region - 1) {
                    point_local[np * 6 + 0] = halo[i * nh2 + 2 * nlhalo + nl_region - k];
                    point_local[np * 6 + 3] = nl2 * i + (j - 1) * nl_region + k;
                }
                else {
                    point_local[np * 6 + 0] = nl2 * i + (j + 1) * nl_region + k;
                    point_local[np * 6 + 3] = nl2 * i + (j - 1) * nl_region + k;
                }
                // Up
                if (k == 0) {
                    point_local[np * 6 + 2] = nl2 * i + j * nl_region + k + 1;
                    if (geo == 6)
                        point_local[np * 6 + 5] = halo[i * nh2 + 3 * nlhalo + nl_region - j - 1];
                    else {
                        if (j == 0)
                            point_local[np * 6 + 5] = -1;
                        else
                            point_local[np * 6 + 5] =
                                halo[i * nh2 + 3 * nlhalo + nl_region - j - 1];
                    }
                }
                else if (k == nl_region - 1) {
                    point_local[np * 6 + 2] = halo[i * nh2 + nlhalo + j];
                    point_local[np * 6 + 5] = nl2 * i + j * nl_region + k - 1;
                }
                else {
                    point_local[np * 6 + 2] = nl2 * i + j * nl_region + k + 1;
                    point_local[np * 6 + 5] = nl2 * i + j * nl_region + k - 1;
                }
                // Front
                if (j == 0) {
                    if (k == 0) {
                        if (geo == 6)
                            point_local[np * 6 + 4] = halo[i * nh2 + 3 * nlhalo + nl_region];
                        else
                            point_local[np * 6 + 4] = halo[i * nh2 + 3 * nlhalo + nl_region - 1];
                        point_local[np * 6 + 1] = nl2 * i + nl_region + k + 1;
                    }
                    else if (k == nl_region - 1) {
                        point_local[np * 6 + 4] = halo[i * nh2 + nl_region - 2];
                        point_local[np * 6 + 1] = halo[i * nh2 + nl_region + 3];
                    }
                    else {
                        point_local[np * 6 + 4] = halo[i * nh2 + k - 1];
                        point_local[np * 6 + 1] = nl2 * i + (j + 1) * nl_region + k + 1;
                    }
                }
                else if (j == nl_region - 1) {
                    if (k == 0) {
                        point_local[np * 6 + 4] = halo[i * nh2 + 3 * nlhalo + 1];
                        point_local[np * 6 + 1] = halo[i * nh2 + 2 * nlhalo + nl_region - 1];
                    }
                    else if (k == nl_region - 1) {
                        point_local[np * 6 + 4] = nl2 * i + (j - 1) * nl_region + k - 1;
                        point_local[np * 6 + 1] = halo[i * nh2 + 2 * nlhalo];
                    }
                    else {
                        point_local[np * 6 + 4] = nl2 * i + (j - 1) * nl_region + k - 1;
                        point_local[np * 6 + 1] = halo[i * nh2 + 2 * nlhalo - k + nl_region - 1];
                    }
                }
                else {
                    if (k == 0) {
                        point_local[np * 6 + 4] = halo[i * nh2 + 3 * nlhalo - j + nl_region];
                        point_local[np * 6 + 1] = nl2 * i + (j + 1) * nl_region + k + 1;
                    }
                    else if (k == nl_region - 1) {
                        point_local[np * 6 + 4] = nl2 * i + (j - 1) * nl_region + k - 1;
                        point_local[np * 6 + 1] = halo[i * nh2 + nl_region + 2 + j + 1];
                    }
                    else {
                        point_local[np * 6 + 4] = nl2 * i + (j - 1) * nl_region + k - 1;
                        point_local[np * 6 + 1] = nl2 * i + (j + 1) * nl_region + k + 1;
                    }
                }
            }
        }
    }
}


void Icogrid::produce_maps(int *maps, int *halo, int nr, int nl_region) {

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
    int nl2        = nl_region + 2;
    int nl22       = nl2 * nl2;
    int nh2        = 4 * (nl_region + 2);
    int nl_region2 = nl_region * nl_region;

    for (int l = 0; l < nr; l++) {
        for (int i = 0; i < nl_region; i++) {
            for (int j = 0; j < nl_region; j++) {
                maps[l * nl22 + (j + 1) * nl2 + i + 1] = l * nl_region2 + j * nl_region + i;
            }
        }
        //Side 1
        for (int i = 0; i < nl_region + 1; i++) {
            maps[l * nl22 + 0 * (nl_region + 2) + i + 1] = halo[l * nh2 + i];
        }
        //Side 2
        for (int i = 0; i < nl_region + 2; i++) {
            maps[l * nl22 + i * (nl_region + 2) + 0] = halo[l * nh2 + 3 * nl2 + nl_region - i];
        }
        //Side 3
        for (int i = 0; i < nl_region + 1; i++) {
            maps[l * nl22 + (nl_region + 1) * (nl_region + 2) + i + 1] =
                halo[l * nh2 + 2 * nl2 + nl_region - i];
        }
        //Side 4
        for (int i = 0; i < nl_region; i++) {
            maps[l * nl22 + (i + 1) * (nl_region + 2) + nl_region + 1] = halo[l * nh2 + nl2 + i];
        }
    }
}

void Icogrid::spring_dynamics(int *   point_local,
                              int *   pent_ind,
                              int     glevel,
                              double  spring_beta,
                              double *xyz,
                              int     point_num) {

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
    double3 *xyz3 = (double3 *)xyz;

    //  Local variables
    double3 W;
    double3 V;
    double3 vel_val;
    double3 F_v_net;
    double3 D;
    double  l, d, H, max_v;

    // Central point
    double3 *xyzi = new double3[12]();
    double3 *vel  = new double3[point_num]();
    double * Fnet = new double[point_num]();
    double3 *F    = new double3[point_num]();
    double3 *P    = new double3[point_num]();

    // First neighbors.
    double3 *F_Nei = new double3[point_num * 6]();
    double3 *P_Nei = new double3[point_num * 6]();

    // Routine parameters.
    double lba   = 2.0 * M_PI / (10.0 * pow(2.0, glevel - 1));
    double dbar  = spring_beta * lba;
    double drag  = 1;      // Friction coefficient.
    double tstep = 2.0E-2; // time step.
    double cond  = 1.0E-4; // Criteria for convergence.

    int lim = 100000; // Limit of iterations.

    log::printf("\n\n Running spring dynamics.\n\n");

    for (int i = 0; i < point_num; i++) {
        vel[i] = make_double3(0.0, 0.0, 0.0);
    }

    for (int i = 0; i < 12; i++) {
        xyzi[i]                          = xyz3[pent_ind[i]];
        point_local[pent_ind[i] * 6 + 5] = point_local[pent_ind[i] * 6];
    }

    // Solving spring dynamics.
    for (int it = 0; it < lim; it++) {

        for (int i = 0; i < point_num; i++) {
            P[i] = xyz3[i];
            for (int j = 0; j < 6; j++) {
                P_Nei[i * 6 + j] = xyz3[point_local[i * 6 + j]];
            }
        }
        for (int i = 0; i < point_num; i++) {
            for (int j = 0; j < 6; j++) {
                W = cross(P[i], P_Nei[i * 6 + j]);

                V = cross(W, P[i]);

                //norm
                l = length(V);

                d = acos(dot(P[i], P_Nei[i * 6 + j]));

                F_Nei[i * 6 + j] = ((d - dbar) * V) / l;
            }
        }

        for (int i = 0; i < 12; i++) {
            F_Nei[pent_ind[i] * 6 + 5] = make_double3(0.0, 0.0, 0.0);
        }

        for (int i = 0; i < point_num; i++) {
            // Calculate net forces in each point.
            F_v_net = F_Nei[i * 6 + 0] + F_Nei[i * 6 + 1] + F_Nei[i * 6 + 2] + F_Nei[i * 6 + 3]
                      + F_Nei[i * 6 + 4] + F_Nei[i * 6 + 5];

            // Used to check convergence.
            Fnet[i] = length(F_v_net) / lba;

            // Drag move.
            F[i] = F_v_net - drag * vel[i];
        }

        for (int i = 0; i < point_num; i++) {

            // Update points
            D = xyz3[i] + vel[i] * tstep;

            xyz3[i] = normalize(D);
        }

        for (int i = 0; i < point_num; i++) {

            // Update vel
            vel_val = vel[i] + F[i] * tstep;

            H = dot(xyz3[i], vel_val);

            // Remove radial component (if any).
            vel[i] = vel_val - H * xyz3[i];
        }

        // Fix Petagon's position.
        for (int i = 0; i < 12; i++) {
            vel[pent_ind[i]]  = make_double3(0.0, 0.0, 0.0);
            Fnet[pent_ind[i]] = 0.0;
            xyz3[pent_ind[i]] = xyzi[i];
        }

        // Check convergence.
        max_v = -1E-10;
        for (int i = 0; i < point_num; i++)
            if (Fnet[i] > max_v)
                max_v = Fnet[i];
        if (max_v < cond)
            break;
    }

    log::printf(" Done!\n\n");

    delete[] xyzi;
    delete[] vel;
    delete[] F;
    delete[] P;
    delete[] Fnet;
    delete[] F_Nei;
    delete[] P_Nei;
}

void Icogrid::find_qpoints(int *   point_local,
                           double *xyzq,
                           double *xyz,
                           int *   pent_ind,
                           int     point_num) {

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

    double3 *xyz3  = (double3 *)xyz;
    double3 *xyzq3 = (double3 *)xyzq;

    double3 vc1, vc2, vc3;
    double3 sum2;

    int geo;

    for (int i = 0; i < point_num; i++) {
        geo = 6; // Hexagons.
        for (int k = 0; k < 12; k++)
            if (i == pent_ind[k])
                geo = 5; // Pentagons.


        for (int j = 0; j < geo - 1; j++) {
            vc1  = normproj(xyz3[point_local[i * 6 + j]], xyz3[i]);
            vc2  = normproj(xyz3[point_local[i * 6 + j + 1]], xyz3[point_local[i * 6 + j]]);
            vc3  = normproj(xyz3[i], xyz3[point_local[i * 6 + j + 1]]);
            sum2 = sort_add3(vc1, vc2, vc3);

            xyzq3[i * 6 + j] = normalize(sum2);
        }

        vc1 = normproj(xyz3[point_local[i * 6 + geo - 1]], xyz3[i]);
        vc2 = normproj(xyz3[point_local[i * 6 + 0]], xyz3[point_local[i * 6 + geo - 1]]);
        vc3 = normproj(xyz3[i], xyz3[point_local[i * 6 + 0]]);

        sum2 = sort_add3(vc1, vc2, vc3);

        xyzq3[i * 6 + geo - 1] = normalize(sum2);
        if (geo == 5)
            xyzq3[i * 6 + 5] = make_double3(0.0, 0.0, 0.0);
    }
}


void Icogrid::relocate_centres(int *   point_local,
                               double *xyzq,
                               double *xyz,
                               int *   pent_ind,
                               int     point_num) {

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
    double3 *xyzq3 = (double3 *)xyzq;
    double3 *xyz3  = (double3 *)xyz;

    // local vectors
    double3 vgc;

    for (int i = 0; i < point_num; i++) {
        geo = 6; // Hexagons.
        for (int k = 0; k < 12; k++)
            if (i == pent_ind[k])
                geo = 5; // Pentagons.

        vgc = make_double3(0.0, 0.0, 0.0);

        // add all projected vectors
        for (int j = 1; j < geo; j++)
            vgc += normproj(xyzq3[i * 6 + j], xyzq3[i * 6 + j - 1]);

        vgc += normproj(xyzq3[i * 6 + 0], xyzq3[i * 6 + geo - 1]);


        xyz3[i] = normalize(vgc);
    }
}

void Icogrid::set_altitudes_uniform(double *Altitude,
                                    double *Altitudeh,
                                    double  Top_altitude,
                                    int     nv) {

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
    Altitudeh[0]    = 0.0;
    for (int lev = 0; lev < nv; lev++)
        Altitudeh[lev + 1] = Altitudeh[lev] + res_vert;
    for (int lev = 0; lev < nv; lev++) {
        Altitude[lev] = (Altitudeh[lev] + Altitudeh[lev + 1]) / 2.0;
        // printf("Vertical layer, half-layer %d = %f %f\n", lev, Altitude[lev], Altitudeh[lev + 1]);
    }
}

void Icogrid::set_altitudes_refined(double *Altitude,
                                    double *Altitudeh,
                                    double  Top_altitude,
                                    int     nv,
                                    int     n_bl_layers) {

    //
    //  Description:
    //
    //  Sets the layers and interfaces altitudes with refinement in the lower atmosphere.
    //  Useful for turbulent boundary layers of terrestrial planets
    //  The first n_bl_layers are logarithmically spaced (coefficients are ad hoc)
    //  Above that, the remaining layers are uniformly spaced
    //  Algorithm written by Pierre Auclair-Desrotours
    //
    //  Input:  - nv - Number of vertical layers.
    //          - top_altitude - Altitude of the top model domain.
    //
    //  Output: - Altitude  - Layers altitudes.
    //          - Altitudeh - Interfaces altitudes.
    //
    int    ntrans  = n_bl_layers + 2;
    double xup     = nv * 1.0 / ntrans;
    double zup     = Top_altitude;
    double a       = -4.4161; //These numbers are completely ad hoc.
    double b       = 12.925;  //May need to be adjusted depending on performance, &c
    double c       = -10.614;
    double err_tar = 1.0E-8;

    // stitch the two domains together (log and linear regions)
    // solve for intersection via secant method
    int    nitmax = 100, j = 0;
    double x1 = 0.5;
    double x2 = 1.0;
    double f1 = exp(a * pow(x1, 2) + b * x1 + c) * (1 + (xup - x1) * (2 * a * x1 + b)) - 1;
    double f2 = exp(a * pow(x2, 2) + b * x2 + c) * (1 + (xup - x2) * (2 * a * x2 + b)) - 1;
    double xnew, fnew, err = 1.0, xbl, xh, d, k;
    int    levbl, levh;

    while (j < nitmax && err > err_tar) {
        xnew = (x1 * f2 - x2 * f1) / (f2 - f1);
        fnew = exp(a * pow(xnew, 2) + b * xnew + c) * (1 + (xup - xnew) * (2 * a * xnew + b)) - 1;
        err  = fabs(xnew - x2);
        x1   = x2;
        x2   = xnew;
        f1   = f2;
        f2   = fnew;
        j++;
    }

    xbl   = x2;
    levbl = floor(xbl * ntrans);
    d     = (2 * a * xbl + b) * exp(a * pow(xbl, 2) + b * xbl + c);
    k     = 1 - d * xup;

    Altitudeh[0] = 0.0;
    for (int lev = 0; lev < nv; lev++) {
        levh = lev + 1;
        xh   = levh * 1.0 / ntrans;
        if (levh > levbl) {
            Altitudeh[levh] = zup * (d * xh + k);
        }
        else {
            Altitudeh[levh] = zup * exp(a * pow(xh, 2) + b * xh + c);
        }
    }
    for (int lev = 0; lev < nv; lev++) {
        Altitude[lev] = (Altitudeh[lev] + Altitudeh[lev + 1]) / 2.0;
        // printf("Vertical layer, half-layer %d = %f %f\n", lev, Altitude[lev], Altitudeh[lev + 1]);
    }
}

void Icogrid::set_altitudes_softplus(double *Altitude,
                                     double *Altitudeh,
                                     double  Top_altitude,
                                     double  lowest_layer_thickness,
                                     double  transition_altitude,
                                     int     nv) {


    //
    //  Description:
    //
    //  Sets the layers and interfaces altitudes using softplus function
    //  (exponential below transition_altitude, linear above)
    //
    //  Input:  - nv - Number of vertical layers.
    //          - top_altitude - Altitude of the top model domain.
    //          - lowest_layer_thickness - Thickness of lowest layer
    //          - transition_altitude - Transition b/w exponential and linear
    //
    //  Output: - Altitude  - Layers altitudes.
    //          - Altitudeh - Interfaces altitudes.
    //

    //parameters controlling shape of softplus function
    double alpha, k, xbl, x1;
    x1 = 1.0 / nv; // fractional index of first layer top
    // Calculate sharpness and amplitude to match top, bottom, transition altitude
    alpha = transition_altitude / log(2);
    k = log((exp(lowest_layer_thickness / alpha) - 1) / (exp(Top_altitude / alpha) - 1)) / (x1 - 1);
    xbl = -log(exp(lowest_layer_thickness / alpha) - 1) / k + x1; //centering the function

    double x     = 0;
    Altitudeh[0] = 0.0;
    for (int lev = 0; lev < nv; lev++) { //interfaces
        x += x1;                         //next fractional index of layer
        Altitudeh[lev + 1] = alpha * log(1 + exp(k * (x - xbl)));
    }
    for (int lev = 0; lev < nv; lev++) { //centers of layers
        Altitude[lev] = (Altitudeh[lev] + Altitudeh[lev + 1]) / 2.0;
    }
    printf("stop here");
}


void Icogrid::cart2sphe(double *lonlat, double *xyz, int point_num) {

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
    for (int i = 0; i < point_num; i++) {
        lonlat[i * 2 + 0] = atan2(xyz[i * 3 + 1], xyz[i * 3 + 0]);
        lonlat[i * 2 + 1] =
            atan2(xyz[i * 3 + 2], sqrt(pow(xyz[i * 3 + 0], 2) + pow(xyz[i * 3 + 1], 2)));
    }
}

void Icogrid::correct_xyz_points(double  A,
                                 double *xyzq,
                                 double *xyz,
                                 int *   pent_ind,
                                 int     point_num) {

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

    double3 *xyzq3 = (double3 *)xyzq;
    double3 *xyz3  = (double3 *)xyz;

    for (int i = 0; i < point_num; i++) {
        geo = 6; // Hexagons.
        for (int k = 0; k < 12; k++)
            if (i == pent_ind[k])
                geo = 5; // Pentagons.
        if (geo == 5) {
            for (int j = 0; j < 5; j++) { // Pentagons.
                xyzq3[i * 6 + j] = xyzq3[i * 6 + j] * A;
            }
            xyz3[i] = xyz3[i] * A;
        }
        else {
            for (int j = 0; j < 6; j++) { // Hexagons.
                xyzq3[i * 6 + j] = xyzq3[i * 6 + j + 0] * A;
            }
            xyz3[i] = xyz3[i] * A;
        }
    }
}

inline double comp_ang(const double3 &a, const double3 &b, const double3 &c) {
    double3 v0 = a;
    double3 v1 = b - a;
    double3 v2 = c - a;

    double t = 1.0 / dot(v0, v0);

    double fac1 = dot(v0, v1) * t;
    double fac2 = dot(v0, v2) * t;

    double3 w1 = v1 - fac1 * v0;

    double3 w2 = v2 - fac2 * v0;

    double ang = dot(w1, w2) / (length(w1) * length(w2));

    if (ang > 1)
        ang = 1;
    if (ang < -1)
        ang = -1;

    return acos(ang);
}

void Icogrid::control_areas(double *areasT,
                            double *areasTr,
                            double *areas,
                            int *   point_local,
                            double *xyzq,
                            double *xyz,
                            int *   pent_ind,
                            int     point_num,
                            double  A) {

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

    // input variables
    double3 *xyzq3 = (double3 *)xyzq;
    double3 *xyz3  = (double3 *)xyz;

    // Local variables.
    int geo;

    double radius;

    // Local arrays
    double3 a, b, c;

    double ang[3];
    double areav;
    double area_tot_exact, area_tot_T, area_tot_Tr, area_tot_subTr;
    area_tot_exact = 4 * M_PI * pow(A, 2.0); //exact surface area
    area_tot_T     = 0.0;                    // surface area summed over control volumes
    area_tot_Tr    = 0.0;                    // surface area summed over triangles
    area_tot_subTr = 0.0;                    // surface area summed over sub-triangles

    for (int i = 0; i < point_num * 6 * 3; i++)
        areas[i] = 0.0;

    for (int i = 0; i < point_num; i++) {
        //
        // Areas - Alpha, Beta, Gamma
        //
        geo = 6; // Hexagons.
        for (int k = 0; k < 12; k++)
            if (i == pent_ind[k])
                geo = 5; // Pentagons.
        if (geo == 5) {
            for (int j = 0; j < 5; j++) { //Pentagons
                for (int k = 0; k < 3; k++) {
                    if (j < 4) {
                        if (k == 0) {
                            a = xyz3[point_local[i * 6 + j + 1]];
                            b = xyz3[point_local[i * 6 + j]];
                            c = xyzq3[i * 6 + j];
                        }
                        else if (k == 1) {
                            a = xyz3[point_local[i * 6 + j + 1]];
                            b = xyz3[i];
                            c = xyzq3[i * 6 + j];
                        }
                        else {
                            a = xyz3[i];
                            b = xyz3[point_local[i * 6 + j]];
                            c = xyzq3[i * 6 + j];
                        }
                    }
                    else {
                        if (k == 0) {
                            a = xyz3[point_local[i * 6 + 0]];
                            b = xyz3[point_local[i * 6 + 4]];
                            c = xyzq3[i * 6 + 4];
                        }
                        else if (k == 1) {
                            a = xyz3[point_local[i * 6 + 0]];
                            b = xyz3[i];
                            c = xyzq3[i * 6 + 4];
                        }
                        else {
                            a = xyz3[i];
                            b = xyz3[point_local[i * 6 + 4]];
                            c = xyzq3[i * 6 + 4];
                        }
                    }

                    ang[0] = comp_ang(a, b, c);
                    ang[1] = comp_ang(b, a, c);
                    ang[2] = comp_ang(c, a, b);

                    radius = length(a);
                    areas[i * 6 * 3 + j * 3 + k] =
                        (ang[0] + ang[1] + ang[2] - M_PI) * pow(radius, 2);
                }
                areasTr[i * 6 + j] = areas[(i * 6 + j) * 3 + 0] + areas[(i * 6 + j) * 3 + 1]
                                     + areas[(i * 6 + j) * 3 + 2];
            }
            areasTr[i * 6 + 5] = 0.0;
        }
        else {
            for (int j = 0; j < 6; j++) { //Hexagons
                for (int k = 0; k < 3; k++) {
                    if (j < 5) {
                        if (k == 0) {
                            a = xyz3[point_local[i * 6 + j + 1]];
                            b = xyz3[point_local[i * 6 + j]];
                            c = xyzq3[i * 6 + j];
                        }
                        else if (k == 1) {
                            a = xyz3[point_local[i * 6 + j + 1]];
                            b = xyz3[i];
                            c = xyzq3[i * 6 + j];
                        }
                        else {
                            a = xyz3[i];
                            b = xyz3[point_local[i * 6 + j]];
                            c = xyzq3[i * 6 + j];
                        }
                    }
                    else {
                        if (k == 0) {
                            a = xyz3[point_local[i * 6 + 0]];
                            b = xyz3[point_local[i * 6 + 5]];
                            c = xyzq3[i * 6 + 5];
                        }
                        else if (k == 1) {
                            a = xyz3[point_local[i * 6 + 0]];
                            b = xyz3[i];
                            c = xyzq3[i * 6 + 5];
                        }
                        else {
                            a = xyz3[i];
                            b = xyz3[point_local[i * 6 + 5]];
                            c = xyzq3[i * 6 + 5];
                        }
                    }

                    ang[0] = comp_ang(a, b, c);
                    ang[1] = comp_ang(b, a, c);
                    ang[2] = comp_ang(c, a, b);

                    radius = length(a);
                    areas[i * 6 * 3 + j * 3 + k] =
                        (ang[0] + ang[1] + ang[2] - M_PI) * pow(radius, 2);
                }
                areasTr[i * 6 + j] = areas[(i * 6 + j) * 3 + 0] + areas[(i * 6 + j) * 3 + 1]
                                     + areas[(i * 6 + j) * 3 + 2];
            }
        }
        if (i < point_num - 2) {
            //sum over all areas (polar triangles will be counted with neighboring vertices)
            for (int j = 0; j < 2; j++) {
                //sum only two per vertex to prevent double/triple counting
                area_tot_Tr += areasTr[i * 6 + j];
                area_tot_subTr += (areas[i * 6 * 3 + j * 3 + 0] + areas[i * 6 * 3 + j * 3 + 1]
                                   + areas[i * 6 * 3 + j * 3 + 2]);
            }
        }
    }

    for (int i = 0; i < point_num; i++) {
        //
        //      Control Volumes
        //
        areasT[i] = 0.0;

        geo = 6; // Hexagons.
        for (int k = 0; k < 12; k++)
            if (i == pent_ind[k])
                geo = 5; // Pentagons.
        for (int j = 0; j < 6; j++) {
            if (geo == 5 && j == 4) {
                a = xyz3[i];
                b = xyzq3[i * 6 + 4];
                c = xyzq3[i * 6 + 0];
            }
            else if (geo == 5 && j == 5) {
                break;
            }
            else if (geo == 6 && j == 5) {
                a = xyz3[i];
                b = xyzq3[i * 6 + 5];
                c = xyzq3[i * 6 + 0];
            }
            else {
                a = xyz3[i];
                b = xyzq3[i * 6 + j];
                c = xyzq3[i * 6 + (j + 1)];
            }

            ang[0] = comp_ang(a, b, c);
            ang[1] = comp_ang(b, a, c);
            ang[2] = comp_ang(c, a, b);

            radius = length(a);
            areav  = (ang[0] + ang[1] + ang[2] - M_PI) * pow(radius, 2);

            areasT[i] += areav;
        }
        area_tot_T += areasT[i];
    }

    //normalize areas to get exact area of sphere
    for (int i = 0; i < point_num; i++) {
        areasT[i] *= area_tot_exact / area_tot_T;
        for (int j = 0; j < 6; j++) {
            areasTr[i * 6 + j] *= area_tot_exact / area_tot_Tr;
            for (int k = 0; k < 3; k++) {
                areas[i * 6 * 3 + j * 3 + k] *= area_tot_exact / area_tot_subTr;
            }
        }
    }
}


void Icogrid::control_vec(double *nvec,
                          double *nvecoa,
                          double *nvecti,
                          double *nvecte,
                          double *areasT,
                          double *xyz,
                          double *xyzq,
                          int *   point_local,
                          int *   pent_ind,
                          int     point_num,
                          double *mvec) {

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
    double3 *nvec3   = (double3 *)nvec;
    double3 *nvecoa3 = (double3 *)nvecoa;
    double3 *nvecti3 = (double3 *)nvecti;
    double3 *nvecte3 = (double3 *)nvecte;

    double3 *mvec3 = (double3 *)mvec;

    double3 *xyz3  = (double3 *)xyz;
    double3 *xyzq3 = (double3 *)xyzq;

    // Local variables.
    double vec_l, l;
    double fac_nv, fac_mv;

    int geo;

    // Local arrays.
    double3 v1, v2;
    double3 nv, mv;

    for (int i = 0; i < point_num; i++) {

        geo = 6; // Hexagons.
        for (int k = 0; k < 12; k++)
            if (i == pent_ind[k])
                geo = 5; // Pentagons.

        if (geo == 5) {

            for (int j = 0; j < 5; j++) {

                if (j < 4) {
                    v2 = xyzq3[i * 6 + (j + 1)];
                    v1 = xyzq3[i * 6 + j];
                }
                else {
                    v2 = xyzq3[i * 6 + 0];
                    v1 = xyzq3[i * 6 + 4];
                }

                vec_l = length(v1);

                l = acos(dot(v1, v2) / (pow(vec_l, 2.0))) * vec_l;

                nv = cross(v1, v2);


                fac_nv = l / length(nv);

                nvec3[i * 6 + j] = nv * fac_nv;

                nvecoa3[i * 6 + j] = nvec3[i * 6 + j] / areasT[i];

                mv               = v1 - v2; //unnormalized vector parallel to cv edges
                fac_mv           = l / length(mv);
                mvec3[i * 6 + j] = mv * fac_mv; //now normalized and l is incorporated
            }

            nvec3[i * 6 + 5] = make_double3(0.0, 0.0, 0.0);

            mvec3[i * 6 + 5] = make_double3(0.0, 0.0, 0.0);

            nvecoa3[i * 6 + 5] = make_double3(0.0, 0.0, 0.0);
        }
        else {

            for (int j = 0; j < 6; j++) {

                if (j < 5) {
                    v2 = xyzq3[i * 6 + (j + 1)];
                    v1 = xyzq3[i * 6 + j];
                }
                else {
                    v2 = xyzq3[i * 6 + 0];
                    v1 = xyzq3[i * 6 + 5];
                }
                vec_l = length(v1);

                l = acos(dot(v1, v2) / (pow(vec_l, 2.0))) * vec_l;

                nv = cross(v1, v2);


                fac_nv = l / length(nv);

                nvec3[i * 6 + j] = nv * fac_nv;

                nvecoa3[i * 6 + j] = nvec3[i * 6 + j] / areasT[i];

                mv               = v1 - v2; //unnormalized vector parallel to cv edges
                fac_mv           = l / length(mv);
                mvec3[i * 6 + j] = mv * fac_mv; //now normalized and l is incorporated
            }
        }

        // nvecti & nvecte! Triangles.
        if (geo == 5) {

            for (int j = 0; j < 5; j++) {
                v1 = xyz3[point_local[i * 6 + j]];
                v2 = xyz3[i];

                vec_l = length(v1);

                l = acos(dot(v1, v2) / (pow(vec_l, 2.0))) * vec_l;

                nv = cross(v1, v2);

                fac_nv = l / length(nv);

                nvecti3[i * 6 + j] = nv * fac_nv;
            }
            nvecti3[i * 6 + 5] = make_double3(0.0, 0.0, 0.0);
        }
        else {
            for (int j = 0; j < 6; j++) {
                v1 = xyz3[point_local[i * 6 + j]];
                v2 = xyz3[i];

                vec_l = length(v1);

                l = acos(dot(v1, v2) / (pow(vec_l, 2.0))) * vec_l;

                nv = cross(v1, v2);

                fac_nv = l / length(nv);

                nvecti3[i * 6 + j] = nv * fac_nv;
            }
        }

        if (geo == 5) {

            for (int j = 0; j < 5; j++) {

                if (j < 4) {
                    v1 = xyz3[point_local[i * 6 + j]];
                    v2 = xyz3[point_local[i * 6 + j + 1]];
                }
                else {
                    v1 = xyz3[point_local[i * 6 + 4]];
                    v2 = xyz3[point_local[i * 6 + 0]];
                }

                vec_l = length(v1);

                l = acos(dot(v1, v2) / (pow(vec_l, 2.0))) * vec_l;

                nv = cross(v1, v2);

                fac_nv = l / length(nv);

                nvecte3[i * 6 + j] = nv * fac_nv;
            }
            nvecti3[i * 6 + 5] = make_double3(0.0, 0.0, 0.0);
        }
        else {
            for (int j = 0; j < 6; j++) {

                if (j < 5) {
                    v1 = xyz3[point_local[i * 6 + j]];
                    v2 = xyz3[point_local[i * 6 + j + 1]];
                }
                else {
                    v1 = xyz3[point_local[i * 6 + 5]];
                    v2 = xyz3[point_local[i * 6 + 0]];
                }

                vec_l = length(v1);

                l = acos(dot(v1, v2) / (pow(vec_l, 2.0))) * vec_l;

                nv = cross(v1, v2);

                fac_nv = l / length(nv);

                nvecte3[i * 6 + j] = nv * fac_nv;
            }
        }
    }
}


void Icogrid::compute_func(double *func_r, double *xyz, int point_num) {

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

    double3 *func_r3 = (double3 *)func_r;
    double3 *xyz3    = (double3 *)xyz;


    for (int i = 0; i < point_num; i++)
        func_r3[i] = normalize(xyz3[i]);
}

void Icogrid::div_operator(double *areasT,
                           double *areas,
                           double *div,
                           double *nvec,
                           int *   pent_ind,
                           int     point_num) {

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
    int    geo;
    double area1, area2, area3, area4, area5, area6;

    double *areasq2;

    areasq2 = new double[6 * 3 * point_num]();

    for (int i = 0; i < point_num; i++) {

        geo = 6; // Hexagons.
        for (int k = 0; k < 12; k++)
            if (i == pent_ind[k])
                geo = 5; // Pentagons.

        area1 = areas[i * 6 * 3 + 0 * 3 + 0] + areas[i * 6 * 3 + 0 * 3 + 1]
                + areas[i * 6 * 3 + 0 * 3 + 2];
        area2 = areas[i * 6 * 3 + 1 * 3 + 0] + areas[i * 6 * 3 + 1 * 3 + 1]
                + areas[i * 6 * 3 + 1 * 3 + 2];
        area3 = areas[i * 6 * 3 + 2 * 3 + 0] + areas[i * 6 * 3 + 2 * 3 + 1]
                + areas[i * 6 * 3 + 2 * 3 + 2];
        area4 = areas[i * 6 * 3 + 3 * 3 + 0] + areas[i * 6 * 3 + 3 * 3 + 1]
                + areas[i * 6 * 3 + 3 * 3 + 2];
        area5 = areas[i * 6 * 3 + 4 * 3 + 0] + areas[i * 6 * 3 + 4 * 3 + 1]
                + areas[i * 6 * 3 + 4 * 3 + 2];

        if (geo == 5)
            area6 = 0;
        else
            area6 = areas[i * 6 * 3 + 5 * 3 + 0] + areas[i * 6 * 3 + 5 * 3 + 1]
                    + areas[i * 6 * 3 + 5 * 3 + 2];

        for (int k = 0; k < 3; k++) {
            areasq2[i * 6 * 3 + 0 * 3 + k] = areas[i * 6 * 3 + 0 * 3 + k] / area1;
            areasq2[i * 6 * 3 + 1 * 3 + k] = areas[i * 6 * 3 + 1 * 3 + k] / area2;
            areasq2[i * 6 * 3 + 2 * 3 + k] = areas[i * 6 * 3 + 2 * 3 + k] / area3;
            areasq2[i * 6 * 3 + 3 * 3 + k] = areas[i * 6 * 3 + 3 * 3 + k] / area4;
            areasq2[i * 6 * 3 + 4 * 3 + k] = areas[i * 6 * 3 + 4 * 3 + k] / area5;
            if (geo == 5)
                areasq2[i * 6 * 3 + 5 * 3 + k] = 0;
            else
                areasq2[i * 6 * 3 + 5 * 3 + k] = areas[i * 6 * 3 + 5 * 3 + k] / area6;
        }
    }

    for (int i = 0; i < point_num; i++) {

        geo = 6; // Hexagons.
        for (int k = 0; k < 12; k++)
            if (i == pent_ind[k])
                geo = 5; // Pentagons.

        for (int k = 0; k < 3; k++) {
            if (geo == 5) { //Pentagons.
                div[i * 7 * 3 + 0 * 3 + k] =
                    nvec[i * 6 * 3 + 0 * 3 + k]
                        * (areasq2[i * 6 * 3 + 0 * 3 + 0] + areasq2[i * 6 * 3 + 1 * 3 + 0])
                    + nvec[i * 6 * 3 + 1 * 3 + k]
                          * (areasq2[i * 6 * 3 + 1 * 3 + 0] + areasq2[i * 6 * 3 + 2 * 3 + 0])
                    + nvec[i * 6 * 3 + 2 * 3 + k]
                          * (areasq2[i * 6 * 3 + 2 * 3 + 0] + areasq2[i * 6 * 3 + 3 * 3 + 0])
                    + nvec[i * 6 * 3 + 3 * 3 + k]
                          * (areasq2[i * 6 * 3 + 3 * 3 + 0] + areasq2[i * 6 * 3 + 4 * 3 + 0])
                    + nvec[i * 6 * 3 + 4 * 3 + k]
                          * (areasq2[i * 6 * 3 + 4 * 3 + 0] + areasq2[i * 6 * 3 + 0 * 3 + 0]);

                div[i * 7 * 3 + 1 * 3 + k] =
                    nvec[i * 6 * 3 + 4 * 3 + k]
                        * (areasq2[i * 6 * 3 + 0 * 3 + 1] + areasq2[i * 6 * 3 + 4 * 3 + 2])
                    + nvec[i * 6 * 3 + 3 * 3 + k] * areasq2[i * 6 * 3 + 4 * 3 + 2]
                    + nvec[i * 6 * 3 + 0 * 3 + k] * areasq2[i * 6 * 3 + 0 * 3 + 1];

                div[i * 7 * 3 + 2 * 3 + k] =
                    nvec[i * 6 * 3 + 0 * 3 + k]
                        * (areasq2[i * 6 * 3 + 1 * 3 + 1] + areasq2[i * 6 * 3 + 0 * 3 + 2])
                    + nvec[i * 6 * 3 + 4 * 3 + k] * areasq2[i * 6 * 3 + 0 * 3 + 2]
                    + nvec[i * 6 * 3 + 1 * 3 + k] * areasq2[i * 6 * 3 + 1 * 3 + 1];

                div[i * 7 * 3 + 3 * 3 + k] =
                    nvec[i * 6 * 3 + 1 * 3 + k]
                        * (areasq2[i * 6 * 3 + 2 * 3 + 1] + areasq2[i * 6 * 3 + 1 * 3 + 2])
                    + nvec[i * 6 * 3 + 0 * 3 + k] * areasq2[i * 6 * 3 + 1 * 3 + 2]
                    + nvec[i * 6 * 3 + 2 * 3 + k] * areasq2[i * 6 * 3 + 2 * 3 + 1];

                div[i * 7 * 3 + 4 * 3 + k] =
                    nvec[i * 6 * 3 + 2 * 3 + k]
                        * (areasq2[i * 6 * 3 + 3 * 3 + 1] + areasq2[i * 6 * 3 + 2 * 3 + 2])
                    + nvec[i * 6 * 3 + 1 * 3 + k] * areasq2[i * 6 * 3 + 2 * 3 + 2]
                    + nvec[i * 6 * 3 + 3 * 3 + k] * areasq2[i * 6 * 3 + 3 * 3 + 1];

                div[i * 7 * 3 + 5 * 3 + k] =
                    nvec[i * 6 * 3 + 3 * 3 + k]
                        * (areasq2[i * 6 * 3 + 4 * 3 + 1] + areasq2[i * 6 * 3 + 3 * 3 + 2])
                    + nvec[i * 6 * 3 + 2 * 3 + k] * areasq2[i * 6 * 3 + 3 * 3 + 2]
                    + nvec[i * 6 * 3 + 4 * 3 + k] * areasq2[i * 6 * 3 + 4 * 3 + 1];

                div[i * 7 * 3 + 6 * 3 + k] = 0.0;
            }
            else {
                div[i * 7 * 3 + 0 * 3 + k] =
                    nvec[i * 6 * 3 + 0 * 3 + k]
                        * (areasq2[i * 6 * 3 + 0 * 3 + 0] + areasq2[i * 6 * 3 + 1 * 3 + 0])
                    + nvec[i * 6 * 3 + 1 * 3 + k]
                          * (areasq2[i * 6 * 3 + 1 * 3 + 0] + areasq2[i * 6 * 3 + 2 * 3 + 0])
                    + nvec[i * 6 * 3 + 2 * 3 + k]
                          * (areasq2[i * 6 * 3 + 2 * 3 + 0] + areasq2[i * 6 * 3 + 3 * 3 + 0])
                    + nvec[i * 6 * 3 + 3 * 3 + k]
                          * (areasq2[i * 6 * 3 + 3 * 3 + 0] + areasq2[i * 6 * 3 + 4 * 3 + 0])
                    + nvec[i * 6 * 3 + 4 * 3 + k]
                          * (areasq2[i * 6 * 3 + 4 * 3 + 0] + areasq2[i * 6 * 3 + 5 * 3 + 0])
                    + nvec[i * 6 * 3 + 5 * 3 + k]
                          * (areasq2[i * 6 * 3 + 5 * 3 + 0] + areasq2[i * 6 * 3 + 0 * 3 + 0]);

                div[i * 7 * 3 + 1 * 3 + k] =
                    nvec[i * 6 * 3 + 5 * 3 + k]
                        * (areasq2[i * 6 * 3 + 0 * 3 + 1] + areasq2[i * 6 * 3 + 5 * 3 + 2])
                    + nvec[i * 6 * 3 + 4 * 3 + k] * areasq2[i * 6 * 3 + 5 * 3 + 2]
                    + nvec[i * 6 * 3 + 0 * 3 + k] * areasq2[i * 6 * 3 + 0 * 3 + 1];

                div[i * 7 * 3 + 2 * 3 + k] =
                    nvec[i * 6 * 3 + 0 * 3 + k]
                        * (areasq2[i * 6 * 3 + 1 * 3 + 1] + areasq2[i * 6 * 3 + 0 * 3 + 2])
                    + nvec[i * 6 * 3 + 5 * 3 + k] * areasq2[i * 6 * 3 + 0 * 3 + 2]
                    + nvec[i * 6 * 3 + 1 * 3 + k] * areasq2[i * 6 * 3 + 1 * 3 + 1];

                div[i * 7 * 3 + 3 * 3 + k] =
                    nvec[i * 6 * 3 + 1 * 3 + k]
                        * (areasq2[i * 6 * 3 + 2 * 3 + 1] + areasq2[i * 6 * 3 + 1 * 3 + 2])
                    + nvec[i * 6 * 3 + 0 * 3 + k] * areasq2[i * 6 * 3 + 1 * 3 + 2]
                    + nvec[i * 6 * 3 + 2 * 3 + k] * areasq2[i * 6 * 3 + 2 * 3 + 1];

                div[i * 7 * 3 + 4 * 3 + k] =
                    nvec[i * 6 * 3 + 2 * 3 + k]
                        * (areasq2[i * 6 * 3 + 3 * 3 + 1] + areasq2[i * 6 * 3 + 2 * 3 + 2])
                    + nvec[i * 6 * 3 + 1 * 3 + k] * areasq2[i * 6 * 3 + 2 * 3 + 2]
                    + nvec[i * 6 * 3 + 3 * 3 + k] * areasq2[i * 6 * 3 + 3 * 3 + 1];

                div[i * 7 * 3 + 5 * 3 + k] =
                    nvec[i * 6 * 3 + 3 * 3 + k]
                        * (areasq2[i * 6 * 3 + 4 * 3 + 1] + areasq2[i * 6 * 3 + 3 * 3 + 2])
                    + nvec[i * 6 * 3 + 2 * 3 + k] * areasq2[i * 6 * 3 + 3 * 3 + 2]
                    + nvec[i * 6 * 3 + 4 * 3 + k] * areasq2[i * 6 * 3 + 4 * 3 + 1];

                div[i * 7 * 3 + 6 * 3 + k] =
                    nvec[i * 6 * 3 + 4 * 3 + k]
                        * (areasq2[i * 6 * 3 + 5 * 3 + 1] + areasq2[i * 6 * 3 + 4 * 3 + 2])
                    + nvec[i * 6 * 3 + 3 * 3 + k] * areasq2[i * 6 * 3 + 4 * 3 + 2]
                    + nvec[i * 6 * 3 + 5 * 3 + k] * areasq2[i * 6 * 3 + 5 * 3 + 1];
            }
        }
        for (int j = 0; j < 7; j++)
            for (int k = 0; k < 3; k++)
                div[i * 7 * 3 + j * 3 + k] = div[i * 7 * 3 + j * 3 + k] / (2.0 * areasT[i]);
    }
    delete[] areasq2;
}

void Icogrid::curlz_operator(double *areasT,
                             double *areas,
                             double *curlz,
                             double *mvec,
                             int *   pent_ind,
                             int     point_num) {

    //
    //  Description:
    //
    //  Computes the vertical curl operator from equation 10 of Tomita+ 2001
    //  (There is probably no clean way of calculating the horizontal components,
    //   which should be small in most situations)
    //
    //  This is currently not used by the code, it is only calculated and output
    //  for convenience in post-processing (calculating vorticity, e.g.)
    //
    //  Input: - areas     - Sub-areas (alpha, beta and gamma).
    //         - areasT    - Control volume area.
    //         - mvec      - Vectors parallel to the edges of the control volume midway between vertices
    //         - pent_ind  - Pentagon' indexes.
    //         - point_num - Number of vertices.
    //
    //  Output: - curlz - vertical curl operator.
    //
    // Local variables:
    int    geo;
    double area1, area2, area3, area4, area5, area6;

    double *areasq2;

    areasq2 = new double[6 * 3 * point_num]();

    for (int i = 0; i < point_num; i++) {

        geo = 6; // Hexagons.
        for (int k = 0; k < 12; k++)
            if (i == pent_ind[k])
                geo = 5; // Pentagons.

        area1 = areas[i * 6 * 3 + 0 * 3 + 0] + areas[i * 6 * 3 + 0 * 3 + 1]
                + areas[i * 6 * 3 + 0 * 3 + 2];
        area2 = areas[i * 6 * 3 + 1 * 3 + 0] + areas[i * 6 * 3 + 1 * 3 + 1]
                + areas[i * 6 * 3 + 1 * 3 + 2];
        area3 = areas[i * 6 * 3 + 2 * 3 + 0] + areas[i * 6 * 3 + 2 * 3 + 1]
                + areas[i * 6 * 3 + 2 * 3 + 2];
        area4 = areas[i * 6 * 3 + 3 * 3 + 0] + areas[i * 6 * 3 + 3 * 3 + 1]
                + areas[i * 6 * 3 + 3 * 3 + 2];
        area5 = areas[i * 6 * 3 + 4 * 3 + 0] + areas[i * 6 * 3 + 4 * 3 + 1]
                + areas[i * 6 * 3 + 4 * 3 + 2];

        if (geo == 5)
            area6 = 0;
        else
            area6 = areas[i * 6 * 3 + 5 * 3 + 0] + areas[i * 6 * 3 + 5 * 3 + 1]
                    + areas[i * 6 * 3 + 5 * 3 + 2];

        for (int k = 0; k < 3; k++) {
            areasq2[i * 6 * 3 + 0 * 3 + k] = areas[i * 6 * 3 + 0 * 3 + k] / area1;
            areasq2[i * 6 * 3 + 1 * 3 + k] = areas[i * 6 * 3 + 1 * 3 + k] / area2;
            areasq2[i * 6 * 3 + 2 * 3 + k] = areas[i * 6 * 3 + 2 * 3 + k] / area3;
            areasq2[i * 6 * 3 + 3 * 3 + k] = areas[i * 6 * 3 + 3 * 3 + k] / area4;
            areasq2[i * 6 * 3 + 4 * 3 + k] = areas[i * 6 * 3 + 4 * 3 + k] / area5;
            if (geo == 5)
                areasq2[i * 6 * 3 + 5 * 3 + k] = 0;
            else
                areasq2[i * 6 * 3 + 5 * 3 + k] = areas[i * 6 * 3 + 5 * 3 + k] / area6;
        }
    }

    for (int i = 0; i < point_num; i++) {

        geo = 6; // Hexagons.
        for (int k = 0; k < 12; k++)
            if (i == pent_ind[k])
                geo = 5; // Pentagons.

        for (int k = 0; k < 3; k++) {
            if (geo == 5) { //Pentagons.
                curlz[i * 7 * 3 + 0 * 3 + k] =
                    mvec[i * 6 * 3 + 0 * 3 + k]
                        * (areasq2[i * 6 * 3 + 0 * 3 + 0] + areasq2[i * 6 * 3 + 1 * 3 + 0])
                    + mvec[i * 6 * 3 + 1 * 3 + k]
                          * (areasq2[i * 6 * 3 + 1 * 3 + 0] + areasq2[i * 6 * 3 + 2 * 3 + 0])
                    + mvec[i * 6 * 3 + 2 * 3 + k]
                          * (areasq2[i * 6 * 3 + 2 * 3 + 0] + areasq2[i * 6 * 3 + 3 * 3 + 0])
                    + mvec[i * 6 * 3 + 3 * 3 + k]
                          * (areasq2[i * 6 * 3 + 3 * 3 + 0] + areasq2[i * 6 * 3 + 4 * 3 + 0])
                    + mvec[i * 6 * 3 + 4 * 3 + k]
                          * (areasq2[i * 6 * 3 + 4 * 3 + 0] + areasq2[i * 6 * 3 + 0 * 3 + 0]);

                curlz[i * 7 * 3 + 1 * 3 + k] =
                    mvec[i * 6 * 3 + 4 * 3 + k]
                        * (areasq2[i * 6 * 3 + 0 * 3 + 1] + areasq2[i * 6 * 3 + 4 * 3 + 2])
                    + mvec[i * 6 * 3 + 3 * 3 + k] * areasq2[i * 6 * 3 + 4 * 3 + 2]
                    + mvec[i * 6 * 3 + 0 * 3 + k] * areasq2[i * 6 * 3 + 0 * 3 + 1];

                curlz[i * 7 * 3 + 2 * 3 + k] =
                    mvec[i * 6 * 3 + 0 * 3 + k]
                        * (areasq2[i * 6 * 3 + 1 * 3 + 1] + areasq2[i * 6 * 3 + 0 * 3 + 2])
                    + mvec[i * 6 * 3 + 4 * 3 + k] * areasq2[i * 6 * 3 + 0 * 3 + 2]
                    + mvec[i * 6 * 3 + 1 * 3 + k] * areasq2[i * 6 * 3 + 1 * 3 + 1];

                curlz[i * 7 * 3 + 3 * 3 + k] =
                    mvec[i * 6 * 3 + 1 * 3 + k]
                        * (areasq2[i * 6 * 3 + 2 * 3 + 1] + areasq2[i * 6 * 3 + 1 * 3 + 2])
                    + mvec[i * 6 * 3 + 0 * 3 + k] * areasq2[i * 6 * 3 + 1 * 3 + 2]
                    + mvec[i * 6 * 3 + 2 * 3 + k] * areasq2[i * 6 * 3 + 2 * 3 + 1];

                curlz[i * 7 * 3 + 4 * 3 + k] =
                    mvec[i * 6 * 3 + 2 * 3 + k]
                        * (areasq2[i * 6 * 3 + 3 * 3 + 1] + areasq2[i * 6 * 3 + 2 * 3 + 2])
                    + mvec[i * 6 * 3 + 1 * 3 + k] * areasq2[i * 6 * 3 + 2 * 3 + 2]
                    + mvec[i * 6 * 3 + 3 * 3 + k] * areasq2[i * 6 * 3 + 3 * 3 + 1];

                curlz[i * 7 * 3 + 5 * 3 + k] =
                    mvec[i * 6 * 3 + 3 * 3 + k]
                        * (areasq2[i * 6 * 3 + 4 * 3 + 1] + areasq2[i * 6 * 3 + 3 * 3 + 2])
                    + mvec[i * 6 * 3 + 2 * 3 + k] * areasq2[i * 6 * 3 + 3 * 3 + 2]
                    + mvec[i * 6 * 3 + 4 * 3 + k] * areasq2[i * 6 * 3 + 4 * 3 + 1];

                curlz[i * 7 * 3 + 6 * 3 + k] = 0.0;
            }
            else {
                curlz[i * 7 * 3 + 0 * 3 + k] =
                    mvec[i * 6 * 3 + 0 * 3 + k]
                        * (areasq2[i * 6 * 3 + 0 * 3 + 0] + areasq2[i * 6 * 3 + 1 * 3 + 0])
                    + mvec[i * 6 * 3 + 1 * 3 + k]
                          * (areasq2[i * 6 * 3 + 1 * 3 + 0] + areasq2[i * 6 * 3 + 2 * 3 + 0])
                    + mvec[i * 6 * 3 + 2 * 3 + k]
                          * (areasq2[i * 6 * 3 + 2 * 3 + 0] + areasq2[i * 6 * 3 + 3 * 3 + 0])
                    + mvec[i * 6 * 3 + 3 * 3 + k]
                          * (areasq2[i * 6 * 3 + 3 * 3 + 0] + areasq2[i * 6 * 3 + 4 * 3 + 0])
                    + mvec[i * 6 * 3 + 4 * 3 + k]
                          * (areasq2[i * 6 * 3 + 4 * 3 + 0] + areasq2[i * 6 * 3 + 5 * 3 + 0])
                    + mvec[i * 6 * 3 + 5 * 3 + k]
                          * (areasq2[i * 6 * 3 + 5 * 3 + 0] + areasq2[i * 6 * 3 + 0 * 3 + 0]);

                curlz[i * 7 * 3 + 1 * 3 + k] =
                    mvec[i * 6 * 3 + 5 * 3 + k]
                        * (areasq2[i * 6 * 3 + 0 * 3 + 1] + areasq2[i * 6 * 3 + 5 * 3 + 2])
                    + mvec[i * 6 * 3 + 4 * 3 + k] * areasq2[i * 6 * 3 + 5 * 3 + 2]
                    + mvec[i * 6 * 3 + 0 * 3 + k] * areasq2[i * 6 * 3 + 0 * 3 + 1];

                curlz[i * 7 * 3 + 2 * 3 + k] =
                    mvec[i * 6 * 3 + 0 * 3 + k]
                        * (areasq2[i * 6 * 3 + 1 * 3 + 1] + areasq2[i * 6 * 3 + 0 * 3 + 2])
                    + mvec[i * 6 * 3 + 5 * 3 + k] * areasq2[i * 6 * 3 + 0 * 3 + 2]
                    + mvec[i * 6 * 3 + 1 * 3 + k] * areasq2[i * 6 * 3 + 1 * 3 + 1];

                curlz[i * 7 * 3 + 3 * 3 + k] =
                    mvec[i * 6 * 3 + 1 * 3 + k]
                        * (areasq2[i * 6 * 3 + 2 * 3 + 1] + areasq2[i * 6 * 3 + 1 * 3 + 2])
                    + mvec[i * 6 * 3 + 0 * 3 + k] * areasq2[i * 6 * 3 + 1 * 3 + 2]
                    + mvec[i * 6 * 3 + 2 * 3 + k] * areasq2[i * 6 * 3 + 2 * 3 + 1];

                curlz[i * 7 * 3 + 4 * 3 + k] =
                    mvec[i * 6 * 3 + 2 * 3 + k]
                        * (areasq2[i * 6 * 3 + 3 * 3 + 1] + areasq2[i * 6 * 3 + 2 * 3 + 2])
                    + mvec[i * 6 * 3 + 1 * 3 + k] * areasq2[i * 6 * 3 + 2 * 3 + 2]
                    + mvec[i * 6 * 3 + 3 * 3 + k] * areasq2[i * 6 * 3 + 3 * 3 + 1];

                curlz[i * 7 * 3 + 5 * 3 + k] =
                    mvec[i * 6 * 3 + 3 * 3 + k]
                        * (areasq2[i * 6 * 3 + 4 * 3 + 1] + areasq2[i * 6 * 3 + 3 * 3 + 2])
                    + mvec[i * 6 * 3 + 2 * 3 + k] * areasq2[i * 6 * 3 + 3 * 3 + 2]
                    + mvec[i * 6 * 3 + 4 * 3 + k] * areasq2[i * 6 * 3 + 4 * 3 + 1];

                curlz[i * 7 * 3 + 6 * 3 + k] =
                    mvec[i * 6 * 3 + 4 * 3 + k]
                        * (areasq2[i * 6 * 3 + 5 * 3 + 1] + areasq2[i * 6 * 3 + 4 * 3 + 2])
                    + mvec[i * 6 * 3 + 3 * 3 + k] * areasq2[i * 6 * 3 + 4 * 3 + 2]
                    + mvec[i * 6 * 3 + 5 * 3 + k] * areasq2[i * 6 * 3 + 5 * 3 + 1];
            }
        }
        for (int j = 0; j < 7; j++)
            for (int k = 0; k < 3; k++)
                curlz[i * 7 * 3 + j * 3 + k] = curlz[i * 7 * 3 + j * 3 + k] / (2.0 * areasT[i]);
    }
    delete[] areasq2;
}

void Icogrid::gra_operator(double *areasT,
                           double *areas,
                           double *grad,
                           double *nvec,
                           int *   pent_ind,
                           int     point_num) {
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
    int    geo;
    double area1, area2, area3, area4, area5, area6;

    double *areasq2;

    areasq2 = new double[6 * 3 * point_num]();

    for (int i = 0; i < point_num; i++) {

        geo = 6; // Hexagons.
        for (int k = 0; k < 12; k++)
            if (i == pent_ind[k])
                geo = 5; // Pentagons.

        area1 = areas[i * 6 * 3 + 0 * 3 + 0] + areas[i * 6 * 3 + 0 * 3 + 1]
                + areas[i * 6 * 3 + 0 * 3 + 2];
        area2 = areas[i * 6 * 3 + 1 * 3 + 0] + areas[i * 6 * 3 + 1 * 3 + 1]
                + areas[i * 6 * 3 + 1 * 3 + 2];
        area3 = areas[i * 6 * 3 + 2 * 3 + 0] + areas[i * 6 * 3 + 2 * 3 + 1]
                + areas[i * 6 * 3 + 2 * 3 + 2];
        area4 = areas[i * 6 * 3 + 3 * 3 + 0] + areas[i * 6 * 3 + 3 * 3 + 1]
                + areas[i * 6 * 3 + 3 * 3 + 2];
        area5 = areas[i * 6 * 3 + 4 * 3 + 0] + areas[i * 6 * 3 + 4 * 3 + 1]
                + areas[i * 6 * 3 + 4 * 3 + 2];

        if (geo == 5)
            area6 = 0;
        else
            area6 = areas[i * 6 * 3 + 5 * 3 + 0] + areas[i * 6 * 3 + 5 * 3 + 1]
                    + areas[i * 6 * 3 + 5 * 3 + 2];

        for (int k = 0; k < 3; k++) {
            areasq2[i * 6 * 3 + 0 * 3 + k] = areas[i * 6 * 3 + 0 * 3 + k] / area1;
            areasq2[i * 6 * 3 + 1 * 3 + k] = areas[i * 6 * 3 + 1 * 3 + k] / area2;
            areasq2[i * 6 * 3 + 2 * 3 + k] = areas[i * 6 * 3 + 2 * 3 + k] / area3;
            areasq2[i * 6 * 3 + 3 * 3 + k] = areas[i * 6 * 3 + 3 * 3 + k] / area4;
            areasq2[i * 6 * 3 + 4 * 3 + k] = areas[i * 6 * 3 + 4 * 3 + k] / area5;
            if (geo == 5)
                areasq2[i * 6 * 3 + 5 * 3 + k] = 0;
            else
                areasq2[i * 6 * 3 + 5 * 3 + k] = areas[i * 6 * 3 + 5 * 3 + k] / area6;
        }
    }

    for (int i = 0; i < point_num; i++) {

        geo = 6; // Hexagons.

        for (int k = 0; k < 12; k++)
            if (i == pent_ind[k])
                geo = 5; // Pentagons.
        for (int k = 0; k < 3; k++) {
            if (geo == 5) { //Pentagons.
                grad[i * 7 * 3 + 0 * 3 + k] =
                    nvec[i * 6 * 3 + 0 * 3 + k]
                        * (areasq2[i * 6 * 3 + 0 * 3 + 0] + areasq2[i * 6 * 3 + 1 * 3 + 0])
                    + nvec[i * 6 * 3 + 1 * 3 + k]
                          * (areasq2[i * 6 * 3 + 1 * 3 + 0] + areasq2[i * 6 * 3 + 2 * 3 + 0])
                    + nvec[i * 6 * 3 + 2 * 3 + k]
                          * (areasq2[i * 6 * 3 + 2 * 3 + 0] + areasq2[i * 6 * 3 + 3 * 3 + 0])
                    + nvec[i * 6 * 3 + 3 * 3 + k]
                          * (areasq2[i * 6 * 3 + 3 * 3 + 0] + areasq2[i * 6 * 3 + 4 * 3 + 0])
                    + nvec[i * 6 * 3 + 4 * 3 + k]
                          * (areasq2[i * 6 * 3 + 4 * 3 + 0] + areasq2[i * 6 * 3 + 0 * 3 + 0]);
                grad[i * 7 * 3 + 0 * 3 + k] =
                    grad[i * 7 * 3 + 0 * 3 + k]
                    - 2
                          * (nvec[i * 6 * 3 + 0 * 3 + k] + nvec[i * 6 * 3 + 1 * 3 + k]
                             + nvec[i * 6 * 3 + 2 * 3 + k] + nvec[i * 6 * 3 + 3 * 3 + k]
                             + nvec[i * 6 * 3 + 4 * 3 + k]);

                grad[i * 7 * 3 + 1 * 3 + k] =
                    nvec[i * 6 * 3 + 4 * 3 + k]
                        * (areasq2[i * 6 * 3 + 0 * 3 + 1] + areasq2[i * 6 * 3 + 4 * 3 + 2])
                    + nvec[i * 6 * 3 + 3 * 3 + k] * areasq2[i * 6 * 3 + 4 * 3 + 2]
                    + nvec[i * 6 * 3 + 0 * 3 + k] * areasq2[i * 6 * 3 + 0 * 3 + 1];

                grad[i * 7 * 3 + 2 * 3 + k] =
                    nvec[i * 6 * 3 + 0 * 3 + k]
                        * (areasq2[i * 6 * 3 + 1 * 3 + 1] + areasq2[i * 6 * 3 + 0 * 3 + 2])
                    + nvec[i * 6 * 3 + 4 * 3 + k] * areasq2[i * 6 * 3 + 0 * 3 + 2]
                    + nvec[i * 6 * 3 + 1 * 3 + k] * areasq2[i * 6 * 3 + 1 * 3 + 1];

                grad[i * 7 * 3 + 3 * 3 + k] =
                    nvec[i * 6 * 3 + 1 * 3 + k]
                        * (areasq2[i * 6 * 3 + 2 * 3 + 1] + areasq2[i * 6 * 3 + 1 * 3 + 2])
                    + nvec[i * 6 * 3 + 0 * 3 + k] * areasq2[i * 6 * 3 + 1 * 3 + 2]
                    + nvec[i * 6 * 3 + 2 * 3 + k] * areasq2[i * 6 * 3 + 2 * 3 + 1];

                grad[i * 7 * 3 + 4 * 3 + k] =
                    nvec[i * 6 * 3 + 2 * 3 + k]
                        * (areasq2[i * 6 * 3 + 3 * 3 + 1] + areasq2[i * 6 * 3 + 2 * 3 + 2])
                    + nvec[i * 6 * 3 + 1 * 3 + k] * areasq2[i * 6 * 3 + 2 * 3 + 2]
                    + nvec[i * 6 * 3 + 3 * 3 + k] * areasq2[i * 6 * 3 + 3 * 3 + 1];

                grad[i * 7 * 3 + 5 * 3 + k] =
                    nvec[i * 6 * 3 + 3 * 3 + k]
                        * (areasq2[i * 6 * 3 + 4 * 3 + 1] + areasq2[i * 6 * 3 + 3 * 3 + 2])
                    + nvec[i * 6 * 3 + 2 * 3 + k] * areasq2[i * 6 * 3 + 3 * 3 + 2]
                    + nvec[i * 6 * 3 + 4 * 3 + k] * areasq2[i * 6 * 3 + 4 * 3 + 1];

                grad[i * 7 * 3 + 6 * 3 + k] = 0.0;
            }
            else {
                grad[i * 7 * 3 + 0 * 3 + k] =
                    nvec[i * 6 * 3 + 0 * 3 + k]
                        * (areasq2[i * 6 * 3 + 0 * 3 + 0] + areasq2[i * 6 * 3 + 1 * 3 + 0])
                    + nvec[i * 6 * 3 + 1 * 3 + k]
                          * (areasq2[i * 6 * 3 + 1 * 3 + 0] + areasq2[i * 6 * 3 + 2 * 3 + 0])
                    + nvec[i * 6 * 3 + 2 * 3 + k]
                          * (areasq2[i * 6 * 3 + 2 * 3 + 0] + areasq2[i * 6 * 3 + 3 * 3 + 0])
                    + nvec[i * 6 * 3 + 3 * 3 + k]
                          * (areasq2[i * 6 * 3 + 3 * 3 + 0] + areasq2[i * 6 * 3 + 4 * 3 + 0])
                    + nvec[i * 6 * 3 + 4 * 3 + k]
                          * (areasq2[i * 6 * 3 + 4 * 3 + 0] + areasq2[i * 6 * 3 + 5 * 3 + 0])
                    + nvec[i * 6 * 3 + 5 * 3 + k]
                          * (areasq2[i * 6 * 3 + 5 * 3 + 0] + areasq2[i * 6 * 3 + 0 * 3 + 0]);
                grad[i * 7 * 3 + 0 * 3 + k] =
                    grad[i * 7 * 3 + 0 * 3 + k]
                    - 2
                          * (nvec[i * 6 * 3 + 0 * 3 + k] + nvec[i * 6 * 3 + 1 * 3 + k]
                             + nvec[i * 6 * 3 + 2 * 3 + k] + nvec[i * 6 * 3 + 3 * 3 + k]
                             + nvec[i * 6 * 3 + 4 * 3 + k] + nvec[i * 6 * 3 + 5 * 3 + k]);

                grad[i * 7 * 3 + 1 * 3 + k] =
                    nvec[i * 6 * 3 + 5 * 3 + k]
                        * (areasq2[i * 6 * 3 + 0 * 3 + 1] + areasq2[i * 6 * 3 + 5 * 3 + 2])
                    + nvec[i * 6 * 3 + 4 * 3 + k] * areasq2[i * 6 * 3 + 5 * 3 + 2]
                    + nvec[i * 6 * 3 + 0 * 3 + k] * areasq2[i * 6 * 3 + 0 * 3 + 1];

                grad[i * 7 * 3 + 2 * 3 + k] =
                    nvec[i * 6 * 3 + 0 * 3 + k]
                        * (areasq2[i * 6 * 3 + 1 * 3 + 1] + areasq2[i * 6 * 3 + 0 * 3 + 2])
                    + nvec[i * 6 * 3 + 5 * 3 + k] * areasq2[i * 6 * 3 + 0 * 3 + 2]
                    + nvec[i * 6 * 3 + 1 * 3 + k] * areasq2[i * 6 * 3 + 1 * 3 + 1];

                grad[i * 7 * 3 + 3 * 3 + k] =
                    nvec[i * 6 * 3 + 1 * 3 + k]
                        * (areasq2[i * 6 * 3 + 2 * 3 + 1] + areasq2[i * 6 * 3 + 1 * 3 + 2])
                    + nvec[i * 6 * 3 + 0 * 3 + k] * areasq2[i * 6 * 3 + 1 * 3 + 2]
                    + nvec[i * 6 * 3 + 2 * 3 + k] * areasq2[i * 6 * 3 + 2 * 3 + 1];

                grad[i * 7 * 3 + 4 * 3 + k] =
                    nvec[i * 6 * 3 + 2 * 3 + k]
                        * (areasq2[i * 6 * 3 + 3 * 3 + 1] + areasq2[i * 6 * 3 + 2 * 3 + 2])
                    + nvec[i * 6 * 3 + 1 * 3 + k] * areasq2[i * 6 * 3 + 2 * 3 + 2]
                    + nvec[i * 6 * 3 + 3 * 3 + k] * areasq2[i * 6 * 3 + 3 * 3 + 1];

                grad[i * 7 * 3 + 5 * 3 + k] =
                    nvec[i * 6 * 3 + 3 * 3 + k]
                        * (areasq2[i * 6 * 3 + 4 * 3 + 1] + areasq2[i * 6 * 3 + 3 * 3 + 2])
                    + nvec[i * 6 * 3 + 2 * 3 + k] * areasq2[i * 6 * 3 + 3 * 3 + 2]
                    + nvec[i * 6 * 3 + 4 * 3 + k] * areasq2[i * 6 * 3 + 4 * 3 + 1];

                grad[i * 7 * 3 + 6 * 3 + k] =
                    nvec[i * 6 * 3 + 4 * 3 + k]
                        * (areasq2[i * 6 * 3 + 5 * 3 + 1] + areasq2[i * 6 * 3 + 4 * 3 + 2])
                    + nvec[i * 6 * 3 + 3 * 3 + k] * areasq2[i * 6 * 3 + 4 * 3 + 2]
                    + nvec[i * 6 * 3 + 5 * 3 + k] * areasq2[i * 6 * 3 + 5 * 3 + 1];
            }
        }
        for (int j = 0; j < 7; j++)
            for (int k = 0; k < 3; k++)
                grad[i * 7 * 3 + j * 3 + k] = grad[i * 7 * 3 + j * 3 + k] / (2.0 * areasT[i]);
    }
    delete[] areasq2;
}

void Icogrid::zonal_mean_tab_f(int *   zonal_mean_tab,
                               double *lonlat,
                               int     nlat,
                               int     point_num,
                               int *   max_count) {

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

    double  des_lat = M_PI / nlat;
    int *   count_num;
    double *lat_array;
    // int     max_count = 0;

    count_num = new int[nlat]();
    lat_array = new double[nlat + 1]();

    for (int i = 0; i < point_num; i++)
        for (int k = 0; k < 2; k++)
            zonal_mean_tab[i * 3 + k] = 0;

    for (int j = 0; j < nlat; j++)
        count_num[j] = 0;
    lat_array[0] = -M_PI / 2.0;
    for (int j = 1; j < nlat + 1; j++)
        lat_array[j] = lat_array[j - 1] + des_lat;
    for (int i = 0; i < point_num; i++) {
        for (int j = 0; j < nlat; j++) {
            if (lonlat[i * 2 + 1] >= lat_array[j] && lonlat[i * 2 + 1] < lat_array[j + 1]) {
                zonal_mean_tab[i * 3]     = j;
                zonal_mean_tab[i * 3 + 2] = count_num[j];
                count_num[j]              = count_num[j] + 1;
                break;
            }
        }
        if (lonlat[i * 2 + 1] >= lat_array[nlat]) {
            zonal_mean_tab[i * 3]     = nlat - 1;
            zonal_mean_tab[i * 3 + 2] = count_num[nlat - 1];
            count_num[nlat - 1]       = count_num[nlat - 1] + 1;
        }
    }
    for (int i = 0; i < point_num; i++) {
        int ind                   = zonal_mean_tab[i * 3];
        zonal_mean_tab[i * 3 + 1] = count_num[ind];
        // log::printf("ind, count = %d, %d\n", ind, count_num[ind]);
        if (count_num[ind] > *max_count)
            *max_count = count_num[ind];
    }
    *max_count = pow(2, ceil(log(*max_count) / log(2)));

    // log::printf("max count = %d\n", max_count);
    delete[] count_num;
    delete[] lat_array;
}

void Icogrid::free_memory() {
    log::printf("Freeing Grid memory.\n");

    // Free memory that's not passed and freed by ESP.

    free(nvec);
    free(point_xyz);
    free(pent_ind);
    free(halo);
    free(point_xyzq);
    free(areas);

    // Allocated data passed over to ESP object
    // Grid serves only for setup, but should live all throughout sim.
    free(point_local);
    free(maps);
    free(lonlat);
    free(Altitude);
    free(Altitudeh);
    free(nvecoa);
    free(nvecti);
    free(nvecte);
    free(areasT);

    free(areasTr);
    free(div);
    free(grad);
    free(func_r);
}


//END OF GRID.CU
