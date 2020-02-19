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
// Operations on vector data types
//
//
// Description:-
//
//
// Known limitations:
//   - None.
//
// Known issues:
//   - None.
//
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

#pragma once

#include <math.h>

// negation
inline __host__ __device__ double3 operator-(const double3 &a) {
    return make_double3(-a.x, -a.y, -a.z);
}

//adition
inline __host__ __device__ double3 operator+(const double3 &a, const double3 &b) {
    return make_double3(a.x + b.x, a.y + b.y, a.z + b.z);
}

inline __host__ __device__ double3 &operator+=(double3 &a, const double3 &b) {
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;

    return a;
}

inline __host__ __device__ double3 operator+(const double3 &a, const double &b) {
    return make_double3(a.x + b, a.y + b, a.z + b);
}
inline __host__ __device__ double3 operator+(const double &b, const double3 &a) {
    return make_double3(a.x + b, a.y + b, a.z + b);
}
inline __host__ __device__ double3 &operator+=(double3 &a, const double &b) {
    a.x += b;
    a.y += b;
    a.z += b;

    return a;
}

// subtraction
inline __host__ __device__ double3 operator-(const double3 &a, const double3 &b) {
    return make_double3(a.x - b.x, a.y - b.y, a.z - b.z);
}

inline __host__ __device__ double3 &operator-=(double3 &a, const double3 &b) {
    a.x -= b.x;
    a.y -= b.y;
    a.z -= b.z;

    return a;
}

inline __host__ __device__ double3 operator-(const double3 &a, const double &b) {
    return make_double3(a.x - b, a.y - b, a.z - b);
}
inline __host__ __device__ double3 operator-(const double &b, const double3 &a) {
    return make_double3(b - a.x, b - a.y, b - a.z);
}
inline __host__ __device__ double3 &operator-=(double3 &a, const double &b) {
    a.x -= b;
    a.y -= b;
    a.z -= b;

    return a;
}

// multiplication
inline __host__ __device__ double3 operator*(const double3 &a, const double &b) {
    return make_double3(a.x * b, a.y * b, a.z * b);
}
inline __host__ __device__ double3 operator*(const double b, const double3 &a) {
    return make_double3(b * a.x, b * a.y, b * a.z);
}
inline __host__ __device__ double3 &operator*=(double3 &a, const double &b) {
    a.x *= b;
    a.y *= b;
    a.z *= b;

    return a;
}

// division
inline __host__ __device__ double3 operator/(const double3 &a, const double &b) {
    return make_double3(a.x / b, a.y / b, a.z / b);
}

inline __host__ __device__ double3 &operator/=(double3 &a, const double &b) {
    a.x /= b;
    a.y /= b;
    a.z /= b;

    return a;
}

inline __host__ __device__ double dot(const double3 &a, const double3 &b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline __host__ __device__ double length(const double3 &a) {
    return sqrt(dot(a, a));
}

inline __host__ __device__ double3 normalize(const double3 &a) {
    return a / length(a);
}

inline __host__ __device__ double3 cross(const double3 &a, const double3 &b) {
    return make_double3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}

__host__ __device__ double3 sort_add3(const double3 &a, const double3 &b, const double3 &c) {
    // Sort three vectors for triangle vertices, then add them in a predictable order
    // Exploit the fact that three vertices cannot all have the same z
    // and if two have the same z, x and y must be different
    // x = r * cos(lon) * cos(lat)
    // y = r * sin(lon) * cos(lat)
    // z = r * sin(lat)

    double3 w1, w2, w3;

    if ((a.z == b.z) || (a.z == c.z) || (b.z == c.z)) {
        // two of the vectors have the same z/latitude (unlikely)
        double3 odd1, same1, same2;
        if (a.z == b.z) { // first, find the odd one out
            odd1  = c;
            same1 = a;
            same2 = b;
        }
        else if (a.z == c.z) {
            odd1  = b;
            same1 = a;
            same2 = c;
        }
        else if (b.z == c.z) {
            odd1  = a;
            same1 = b;
            same2 = c;
        } // other choices should not happen

        if (odd1.z < same1.z) { //odd1 < same1/same2
            w1 = odd1;
            if (same1.x < same2.x) { //odd1 < same1 < same2
                w2 = same1;
                w3 = same2;
            }
            else { //odd1 < same2 < same1
                w2 = same2;
                w3 = same1;
            }
        }
        else { //same1/same2 < odd1
            w3 = odd1;
            if (same1.x < same2.x) { //same1 < same2 < odd1
                w1 = same1;
                w2 = same2;
            }
            else { //same2 < same1 < odd1
                w1 = same2;
                w2 = same1;
            }
        }
    }


    else {
        // all z/latitude values are unique
        if (a.z < b.z) {         // abc or acb or cab
            if (a.z < c.z) {     // abc or acb
                if (b.z < c.z) { //a b c
                    w1 = a;
                    w2 = b;
                    w3 = c;
                }
                else { //a c b
                    w1 = a;
                    w2 = c;
                    w3 = b;
                }
            }
            else { //c a b
                w1 = c;
                w2 = a;
                w3 = b;
            }
        }
        else {               // bac or bca or cba
            if (a.z < c.z) { // b a c
                w1 = b;
                w2 = a;
                w3 = c;
            }
            else {               // bca or cba
                if (b.z < c.z) { //b c a
                    w1 = b;
                    w2 = c;
                    w3 = a;
                }
                else { // c b a
                    w1 = c;
                    w2 = b;
                    w3 = a;
                }
            }
        }
    }

    double3 sum3;
    sum3 = w1 + w2 + w3;
    return sum3;
}
