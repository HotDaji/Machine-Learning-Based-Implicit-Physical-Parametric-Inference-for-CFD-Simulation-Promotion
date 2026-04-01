#ifndef PARAMETER_H
#define PARAMETER_H
#include <cmath>
#include<vector>
#include <iostream>

typedef __time_t time_t;
//#define SINGLE
#define D2Q9
#define DirectForce 
//#define Stiffness
//#define Dirac_Delta
#define SRTLBM
#define SmagorinskyByHybrid
//#define RRLBM
// #define MRT
#ifdef SINGLE
    typedef float REAL;
#else
    typedef double REAL;
#endif

// index:   0  1  2  3  4  5  6  7  8
//  ----------------------------------
//  x:       0 +1 -1  0  0 +1 -1 +1 -1
//  y:       0  0  0 +1 -1 +1 -1 -1 +1

//  8 3 5  ^y
//   \|/   |   x
//  2-0-1   --->
//   /|\
//  6 4 7






#endif