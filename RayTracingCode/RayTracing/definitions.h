/*
*****************************************************************
* File with useful definitions
***************************************************************** */

/* Headers */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/* index to coordinate correspondence for radial vectors and momenta
   (but not for the Runge-Kutta arrays */
#define tco 0
#define rco 1
#define theco 2
#define phico 3

/* Global Parameters for the Image Calculation */
                // Nx and Ny define the size of the coarse grid
#define Nx 385   // # of grid points along the x-axis on the image plane [default 385]
#define Ny 385  // # of grid points along the y-axis on the image plane
                // NB!!!: Nx and Ny have to be odd integers
                // Nx_h and Ny_h define the size of the fine grid
#define Nx_h 256 // # of grid points along the x-axis on the image plane [default 512]
#define Ny_h 256 // # of grid points along the y-axis on the image plane

#define frac 1.2 // the fine grid covers "frac" times the spot size

#define Xmax 18 // Size of calculated image along the x-axis in M [default 18]
#define Ymax 18 // Size of calculated image along the y-axis in M
#define Dobs 100028.0// Distance of the observer in M

/* Parameters for the accuracy of the calculation */
#define step_factor 32.0 // fraction of minimum timescale for RK timestep
                          // even 32.0 is mostly adequate

/* Parameters about the metric */
#define RMIN_TROUBLE 2.0     // radius inside which the metric is ill-behaved

/* Intersection-not-found code */
#define NOT_FOUND -9999.0

/* Parameters about the spectra */
/* all energies are in units of the effective temperature on the surface */
#define Nen 128          // number of energy grid points for spectra [default 128]
#define Emin 0.01        // minimum photon Energy for grid
#define Emax 5.0          // maximum photon Energy for grid

/* Parameters about time-dependence spectra */
#define Nt 65            // number of time grid points [default 65]
#define Time_max 1         // maximum number of spin periods to be used

#define Dist 0.2           //distance to the source
