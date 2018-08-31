/* 
************************************************************
             File: global.h
   Header file with "extern" decleration of all global
   variables. Global variables are defined in the main
   source file. This header file is to be included in all
   other source code files.
************************************************************
*/

// Input parameters
extern double a_bh;                      // neutron-star angular momentum
extern double epsilon;                   // deviation from Kerr quadrupole
extern double rns;                       // neutron star radius (in km)
extern double mns;                       // neutron star mass (in Msolar)
extern double fns;                       // neutron star spin (1/P in Hz)
extern double theta0;                    // inclination of observer
extern double the_spot;                  // lattitude of spot
extern double rho_spot;                  // half opening angle of spot
extern double Tbb;                       // Temperature on NS surface

// Other useful parameters
extern double Omegans;                   // angular velocity of NS (in 1/M)
extern double rcirc_ns;                  // the circumferential radius of the star
