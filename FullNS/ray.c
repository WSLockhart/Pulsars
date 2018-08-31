/*
*****************************************************************
* 03/25/10: =====Version 1=====
*   (DP)    This program calculates the trajectories of photons
*           in a quasi-Kerr spacetime. The metric is that of
*           Glampedakis & Babak (2006, Class Quant Grav 23, 4167).
*           In this first incarnation, the emission is assumed to
*           originate from a geometrically thin accretion disk,
*           which lies on the equatorial plane of the spacetime and
*           extends down to the ISCO or another user-defined inner
*           disk radius.
*
* 12/06/10: =====Version 2====
*   (DP)    This is a different direction for the code, where the
*           emission is assumed to originate from the surface of
*           a moderately rotating neutron star.
*
* 12/04/11: =====Version 2.1====
*   (DP)    This version of the code calculates the lightcurve from
*           a hot spot on the spinning neutron star
*
* 18/06/12: =====Version 2.2====
*   (DP)    Calculating the lightcurves with interpolation proved
*           very noise for small spots. The code now find the region
*           on the image plane where the emission appears and then,
*           instead of interpolating, performs ray tracing with a fine
*           resolution.
*
* 29/07/12: =====Version 2.3===
*   (DP)    Changed two inputs to the code. One is the quadrupole moment
*           of the neutron star, which is now given as eta, so that the
*           deviation from Kerr is epsilon=eta*a_bh*a_bh. The second is
*           the equatorial radius of the neutron star, which is now
*           given as the circumferential radius, from which the coordinate
*           radius is given. Also, the calculation of the beaming angle is
*           redone, to be consistent with the Hartle-Thorne metric.
***************************************************************** */

#include "definitions.h"                // File with useful constants

/* ****************************************************************
 Global Variable and Parameter Declarations
     repeated in header file global.h to be included in all
     files linked together with this main */
// Input parameters
double a_bh;                            // neutron-star angular momentum
double epsilon;                         // deviation from Kerr quadrupole
double rns;                             // neutron star radius (in km)
double mns;                             // neutron star mass (in Msolar)
double fns;                             // neutron star spin (1/P in Hz)
double theta0;                          // inclination of observer
double the_spot;                        // lattitude of spot
double rho_spot;                        // half opening angle of spot
double Tbb;                             // Temperature on NS surface
// Other useful parameters
double Omegans;                         // angular velocity of NS (in 1/M)
double rcirc_ns;                        // the circumferential radius of the star
/* **************************************************************** */

/*
****************************************************************
   Protypes of external functions called by the main program
**************************************************************** */
extern void readinput();                          // from io.c
extern void set_energy_grid(double Energy[Nen]);  // from integrals.c
extern void set_time_grid(double Time[Nt]);       // from time.c
extern void time_spectrum(double Energy[], double Time[]);
                                                  // from time.c
//extern void output(double Energy[], double Time[], double dyn_spectra[Nt][Nen]);
                                                  // from io.c

//extern void output(double Energy[], double Transfer[]); // prints the result
//extern void debug(void);                 // debugging subroutine
//extern double v_ns(double theta, double velocity[5]);
extern double r_shape(double theta);
extern void metric(double r, double theta, double gmunu[5]);
extern double Circ_to_Coord(double rcirc_ns);         // from file metric.c
/* **************************************************************** */



/* ****************************************************************
    Main Program
    *************************************************************** */
int main(void)
{
  double Energy[Nen];                   // grid of photon energies for spectrum
  double Time[Nt];                      // grid of times for lightcurve
  //  double dyn_spectra[Nt][Nen];      // dynamical energy spectra
  int index;
   
  // Reads the input parameters from file <in>
  readinput(); //from io.c
  printf("! Read input data\n");

  //  Set the energy grid ()
  set_energy_grid(Energy); //from integrals.c
  printf("! Set energy grid\n");

  //  Set the time grid
  set_time_grid(Time); //from time.c
  printf("! Set time grid\n");

  time_spectrum(Energy,Time);

   // Print the output
  //  output(Energy,Time,dyn_spectra);

  return 0;
}
