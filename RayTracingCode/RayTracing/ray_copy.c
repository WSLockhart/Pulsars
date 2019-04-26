

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
