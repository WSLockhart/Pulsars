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
* 12/06/10: =====Version 2=====
*   (DP)    This is a different direction for the code, where the
*           emission is assumed to originate from the surface of
*           a moderately rotating neutron star.
*
* 12/04/11: =====Version 2.1=====
*   (DP)    This version of the code calculates the lightcurve from
*           a hot spot on the spinning neutron star
*
* 18/06/12: =====Version 2.2=====
*   (DP)    Calculating the lightcurves with interpolation proved
*           very noise for small spots. The code now find the region
*           on the image plane where the emission appears and then,
*           instead of interpolating, performs ray tracing with a fine
*           resolution.
*
* 29/07/12: =====Version 2.3=====
*   (DP)    Changed two inputs to the code. One is the quadrupole moment
*           of the neutron star, which is now given as eta, so that the
*           deviation from Kerr is epsilon=eta*a_bh*a_bh. The second is
*           the equatorial radius of the neutron star, which is now
*           given as the circumferential radius, from which the coordinate
*           radius is given. Also, the calculation of the beaming angle is
*           redone, to be consistent with the Hartle-Thorne metric.


* 15/03/19: NB:=====Version 3=====
*   (WL)    Will Lockhart here. I have modified the code to make the hotspots
            more realistic, incorporating recent advances in magnetospheric
            theory from Gralla (2017). I've added a new script, 'jmath.c'.
            This contains the function current_density() which returns the
            norm of the 4-current at a point (theta,phi) on the surface of
            the neutron star. emissivity() has been modified to call this
            script and return a specific intensity that is a function of the
            current density.

            The script 'jtemp.c' is also new. This is for a specific calculation
            for my first paper (separate from the main program) to compute the
            ratio of the average temperatures between the two hemispheres.

            I have also changed how the code receives input parameters.
            Rather than reading them from a file ('in'), the inputs are
            passed to ray.c when the program is run. This allows us to
            explore the parameter space with a computer cluster. We no longer
            use the parameter rho_spot, since the size and shape of the hotspot
            are now determined by jmath instead. We've introduced one new input,
            'qns', the magnetic quadrupole-to-dipole ratio.



***************************************************************** */

#include "definitions.h"                // File with useful constants
#include "global.h"                     // Global variable headers

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
double qns;                             // quadrupole-to-dipole ratio on NS surface

// Other useful parameters
double Omegans;                         // angular velocity of NS (in 1/M)
double rcirc_ns;                        // the circumferential radius of the star
/* **************************************************************** */

/*
****************************************************************
   Protypes of external functions called by the main program
**************************************************************** */
// extern void readinput();                          // from io.c
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
extern double Circ_to_Coord(double rcirc_ns);         // from file ns.c
/* **************************************************************** */



/* ****************************************************************
    Main Program
    *************************************************************** */
int main(int argc, char **argv)
{
  if (argc==1) {
    fprintf(stderr, "! No arguments were passed!\n");
  } else {
    if (argc!=8){
      fprintf(stderr, "! Wrong number of arguments passed!\n");
    }
  }
  int i;
  double input_values[argc-1];
  for(i=1; i<argc; i+=1) {
    input_values[i] = strtod(argv[i], NULL);
  }
  // Reads the input parameters from file <in>
  //readinput(); //from io.c

  // read inputs
  mns=input_values[1];
  rcirc_ns=input_values[2];
  fns=input_values[3];
  theta0=input_values[4]*(M_PI/180.);
  the_spot=input_values[5]*(M_PI/180.);
  Tbb=input_values[6];
  qns=input_values[7];
  // calculate the coordinate ns radius, given the value of the circumferential ns radius
  rns=Circ_to_Coord(rcirc_ns);
  // calculate the dimensionless angular frequency of the neutron star
  Omegans=2.*M_PI*4.92549e-6*fns*mns;
  // the coefficient is 2*pi*G*Msolar/c^3 to convert Hz to units of M

  double Energy[Nen];      // grid of photon energies for spectrum
  double Time[Nt];         // grid of times for lightcurve
  set_energy_grid(Energy); // from integrals.c
  // printf("! Set energy grid\n");
  set_time_grid(Time); //from time.c
  // printf("! Set time grid\n");

  printf("-----------------------------------\n");
  printf("INPUT PARAMS:\n mns=%.1f Msolar\n rns=%.1f km\n fns=%.1f Hz\n obs_angle=%dº\n spot_angle=%dº\n Tbb=%.2f KeV\n q=%.1f\n", \
      mns,rns,fns,(int)(theta0*(180./M_PI)),(int)(the_spot*(180./M_PI)),Tbb,qns);
  printf("CODE PARAMS:\n Nt=%d\n Nen=%d\n Nx_h=%d,Ny_h=%d\n",Nt,Nen,Nx_h,Ny_h);
  printf("-----------------------------------\n");

  time_spectrum(Energy,Time);  //main program

  return 0;
}
