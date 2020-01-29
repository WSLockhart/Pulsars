/*
************************************************************
             File: ns.c
   Functions that calculate the properties of the neutron star
************************************************************
*/

#include "definitions.h"      // File with useful definitions and headers
#include "global.h"           // Global variable headers

/*
****************************************************************
   Protypes of external functions used in this file
**************************************************************** */
extern void metric(double r, double theta, double gmunu[5]);           //from metric.c
extern double rootbis(double xmin, double xmax,double (*foo)(double));   //from math.c
double current_density(double theta, double phi);  //from jmath.c
/* **************************************************************** */

double eq_circ_to_coord(double rns)
/* July 29, 2012 (DP)
   this is the equation
   sqrt(g_phi_phi(r=rns,theta=PI/2))-rcirc_ns=0
   to be solved by the Circ_to_Coord subroutine. The circumferential
   radius rcirc_ns is given as a global variable
*/
{
  double gmunu[5];

  // calculate the metric on the equatorial plane
  metric(rns,M_PI/2,gmunu);

  // return the LHS of the equation
  return sqrt(gmunu[phico])-rcirc_ns;

}

double Circ_to_Coord(void)
/* July 29, 2012 (DP)
   Given the circumferential radius of the neutron star, rcirc_ns,
   which is stored as a global variable, the routine
   solves the equation
      sqrt(g_phi_phi(r=rns,theta=PI/2))=rcirc_ns
   to find the equatorial coordinate radius of the star. It uses  a
   bisection technique. The parameters fmin and fmax define the
   fractions of the circumferential radius between which the bisection
   routine searches for a solution.
*/
{
  double radius;         // the coordinate neutron star radius
#define fmin 0.95
#define fmax 1.05

  //  radius=rootbis(fmin*rcirc_ns,fmax*rcirc_ns,eq_circ_to_coord);

  //  return radius;

  // spherical star
  return rcirc_ns;

}

double r_shape(double theta)
/* December 7 2010 (DP)
   Returns the radius of a neutron star (in M) as a function of the polar
   angle theta with respect to the rotation axis. The equatorial
   radius of the neutron star is given in the global variable rns.
   The shape of the rotating neutron star surface is calculated using
   the approximation of Morsink et al. 2007, ApJ, 663, 1244, truncated to
   second order, in order to be consistent with the truncation of the
   spacetime to quadrupole order

   December 6, 2011 (DP)
   Corrects the fact that Morsink et al. 2007 use an isotropic coordinate
   sysem and not the quasi-Schwarzschild that the code is using.
   Note: Does it?
*/
{
  double zeta, epsilon_ns;                   // auxiliary parameters
  double a0,a2;
  double radius;                          // polar angle dependent radius
  double x;                               // argument of Legendre polynomials

  // zeta=GM/Rc^2
  zeta=1.47663*mns/rns;

  // epsilon=omega^2*r^2/c^2/zeta
  epsilon_ns=4.39257e-10*fns*fns*rns*rns/zeta;
  //  a0=-0.18*epsilon_ns+0.23*zeta*epsilon_ns-0.05*epsilon_ns*epsilon_ns;
  //a2=-0.39*epsilon_ns+0.29*zeta*epsilon_ns+0.13*epsilon_ns*epsilon_ns;

  // zero oblateness. Change if you want it oblate!!!
  a0=0.0;
  a2=0.0;
  //x=cos(theta);                          // for Legendre polynomials
  //radius=(1.0+a0+0.5*a2*(3.*x*x-1.))/zeta; // approximate radius in GM/c^2
  radius=1./zeta;

  return radius;
}

void vel_ns(double theta, double velocity[])
/* December 7, 2010 (DP)
   Returns in the array velocity[4] the four components of the
   contravariant 4-velocity of matter on the stellar surface.
   The ordering of the components in the array are 0:t, 1:r,
   2:theta, 3:phi. Only the t- and phi- components are non-zero.
   The Glampedakis & Babak (2006) metric are assumed and
   the parameters of the neutron star are in global variables.
*/
{
  double radius;                         // coordinate radius of surface (in M)
  double gmunu[5];                       // metric elements
  double norm;                           // normalization

  // find the coordinate radius of the NS surface at that polar angle
  radius=r_shape(theta);

  // calculate the metric elements
  metric(radius,theta,gmunu);

  norm=pow(-gmunu[tco]-Omegans*(gmunu[phico]*Omegans+2.*gmunu[4]),-0.5);

  velocity[tco]=norm;
  velocity[rco]=0.0;
  velocity[theco]=0.0;
  velocity[phico]=Omegans*norm;

  return;
}

double ksi(double theta, double phi, double phi0)
/* December 6, 2010 (DP)
   Returns the angular distance between the location on the neutron star
   surface with lattitude and longitude theta and phi from the center of the
   hot spot. phi0 is the azimuth of the center of the hot spot. The
   latitude of the hot spot the_spot is a global variable
*/
{
  double cosksi;       // cosine of ksi


  // calculate the angle between the center of the spot and the current location
  cosksi=cos(phi0)*cos(the_spot)*cos(phi)*sin(theta)+
    sin(phi0)*cos(the_spot)*sin(phi)*sin(theta)+
    sin(the_spot)*cos(theta);

  return acos(cosksi);

}

double emissivity(double energy, double angle_phi, double angle_theta,
		  double beamangle, double phi0)
/* December 6, 2010 (DP)
   Returns the monoenergetic specific intensity that emerges
   from each point on the surface of the neutron star. theta and
   phi are the polar angle and azimuth of the location with respect
   to the rotation axis of the star. beamangle is the locally
   measured angle of emission with respect to the local normal.
   everything is now isotropic

   December 5, 2011 (DP)
   Added a prescription for a uniform hot spot. The half-opening angle,
   rho_spot, of the hot spot is in the global variables.  phi0 is the
   azimuth of the fiducial plane

   July 15, 2012 (DP)
   The emissivity is now photon energy dependent. It is a blackbody,
   with all the energies in keV. The blackbody temperature, Tbb, is in the
   global variables and is read from the input file.

   July 22, 2012 (DP)
   The beaming can be either isotropic or given by the Hopf function;
   for now, the switch is performed by adding/removing appropriately
   the relevant lines of code.
*/
{

  double intensity;
  // double distance=ksi(angle_theta,angle_phi,phi0); // angle between center of spot and location
  // // figure out if it is within the spot
  // rho_spot = 15*(M_PI/180.0);
  // if (distance<rho_spot)
    {
    /*
       July 11, 2018 (WL)
       This is where Sam's model goes! We will set the surface intensity
       to be a blackbody spectrum whose temperature, Tbbsur, scales with
       the Gralla current density as j~T^(1/4).
       This will override rho_spot in determining the hotspot shape and size.
    */

      double Jsquared = current_density(angle_theta,angle_phi-phi0);
      //NOTE: must pass phi-phi0 to account for the star's rotation!

      if (Jsquared>0)
      {
          // (1) Blackbody:
          double Tbbsur = Tbb*8.09*pow(Jsquared,(1/8.)); //the factor of ~8 comes about
          //because the max value of J^1/4 for an aligned dipole is ~0.12. Putting in this factor
          //callibrates it so that temperatures are actually Tbb;
          intensity=energy*energy*energy/(exp(energy/Tbbsur)-1.0);
      }
      else{
          intensity=0;
      }

      //BEAMING FACTOR
      //NB: beamangle is actually the COSINE of the beaming angle (see image.c -> geodesic)
      double beamh = 0;  // beamh=0 means isotropic
      intensity *= (1 + beamh*beamangle*beamangle);

      // gaussian
      //intensity=exp(-(energy-Tbb)*(energy-Tbb)/2./1.e-6)/sqrt(2.*M_PI*1.e-6);
      return intensity;
    }
  // else
  //   {
  //     return 0.0;
  //   }
}
