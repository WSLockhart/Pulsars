/* 
************************************************************
             File: metric.c
   Functions that return the various metric elements
   and Christoffel symbols for the spacetime considered.
   This version incorporates the quasi-Kerr metric of 
   Glampedakis & Babak (2006, Class Quant Grav 23, 4167)
************************************************************

AS A DEBUG TOOL, the Schwarzschild Metric is hardwired
*/

#include "definitions.h"      // File with useful definitions and headers
#include "global.h"           // Global variable headers

void metric(double r, double theta, double gmunu[5])
/* March 25, 2010 (DP)
   Function that returns the metric components of the
   spacetime. The coordinates are 'r' and 'theta', while 
   'phi' does not enter because of azimuthal symmetric. The
   properties of the black hole are encoded in the global 
   variables a_bh and epsilon. The non-zero metric components 
   are returned in the array gmunu[5] ordered as 
      gmunu[0]=g_tt
      gmunu[1]=g_rr
      gmunu[2]=g_thth
      gmunu[3]=g_phiphi
      gmunu[4]=g_phit
*/
{
  // Shorthands to be used in the metric elements
  double Sigma,Delta;
  
  // This is g_tt
  gmunu[0]=-(1.-2./r);
  
  // This is g_rr
  gmunu[1]=1./(1.-2./r);
  
  // This is g_thth
  gmunu[2]=r*r;
  
  // This is g_phiphi
  gmunu[3]=r*r*sin(theta)*sin(theta);
  
  // This is g_tphi
  gmunu[4]=0.0;
  
  return;
  // all done
}

void gammas(double r, double theta, double Gamma[4][4][4])
/* March 30, 2010 (DP), based on the subroutine by TJ.
   Function that returns the Christoffel symbols for the
   spacetime. The coordinates are 'r' and 'theta', while 
   'phi' does not enter because of azimuthal symmetric. The
   properties of the black hole are encoded in the global
   variable a_bh and epsilon. The Christoffel symbols are 
   returned in the array Gamma[4][4][4]. The ordering of 
   the indices for the symbols is
   0=t, 1=r, 2=theta, 3=phi */
{
  // Shorthands to be used for the Christoffel symbols
  double sintheta, costheta, sintheta2, costheta2,sincos;
  double r2, a2, Delta, Sigma, Sigma2, mSigma;
  double ra2, r3, term_2M_r, term2, Logterm, c2s, poly5, poly6, poly7, poly8;

  // The following shorthands should be self explanatory
  // ... about angles
  
  Gamma[0][1][0] = 1./(r*r-2.*r);
  Gamma[0][0][1] = Gamma[0][1][0];    // Symmetry of lower two indices

  Gamma[0][2][0] = 0.;
  Gamma[0][0][2] = Gamma[0][2][0];    // Symmetry of lower two indices
  
  Gamma[0][3][1] = 0.;
  Gamma[0][1][3] = Gamma[0][3][1];    // Symmetry of lower two indices

  Gamma[0][3][2] = 0.;
  Gamma[0][2][3] = Gamma[0][3][2];    // Symmetry of lower two indices

  Gamma[1][0][0] = (-2.+r)/(r*r*r);

  Gamma[1][3][0] = 0.;
  Gamma[1][0][3] = Gamma[1][3][0];    // Symmetry of lower two indices

  Gamma[1][1][1] = 1./(2.*r-r*r);
  
  Gamma[1][2][1] = 0.0;
  Gamma[1][1][2] = Gamma[1][2][1];    // Symmetry of lower two indices

  Gamma[1][2][2] = 2.0-r;
  
  Gamma[1][3][3] = (2.0-r)*sin(theta)*sin(theta);

  Gamma[2][0][0] = 0.0;

  Gamma[2][3][0] = 0.0;
  Gamma[2][0][3] = Gamma[2][3][0];    // Symmetry of lower two indices
  
  Gamma[2][1][1] = 0.0;

  Gamma[2][2][1] = 1./r;
  Gamma[2][1][2] = Gamma[2][2][1];    // Symmetry of lower two indices
  
  Gamma[2][2][2] = 0.0;
  
  Gamma[2][3][3] = -sin(theta)*cos(theta);

  Gamma[3][1][0] = 0.0;
  Gamma[3][0][1] = Gamma[3][1][0];    // Symmetry of lower two indices
  
  Gamma[3][2][0] = 0.0;
  Gamma[3][0][2] = Gamma[3][2][0];    // Symmetry of lower two indices

  Gamma[3][3][1] = 1./r;
  Gamma[3][1][3] = Gamma[3][3][1];    // Symmetry of lower two indices

  Gamma[3][3][2] = cos(theta)/sin(theta);
  Gamma[3][2][3] = Gamma[3][3][2];    // Symmetry of lower two indices

  return;
}

