/*
************************************************************
             File: io.c
   Functions that handle the input of parameters and the
   output of the results
************************************************************
*/

#include "definitions.h"      // File with useful definitions and headers
#include "global.h"           // Global variable headers

/*
****************************************************************
   Protypes of external functions used in this file
**************************************************************** */
extern double Circ_to_Coord(double rcirc_ns);         // from file metric.c
/* **************************************************************** */


void readinput(void)
/* December 8, 2010 (DP)
   Function that reads the input from a file named "in"
   in the current directory. The file is expected to have the following
   set of parameters (in order):
      a_bh   : dimensionless spin of neutron star (0<=a<1)
      eta    : the dimensionless deviation of the spacetime quadrupole from
               the Kerr value is q=-(1+eta)*a_bh^2
      rns    : equatorial radius of neutron star in km
      mns    : gravitational mass of neutron star in solar masses
      fns    : spin of neutron star (1/Period) in Hz
      theta0 : inclination of the observer (in degrees) with respect
               to the spin axis of the neutron star
      the_spot : lattitude of spot
      rho_spot : half-opening angle of spot
      Tbb      : Temperature of the surface emission in keV

    The function also calculates the spin angular velocity of the neutron star
    Omegans. All parameters are then stored in global variables
*/
{
  FILE *in;                            // input file

  float readp;                         // float to read single precession
  double eta;                          // second input parameter that is only use to calculate
                                       // epsilon

  in=fopen("in","r");                  // open file for input parameters

  if (in==NULL)                        // error handling
    {
      printf("Error in opening input file <in>\n");
      return;
    }

  fscanf(in,"%e",&readp);              // first parameter is spin
  a_bh=readp;
  fscanf(in,"%e",&readp);              // second parameter is eta
  eta=readp;
  epsilon=eta*a_bh*a_bh;               // use eta to calculate epsilon
  fscanf(in,"%e",&readp);              // third parameter is circumferential ns radius in km
  rcirc_ns=readp;
  fscanf(in,"%e",&readp);              // fourth parameter is ns mass in Msun
  mns=readp;
  fscanf(in,"%e",&readp);              // fifth parameter is ns spin
  fns=readp;
  fscanf(in,"%e",&readp);              // sixth parameter is theta0
  theta0=readp;
  theta0*=M_PI/180.0;                  // convert theta0 to radians
  fscanf(in,"%e",&readp);              // seventh parameter is latitude of spot
  the_spot=readp;
  the_spot*=M_PI/180.0;                // convert theta0 to radians
  fscanf(in,"%e",&readp);              // eightth parameter is half-opening angle of spot
  rho_spot=readp;
  rho_spot*=M_PI/180.0;                // convert rho_spot to radians
  fscanf(in,"%e",&readp);              // ninth parameter is Tbb on surface
  Tbb=readp;
  
  rns=Circ_to_Coord(rcirc_ns);         // calculate the coordinate ns radius, given the
                                       // value of the circumferential ns radius

  // calculate the dimensionless angular frequency of the neutron star
  // the coefficient below is 2*pi*G*Msolar/c^3 to convert Hz to units of M
  Omegans=2.*M_PI*4.92549e-6*fns*mns;

  if (fclose(in)!=0)                  // close input file
    printf("Error in closing input file <in>\n");

  return;                       // all done
}

/*
void output(double Energy[], double Time[], double dyn_spectra[Nt][Nen])
July 17, 2012 (DP)
   Generates the output file at the and of the calculation and prints it
   on the screen (which can be redirected in a file using the unix
   redirection commands). The output contains three columns: the time
   as a fraction of the rotation phase, the photon energy, the flux at
   infinity. For now, the temperature of the blackbody is 2keV and the distance
   to the source is 10kpc.
{
  int time_index,energy_index;
  double en;           // photon energy in keV
  double flux;         // photon flux
  double bol_flux;     // bolometric photon flux
  double flow, fhigh, color;        // for the X-ray color

  printf("! a=%5.3f e=%5.3f R=%5.3f M=%5.3f\n",a_bh,epsilon,rns,mns);
  printf("! fns=%5.3f theta0=%5.3f the_spot=%5.3f rho_spot=%5.3f Tbb=%5.3f\n",
  	 fns,theta0*180.0/M_PI,the_spot*180.0/M_PI,rho_spot*180.0/M_PI,Tbb);
  for (time_index=1;time_index<=Nt-1;time_index++)
    {
      bol_flux=0.0;        // initialize bolometric flux
      flow=0.0;
      fhigh=0.0;
      for (energy_index=1;energy_index<=Nen;energy_index++)
	{
	  en=Energy[energy_index-1];
	  // this is for photon flux
	  //flux=0.07203*mns*mns/Dist/Dist*dyn_spectra[time_index-1][energy_index-1]/en;
	  flux=0.0720429*mns*mns/Dist/Dist*dyn_spectra[time_index-1][energy_index-1]/en;
	  // add to bolometric flux
	  bol_flux+=flux*(Energy[energy_index]-Energy[energy_index-1]);
	  if (en>=2.*Tbb)
	    {
	      fhigh+=flux*(Energy[energy_index]-Energy[energy_index-1]);
	    }
	  else
	    {
	      flow+=flux*(Energy[energy_index]-Energy[energy_index-1]);
	    }
	  // print flux
	  if (energy_index>1 && energy_index<Nen)
	    printf("%e %e %e\n",Time[time_index-1]*Omegans/2./M_PI,en,flux);
	}
      if (flow>0.0)
	{
	  color=fhigh/flow;
	}
      else
	{
	  color=0.0;
	}
      //printf("%e %e %e \n",Time[time_index-1]*Omegans/2./M_PI,flow,fhigh);
      // printf("%e %e %e %e\n",fns,Time[time_index-1]*Omegans/2./M_PI,bol_flux,color);
    }

  return;
}
*/
