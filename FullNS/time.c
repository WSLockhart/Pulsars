/*
************************************************************
             File: time.c
   Functions that calculate the time-dependent spectrum of
   a rotating hot spot on a neutron star
************************************************************
*/

#include "definitions.h"      // File with useful definitions and headers
#include "global.h"           // Global variable headers

#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>

/*
****************************************************************
   Protypes of external functions used in this file
**************************************************************** */
double ksi(double theta, double phi, double phi0);                    // from file ns.c
double emissivity(double energy, double theta, double phi,
		  double beamangle, double phi0);                     // from file ns.c
int geodesic(double X0, double Y0, double* phi_ns, double* theta_ns,
	     double* beam_ns, double* Einfty_E0, double* Time_delay); // from file image.c
double r_shape(double theta);                                         // from file ns.c
//****************************************************************

void set_time_grid(double Time[])
/* December 4, 2011 (DP)
   Sets a linear grid over times and stores it in the array Time[Nt]. The
   minimum value of the time is 0, the maximum value of the time is
   Time_max times the spin period. This and the number of grid points Nt
   can be found in the global variables (definitions.h).
*/
{
  int index;

  for (index=1;index<=Nt;index++)
    {
      Time[index-1]=Time_max*2.*M_PI/Omegans/(Nt-1.0)*(index-1.0);
    }

  return;
}

void time_spectrum(double Energy[], double Time[])
/* June 18, 2012 (DP)
   Uses an integration of the image plane in order to calculate
   the time-dependent energy spectrum of emission from a rotating hot spot,
   at a location rspot provided as a global variable. The spectrum is
   calculated between time=0 and Time_max*Omega_ns (in definitions.h).
   The array Energy[Nen] has the energy grid (in units of the energy
   of emission) and the array Time[Nt] has the time grid. For each phase point,
   the subroutine first finds the region in the image plane where there is
   significant emission and then sets up a fine grid to perform ray tracing
   for the important parameters. On return, it fills the array
   dyn_spectra[Nt][Nen] with the number of photons at different energy bins,
   as a function of time.
*/
{
  int index;                      // general indexing variable
  double X0,Y0;                   // image-plane coordinates
  int ix,iy;                      // indices of image-plane array that correspond
                                  // to the cell where X0,Y0 resides
  double stepX,stepY;             // steps on the image plane

  double slope;                   // slopes for various linear interpolations
  double aux1,aux2;               // intermediate values for 2D linear interpolation

  double Intensity;               // storage of intensity on image plane

  int result;                     // result flag for geodesic integration

  // array for the spectrum per time bin
  double dyn_spectra[Nen];

  // Arrays with low-resolution values
  double phi_ns, theta_ns, beam_ns, Einfty_E0, Time_delay;

  // Arrays with low-resolution values
  double phi_ns_l[Nx][Ny];
  double theta_ns_l[Nx][Ny];
  double Einfty_E0_l[Nx][Ny];
  double Time_delay_l[Nx][Ny];

  double Xmin_h, Xmax_h;          // min and max values where significant emission exists
  double Ymin_h, Ymax_h;
  double X_h,Y_h;
  double stepX_h,stepY_h;         // steps on fine find grid
  int ix_h, iy_h;                 // indices for fine grid on image plane
  double distance;                // angular distance from a particular point to center of spot

  double theta,phi,beam;          // positions and beaming on neutron star
  double eshift;                  // energy shift along each geodesic
  double delay;                   // time delay along each geodesic

  double time0;                   // zero time defined as the time of phi_0=0 (fiducial plane)
  double phi0;                    // azimuth of fiducial plane
  int energy_index;               // energy index for corresponding bin
  double stepEn,energy;           // step and current value in Energy grid
  int time_index;                 // time index for corresponding bin
  double stepT;                   // step in Time grid
  double integral1,integral2;     // generic variables for integrals
  int didItHit=0;       // flag to see if it hit the neutron star
  int ttyfd = open("/dev/tty", O_RDWR);

#define RADIUS_ACC 1.e-5   // how close to make it to the radius
  double NSradius= r_shape(M_PI/2.);        // equatorial radius

  write(ttyfd, "setting up coarse grid...\n", 27);

  // first initialize the dynamical spectra

  stepX=Xmax/(Nx-1.0);            // step for Nx points along the x-axis
  stepY=Ymax/(Ny-1.0);            // step for Ny points along the y-axis

  // first fill the coarse grid

  for (ix=1;ix<=Nx;ix++)          // for all points on the x-axis
  {
    X0=-0.5*Xmax+stepX*(ix-1.0);// along x-axis, -Xmax/2<=X0<=Xmax/2
    for (iy=1;iy<=Ny;iy++)      // for all points on the y-axis
	  {
	    Y0=-0.5*Ymax+stepY*(iy-1.0);//along y-axis -Ymax/2<=Y0<=Ymax/2

  	  if (!(ix==(Nx-1)/2+1))
  	  {
  	      result=geodesic(X0,Y0,&phi_ns,&theta_ns,&beam_ns,&Einfty_E0,&Time_delay);
  	      // fill the coarse_grid with the result
  	      phi_ns_l[ix-1][iy-1]=phi_ns;
  	      theta_ns_l[ix-1][iy-1]=theta_ns;
  	      Einfty_E0_l[ix-1][iy-1]=Einfty_E0;
  	      Time_delay_l[ix-1][iy-1]=Time_delay;
  	  }
  	  else
  	  {
          phi_ns_l[ix-1][iy-1]=NOT_FOUND;
  	      theta_ns_l[ix-1][iy-1]=NOT_FOUND;
  	      Einfty_E0_l[ix-1][iy-1]=NOT_FOUND;
  	      Time_delay_l[ix-1][iy-1]=NOT_FOUND;
  	      result=0.0;
  	  }
  	  //	  if (result>0.)
  	  // printf("%e %e %e %e %e %e\n",X0,Y0,phi_ns*180./M_PI,180./M_PI*theta_ns,Einfty_E0,
  	  // Time_delay);
	   }
  }

  printf("! finished with coarse grid\n");
  write(ttyfd, "finished with coarse grid!\n", 28);
  // Calculate the time delay at alpha0=0, beta0=0, which is set as the fiducial time
  time0=Time_delay_l[(Nx-1)/2-1][(Ny-1)/2-1];

  // Scan all spin phases
  for (time_index=1;time_index<=Nt-1;time_index++)
  {
    // void* buf = String.Concat("scanning spin phase",time_index.ToString());
    write(ttyfd, "scanning next spin phase\n", 26);
    for (energy_index=1;energy_index<=Nen;energy_index++)
	  {
	     dyn_spectra[energy_index-1]=0.0;
	  }

    // first find limits on image plane where there is significant intensity, i.e., inside spot
    // for all the spin phases. Because the star spins, this is a "zone" of latitudes

    // initialize the limits with very wrong numbers
    Xmin_h=1.e5;
    Xmax_h=-1.e5;
    Ymin_h=1.e5;
    Ymax_h=-1.e5;

    didItHit=0;       // flag to see if it hit the neutron star
    // use the coarse grid to figure out for which range of
    // X0 and Y0 is the location inside the hot spot
    for (ix=1;ix<=Nx;ix++)          // for all points on the x-axis
	  {
        X0=-0.5*Xmax+stepX*(ix-1.0);// along x-axis, -Xmax/2<=X0<=Xmax/2
	      for (iy=1;iy<=Ny;iy++)      // for all points on the y-axis
	      {
	          Y0=-0.5*Ymax+stepY*(iy-1.0);//along y-axis -Ymax/2<=Y0<=Ymax/2
	          // check if it hit the neutron star
	          if (Einfty_E0_l[ix-1][iy-1]>0.0)
		        {
		            theta_ns=theta_ns_l[ix-1][iy-1];
		            phi_ns=phi_ns_l[ix-1][iy-1];

		            // azimuth of center of spot spins around the neutron star
		            phi0=Omegans*(Time[time_index-1]-Time_delay_l[ix-1][iy-1]+time0);

		            // to avoid problems with the pole
		            if (!(X0==0.0 && fabs(Y0-NSradius)>RADIUS_ACC))
		            {
		                distance=ksi(theta_ns,phi_ns,phi0);
		            }
		            else
		            {
		                distance=1.e5;
		            }

		            // check if the distance is within "frac" percent of the size of the
		            // hot spot or of at least 2 degrees=0.175/5 rad
	              //	  if (distance<=frac*rho_spot || distance<=0.175/3.0)
		            if (1<2)    // it always hits
		            {
  		              didItHit=1;
		                if (X0<Xmin_h) Xmin_h=X0;
		                if (X0>Xmax_h) Xmax_h=X0;
		                if (Y0<Ymin_h) Ymin_h=Y0;
		                if (Y0>Ymax_h) Ymax_h=Y0;
		            }
		       }
	     }
	}  //end of x,y grid scan

  // check to see if the spot is visible
  //  if (Xmax_h-Xmin_h>0.0)
  if (didItHit>0)
	{
	  // make sure to enlarge the size by one step in both X and Y
	  // in order to ensure to enclose the spot
	  Xmin_h-=stepX;
	  Xmax_h+=stepX;
	  Ymin_h-=stepY;
	  Ymax_h+=stepY;

	  // Having figured this out, set up the refined grid

	  // set the step size for the fine grid
	  stepX_h=(Xmax_h-Xmin_h)/(Nx_h-1.0);
	  stepY_h=(Ymax_h-Ymin_h)/(Ny_h-1.0);
	  // and perform the ray tracing for the refined grid
	  for (ix_h=1;ix_h<=Nx_h;ix_h++)
    {
	     X_h=Xmin_h+stepX_h*(ix_h-1.0);       // x-position in fine grid

	     for (iy_h=1;iy_h<=Ny_h;iy_h++)
		   {
		       Y_h=Ymin_h+stepY_h*(iy_h-1.0);     // y-position in fine grid

		       // perform a ray tracing
		       result=geodesic(X_h,Y_h,&phi_ns,&theta_ns,&beam_ns,&Einfty_E0,&Time_delay);

		       if (result>0)           // if the new "refined" ray has hit the neutron star
		       {
		           // store for simplicity the redshift
		           eshift=Einfty_E0;
		           // calculate the phase of the spot at the time it emitted
		           phi0=Omegans*(Time[time_index-1]-Time_delay+time0);

		           // add these photons to all the corresponding energies on the grid
		           for (energy_index=2;energy_index<=Nen-1;energy_index++)
               {
                    // if (energy_index==(Nen/2)) //use the middle energy as a typical example //NOTE: PRINT OUT (b): HOTSPOTS
                    {
                      // calculate the boosted, energy-dependent emissivity
			                Intensity=emissivity(Energy[energy_index-1]/eshift,phi_ns,theta_ns,beam_ns,phi0)*(eshift*eshift*eshift);
                      // printf("%e %e %e\n",X_h,Y_h,Intensity); //NOTE: PRINT OUT (b): HOTSPOTS
                    }
  			       // NB, since the limits of the refined grid always have zero
  			       // emissivity, this is equivalent to a trapezoid integration
			              dyn_spectra[energy_index-1]+=Intensity*stepX_h*stepY_h;
			         }
		       }
           //This is important for plotting so that we have a full Nx_h by Ny_h grid
          //  else printf("%e %e %e\n",X_h,Y_h,-0.01);  //NOTE: PRINT OUT (b): HOTSPOTS
           // return -1 for jerry-rigged plot w/ dark sky :)
		    }
	   }  //end of setting up fine grid

   }// if statement for whether the spot is within view

   // calculate the flux and print output
   double flux;
   double totalFlux=0;
   for (energy_index=2;energy_index<=Nen-1;energy_index++)
   {
     flux=0.0720429*mns*mns/Dist/Dist*dyn_spectra[energy_index-1]/Energy[energy_index-1];
     totalFlux+=0.5*flux*(Energy[energy_index]-Energy[energy_index-2]);
    //  printf("%e %e %e\n",Time[time_index-1]*Omegans/2./M_PI,Energy[energy_index-1],flux); //NOTE: PRINT OUT (a): LIGHTCURVES / SPECTRA
     //    if (energy_index==64) printf("%e %e %e\n",Time[time_index-1]*Omegans/2./M_PI,Energy[energy_index-1],flux);
   }
   printf("%e %e %e\n",Time[time_index-1]*Omegans/2./M_PI,0.0,totalFlux);  //NOTE: PRINT OUT (c): PULSEPROFILES


}// end of "time" loop

write(ttyfd, "done!\n", 7);

  //for (time_index=1;time_index<=Nt;time_index++)
  //  {
  //   for (energy_index=1;energy_index<=Nen;energy_index++)
  //	{
    // printf("! %d %e %e\n",time_index,Energy[1],0.0720429*mns*mns/Dist/Dist*dyn_spectra[time_index-1][1]);
  //	}
  //  }

  return;

}  // the end of time
