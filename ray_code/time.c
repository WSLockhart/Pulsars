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
// double ksi(double theta, double phi, double phi0);                    // from file ns.c
double emissivity(double energy, double theta, double phi,
		  double beamangle, double phi0);                     // from file ns.c
int geodesic(double X0, double Y0, double* phi_ns, double* theta_ns,
	     double* beam_ns, double* Einfty_E0, double* Time_delay); // from file image.c
double r_shape(double theta);                                         // from file ns.c
double current_density(double theta, double phi);  //from jmath.c
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

/*
    March 2019 (WL)
    Some small changes in the new version of the code which uses more realistic
    hotspots. First, we treat each hemisphere separately (iterating through
    pole=0 and pole=1). This is so that we can see how each emitting region
    contributes to the total pulse profile. We also use the current density
    to set the bounds of the fine grid.
*/

{

  // create output files
  FILE *output_LC;   //lightcurve output
  FILE *output_LC2;
  // FILE *output_HS;   //hotspot output

  char filename[sizeof "../outputs/LC_m.ns_rn.s_fns._tho_ths_T.bb_qns.dat"];
  sprintf(filename, "../outputs/LC_%1.2f_%2.1f_%3.0f_%03d_%03d_%1.2f_%1.1f.dat", \
    mns,rns,fns,(int)(theta0*(180./M_PI)),(int)(the_spot*(180./M_PI)),Tbb,qns);
  output_LC = fopen(filename, "w");
  char filename2[sizeof "../outputs/LC_m.ns_rn.s_fns._tho_ths_T.bb_qns_analytic.dat"];
  sprintf(filename2, "../outputs/LC_%1.1f_%2.1f_%3.0f_%03d_%03d_%1.2f_%1.1f_analytic.dat", \
    mns,rns,fns,(int)(theta0*(180./M_PI)),(int)(the_spot*(180./M_PI)),Tbb,qns);
  output_LC2 = fopen(filename2, "w");
  // char filename2[sizeof filename];
  // sprintf(filename2, "../outputs/HS_80-60-1.5_data/HS_%1.1f_%2.1f_%3.0f_%03d_%03d_%1.2f_%1.1f.dat", \
  //   mns,rns,fns,(int)(theta0*(180./M_PI)),(int)(the_spot*(180./M_PI)),Tbb,qns);
  // output_HS = fopen(filename2, "w");

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
  double phi_ns_l[Nx][Ny];
  double theta_ns_l[Nx][Ny];
  double Einfty_E0_l[Nx][Ny];
  double Time_delay_l[Nx][Ny];

  double Xmin_h, Xmax_h;          // min and max values where significant emission exists
  double Ymin_h, Ymax_h;
  double X_h,Y_h;
  double stepX_h,stepY_h;         // steps on fine find grid
  int ix_h, iy_h;                 // indices for fine grid on image plane
  // double distance;                // angular distance from a particular point to center of spot

  // double theta,phi,beam;          // positions and beaming on neutron star
  double eshift;                  // energy shift along each geodesic
  // double delay;                   // time delay along each geodesic
  double time0;                   // zero time defined as the time of phi_0=0 (fiducial plane)
  double phi0;                    // azimuth of fiducial plane
  int energy_index;               // energy index for corresponding bin
  // double stepEn,energy;           // step and current value in Energy grid
  int time_index;                 // time index for corresponding bin
  double stepT;                   // step in Time grid
  int didItHit=0;                 // flag to see if it hit the neutron star
  int pole;                       //keep track of which hemisphere for tracing spots separately
  double JJ=0;                    //current density, as calculated from Gralla's formulae
  double theta_prime;             // polar angle relative to magnetic axis
  double temperature;             // in case we want to make temperature plots

  #define RADIUS_ACC 1.e-5   // how close to make it to the radius
  double NSradius= r_shape(M_PI/2.);        // equatorial radius

  printf("initializing coarse grid... "); fflush(stdout);

  stepX=Xmax/(Nx-1.0);            // step for Nx points along the x-axis
  stepY=Ymax/(Ny-1.0);            // step for Ny points along the y-axis

  // first fill the coarse grid
  for (ix=1;ix<=Nx;ix++)          // for all points on the x-axis
  {
      X0=-0.5*Xmax+stepX*(ix-1.0);// along x-axis, -Xmax/2<=X0<=Xmax/2
      for (iy=1;iy<=Ny;iy++)      // for all points on the y-axis
  	  {
    	    Y0=-0.5*Ymax+stepY*(iy-1.0);//along y-axis -Ymax/2<=Y0<=Ymax/2
      	  if (!(ix==(Nx-1)/2+1)){
      	      result=geodesic(X0,Y0,&phi_ns,&theta_ns,&beam_ns,&Einfty_E0,&Time_delay);
      	      // fill the coarse_grid with the result
      	      phi_ns_l[ix-1][iy-1]=phi_ns;
      	      theta_ns_l[ix-1][iy-1]=theta_ns;
      	      Einfty_E0_l[ix-1][iy-1]=Einfty_E0;
      	      Time_delay_l[ix-1][iy-1]=Time_delay;
      	  }
      	  else{
              phi_ns_l[ix-1][iy-1]=NOT_FOUND;
      	      theta_ns_l[ix-1][iy-1]=NOT_FOUND;
      	      Einfty_E0_l[ix-1][iy-1]=NOT_FOUND;
      	      Time_delay_l[ix-1][iy-1]=NOT_FOUND;
      	      result=0.0;
      	  }
  	  }
  }
  // Calculate the time delay at alpha0=0, beta0=0, which is set as the fiducial time
  time0=Time_delay_l[(Nx-1)/2-1][(Ny-1)/2-1];

  printf("coarse grid complete."); fflush(stdout);

/********************************************************
                       MAIN LOOP
 ********************************************************/

  // Scan all spin phases
  for (time_index=1;time_index<=Nt-1;time_index++)
  // time_index=1;
  {
    printf("\nspin phase %d/%d;",time_index,Nt-1);
    for (energy_index=1;energy_index<=Nen;energy_index++){
	     dyn_spectra[energy_index-1]=0.0;
	  }

  /* We treat the north and south hemispheres (w.r.t. the magnetic axis) separately
  so that we can see how each contributes to the total light curve. */
  // for (pole=0;pole<=1;pole++)
  pole=0;
  {
    // first find a bounding rectangle on the image plane where there is
    //significant intensity, i.e., inside the spot

    // initialize the limits with very wrong numbers
    Xmin_h=1.e5;
    Xmax_h=-1.e5;
    Ymin_h=1.e5;
    Ymax_h=-1.e5;
    didItHit=0;  // flag to see if it hit the neutron star

    // use the coarse grid to figure out for which range of
    // X0 and Y0 is the location inside the hot spot
    for (ix=1;ix<=Nx;ix++)          // for all points on the x-axis
	  {
        X0=-0.5*Xmax+stepX*(ix-1.0);  // along x-axis, -Xmax/2<=X0<=Xmax/2
	      for (iy=1;iy<=Ny;iy++)        // for all points on the y-axis
	      {
	          Y0=-0.5*Ymax+stepY*(iy-1.0);//along y-axis -Ymax/2<=Y0<=Ymax/2
	          // check if it hit the neutron star
	          if (Einfty_E0_l[ix-1][iy-1]>0.0)
		        {
		            theta_ns=theta_ns_l[ix-1][iy-1];
		            phi_ns=phi_ns_l[ix-1][iy-1];
		            // azimuth of center of spot spins around the neutron star
		            phi0=Omegans*(Time[time_index-1]-Time_delay_l[ix-1][iy-1]+time0);

              // check if we're in the right hemisphere
              theta_prime = acos(cos(-the_spot)*cos(theta_ns) -
                cos(phi_ns-phi0)*sin(-the_spot)*sin(theta_ns));
              if ( (2*pole-1)*(theta_prime-M_PI/2.) >= 0 ) // +/-(theta-pi/2)<0 selects N/S hemisphere
              {
		            // // to avoid problems with the pole
		            // if (!(X0==0.0 && fabs(Y0-NSradius)>RADIUS_ACC))

                //use J^2 to determine the bounds of the fine grid
                //NOTE: must pass (phi_ns-phi0) to account for the star's rotation!
                JJ = current_density(theta_ns,phi_ns-phi0);
                if (JJ!=0)
		            {
  		              didItHit=1;
		                if (X0<Xmin_h) Xmin_h=X0;
		                if (X0>Xmax_h) Xmax_h=X0;
		                if (Y0<Ymin_h) Ymin_h=Y0;
		                if (Y0>Ymax_h) Ymax_h=Y0;
		            }
              }
		       }  //end if it hit the star
	     }
	}  //end of x,y grid scan

  // check to see if the spot is visible
  //  if (Xmax_h-Xmin_h>0.0)
  if (didItHit>0)
	{
    // printf("fine grid set (%.2fx%.2f). ",(Xmax_h-Xmin_h)/Xmax,(Ymax_h-Ymin_h)/Ymax);
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
    printf(" ray-tracing..."); fflush(stdout);

	  for (ix_h=1;ix_h<=Nx_h;ix_h++)
    {
	     X_h=Xmin_h+stepX_h*(ix_h-1.0);       // x-position in fine grid
	     for (iy_h=1;iy_h<=Ny_h;iy_h++)
		   {
		       Y_h=Ymin_h+stepY_h*(iy_h-1.0);     // y-position in fine grid

		       // perform a ray tracing
		       result=geodesic(X_h,Y_h,&phi_ns,&theta_ns,&beam_ns,&Einfty_E0,&Time_delay);
		       if (result>0)  // if the new "refined" ray has hit the neutron star
		       {
		           // store for simplicity the redshift
		           eshift=Einfty_E0;
		           // calculate the phase of the spot at the time it emitted
		           phi0=Omegans*(Time[time_index-1]-Time_delay+time0);

               //NB: NEED to check here too so we don't include the other pole in our fine grid!
               theta_prime = acos(cos(-the_spot)*cos(theta_ns) -
                 cos(phi_ns-phi0)*sin(-the_spot)*sin(theta_ns));
               if ( (2*pole-1)*(theta_prime-M_PI/2.) >= 0)
               {
	               // add these photons to all the corresponding energies on the grid
  		           for (energy_index=2;energy_index<=Nen-1;energy_index++)
                 {
                    // calculate the boosted, energy-dependent emissivity
  	                Intensity=emissivity(Energy[energy_index-1]/eshift,phi_ns,theta_ns,beam_ns,phi0)*(eshift*eshift*eshift);
                    // NB, since the limits of the refined grid always have zero
                    // emissivity, this is equivalent to a trapezoid integration
                    dyn_spectra[energy_index-1]+=Intensity*stepX_h*stepY_h;
                 }

                //  fprintf(output_HS,"%e %e %e\n",X_h,Y_h,current_density(theta_ns,phi_ns-phi0)); //PRINT OUT: HOTSPOTS

                //  //print temperature
                //  JJ = current_density(theta_ns,phi_ns-phi0);
                //  if(JJ>0){
                //    temperature = Tbb*8.09*pow(JJ,(1/8.));
                //    fprintf(output_HS,"%e %e %e\n",X_h,Y_h,temperature); //PRINT OUT: HOTSPOTS
                //  } else if(JJ<0){  //jenky bs to plot the negative J^2 region without including the whole star
                //    // temperature = -Tbb*8.09*pow(-JJ,(1/8.));
                //    temperature=0;
                //    fprintf(output_HS,"%e %e %e\n",X_h,Y_h,temperature); //PRINT OUT: HOTSPOTS
                //  }

			         }

		       } // end if hit
		    } //end fine grid y
     }  //end fine grid x

   } // if statement for whether the spot is within view

 } //end hemisphere loop


/****************************************************/

   // calculate the total flux
   double flux;
   double totalFlux=0.0;
   for (energy_index=2;energy_index<=Nen-1;energy_index++)
   {
     flux=0.0720429*mns*mns/Dist/Dist*dyn_spectra[energy_index-1]/Energy[energy_index-1];
     totalFlux+=0.5*flux*(Energy[energy_index]-Energy[energy_index-2]);
    //  fprintf(output_LC,"%e %e %e\n",Time[time_index-1]*Omegans/2./M_PI,Energy[energy_index-1],flux); //PRINT OUT: LIGHTCURVES / SPECTRA
   }
   fprintf(output_LC,"%e %e %e\n",Time[time_index-1]*Omegans/2./M_PI,0.0,totalFlux);  //PRINT OUT: PULSEPROFILES


   //NB: compare with simplified analytic model:

   double zeta=the_spot;
   double phase = Time[time_index-1]*Omegans;

   double h = 0; //beaming factor, I = I0 (1 + h*cos(Theta)). h=0 means isotropic beaming
   double u = 2*(1.47662)*mns/rns; // numerical value = (G*Msolar/c^2)/1000
   //to make u the dimensionless compactness. //NB: double-check!
   double q = u + (1-u)*cos(zeta)*cos(theta0);
   double v = (1-u)*sin(zeta)*sin(theta0);
   double A = 1; //Normalization A = (I_0*G*dS/D^2) sets an order of magnitude
   //for how many photons we expect in each time bin
   double Q = (q + h*(q*q + v*v/2));   // Q=q  for h=0
   double V1 = (1 + 2*h*q)*v;          // V1=v for h=0
   double V2 = h*v*v/2;                // V2=0 for h=0

   double epsilon = 0.07; //NB: 2PiR/cP in the paper is just epsilon,
   //the dimensionless surface velocity
   double onoff=1; //boolean toggle for velocity-dependent factors phi1 & phi2
   double Gamma=0; //photon spectral index. what is a reasonable guess?
   double phi1 = onoff * atan( (((3+Gamma)*q + (4+Gamma)*h*(q*q + v*v/4))/u/(1+2*h*q))*epsilon*sin(zeta)*sin(theta0) );
   double phi2 = onoff * atan( (((4+Gamma)*(1+2*h*q)-1)/h/u)*epsilon*sin(zeta)*sin(theta0) );

   // consider just the first two harmonics
   double aflux = A * (Q + V1*cos(phase+phi1) + V2*cos(2*phase+phi2)); //flux F(t)
   // printf("result=%e\n",result);

   if(aflux<0) aflux=0; //NB: negative flux corresponds to
   //when the spot is hidden behind the star (I THINK!) *double-check this w/ D&F
   fprintf(output_LC2,"%e %e %e\n",phase/2./M_PI,0.0,aflux);


}// end of phase loop

printf("\nsimulation complete!\n");
fclose(output_LC);
fclose(output_LC2);
// fclose(output_HS);


return;

}  // the end of time
