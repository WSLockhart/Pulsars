/* 
************************************************************
             File: image.c
   Functions that set the grid on the image plane and 
   perform the ray tracing back to the emitting region
************************************************************
*/

#include "definitions.h"      // File with useful definitions and headers
#include "global.h"           // Global variable headers

/* 
****************************************************************
   Protypes of external functions used in this file
**************************************************************** */
extern void metric(double r, double theta, double gmunu[]);
extern void gammas(double r, double theta, double Gamma[4][4][4]);
extern void vel_orb(double r, double velocity[]);
/* **************************************************************** */

double adapt_step(double r[],double k[], double derivs[])
/* April 19, 2010 (DP)
   Calculates the minimum timescale for the change of the variables
   during the integration of the geodesics and sets a timestep equal
   to the minimum timescale divided by the variable "step_factor".
   The minimum timescale is not allowed to be smaller than a
   lightcrossing time for 1M. This allows the integration to go
   through (albeit with a small error) for the small number of rays
   that go close to the pole of the coordinate system.  Checks for
   divisions by zeros etc have been removed in order to increase the
   speed of the calculation. For now, only the change of the
   components of the position vector is considered and not that of the
   momentum vector. At input, r[4] and k[4] are the components of the
   4-vectors for the position and momentum of the photon and derivs[6]
   are the right-hand sides of the equations solved by the Runge-Kutta
   routine.
*/
{
  double mintimescale;

  mintimescale=fabs(r[rco]/derivs[2]);  // try the radius-timescale
                                        // but check if the other two are smaller
  mintimescale=fmin(mintimescale,fabs(r[theco]/derivs[4]));
  mintimescale=fmin(mintimescale,fabs(r[phico]/derivs[1]));

  mintimescale=fmax(mintimescale, 1.0); // mintimescale cannot be smaller than 
                                        // lightcrossing time of 1M
  return mintimescale/step_factor;
}

void RHS(double r[],double k[], double derivs[], double variableadd[], 
	 double factor,double bimpact)
/* April 19, 2010 (DP)
   Returns the right-hand sides of the 6 ODEs to be solved 
   by the Runge-Kutta method. At input, r[4] and k[4] are the
   elements of the 4-vectors for the position and momentum of the
   photon. Because of the particular half steps of the RK4 method,
   the right-hand sides are evaluated (for each quantity) at
   r[]+factor*variableadd[] or k[]+factor*variableadd[]. bimpact is
   the impact parameter of the photon, which could be calculated in 
   principle in this function (slowing down the code). The resulting
   right-hand sides are stored in the array derivs[]. The ordering 
   of the indices for the various arrays is the same as in the function
   trace_rays(). 
*/
{
  double gmunu[5];                       // the metric
  double Gamma[4][4][4];                 // Christoffel symbols
  double radius, theta, phi;             // current coordinates
  double kupt, kupr, kupthe, kupphi;     // components of the photon momentum
  
  // First calculate the "current" position
  radius=r[rco]+factor*variableadd[2];
  theta=r[theco]+factor*variableadd[4];
  phi=r[phico]+factor*variableadd[1];

  // Evaluate the metric elements
  metric(radius,theta,gmunu);
  // Calculate the Christoffel symbols
  gammas(radius,theta,Gamma);

  // Contravariant t-component of the photon 4-momentum from Killing vector
  kupt=(gmunu[3]+bimpact*gmunu[4])/
	(gmunu[3]*gmunu[0]-gmunu[4]*gmunu[4]);
  // Contravariant r-component of the photon 4-momentum
  kupr=k[rco]+factor*variableadd[3];
  // Contravariant theta-component of the photon 4-momentum
  kupthe=k[theco]+factor*variableadd[5];
  // Contravariant phi-component of the photon 4-momentum from Killing vector
  kupphi=(bimpact*gmunu[0]+gmunu[4])/
	(gmunu[3]*gmunu[0]-gmunu[4]*gmunu[4]);

  // This is dt/dlambda=k^t by definition
  derivs[0]=kupt;

  // This is dphi/dlambda=k^phi by definition
  derivs[1]=kupphi;

  // This is dr/dlambda=k^r by definition
  derivs[2]=kupr;
  //  printf("here %e\n",derivs[2]);

  // This is dk^r/dlambda from the geodesic equation
  derivs[3]=-Gamma[rco][tco][tco]*kupt*kupt
    -Gamma[rco][rco][rco]*kupr*kupr
    -2.0*Gamma[rco][theco][rco]*kupthe*kupr
    -Gamma[rco][theco][theco]*kupthe*kupthe
    -2.0*Gamma[rco][phico][tco]*kupphi*kupt
    -Gamma[rco][phico][phico]*kupphi*kupphi;
  //  printf("Here %e\n",derivs[3]);
  // This is dtheta/dlambda=k^theta by definition
  derivs[4]=kupthe;

  // This is dk^theta/dlambda from the geodesic equation
 derivs[5]=-Gamma[theco][tco][tco]*kupt*kupt
    -Gamma[theco][rco][rco]*kupr*kupr
    -2.0*Gamma[theco][theco][rco]*kupthe*kupr
    -Gamma[theco][theco][theco]*kupthe*kupthe
    -2.0*Gamma[theco][phico][tco]*kupphi*kupt
    -Gamma[theco][phico][phico]*kupphi*kupphi;

  return;
}

double trace_rays(double r[],double k[], double bimpact)
/* April 19, 2010 (DP) 

   Traces photon rays from their initial positions to the equatorial
   plane of the metric, on which an accretion disk is assumed to
   exist.  It returns the maximum error in the calculation as inferred
   from the integral of motion k.k=0, where k is the photon
   4-momentum.  The metric of the axisymmetric spacetime is provided
   by the external function
        void metric(double r, double theta, double gmunu[5]) 
   where 'r' and 'theta' are the coordinates, gmunu[0] is the
   tt-component, gmunu[1] is the rr-component, gmunu[2] is the
   theta-theta-component, gmunu[3] is the phi-phi-component, and
   gmunu[4] is the t-phi-component of the metric.  The Christoffel
   symbols are provided by the external function
        void gammas(double r,double theta, double Gamma[4][4][4]) 
   with the indexes for the Gamma array corresponding to the
   parameters {tco,rco,theco,phico}={0,1,2,3}->{t,r,theta,phi}.  The
   initial 4-positions of the photons are provided in the array r_i[4]
   and the initial 4-momenta of the photons are provided in array
   k_i[4]. The impact parameter of the photon bimpact is also
   provided, even though it could be calculated in principle localy.
   This function uses a 4th order Runge-Kutta method with adaptive
   stepsize to solve 6 equations of motions for the photons. For the
   t- and phi-components of the photon position 4-vector it solves the
   1st order equations that come from the Killing vectors. For the r-
   and theta-components of the photon position 4-vector it solves the
   2nd order geodesic equations.  The order of the variables for the
   kRK[] auxiliary arrays is 0:t, 1:phi, 2:r, 3:k^r, 4: theta,
   5:k^theta Upon completion, the function returns the maximum error
   that occurred along the geodesic (as inferred from the deviation of
   the norm of the photon 4-momentum for zero). Also, at exit, the
   arrays r[] and k[] contain the position and momentum 4-vectors of
   the photon.
*/
{
  int index;                             // indexing variable
  double lambda;                         // affine parameter 
  double step;                           // adaptive stepsize
  double kRK1[6],kRK2[6],kRK3[6],kRK4[6];// auxiliary variables for Runge-Kutta
  double null[6]={0.0,0.0,0.0,0.0,0.0,0.0};// null 6-element vector 
  double gmunu[5];                       // the metric
  double integral;                       // integral of motion
  double r_prev[4],k_prev[4];     // previous position and momentum of photon
  double mu1, mu2;                // previous and current cos(theta)
  double Delta;                   // aux variable for quadratic discriminant

  integral=0.0;                          // start with zero error
  lambda=0.0;                            // initialize affine parameter
  do
    {
      // keep a record for current position and momentum of photon   
      for (index=0;index<=3;index++) 
	{
	  r_prev[index]=r[index];
	  k_prev[index]=k[index];
	}
      
      // evaluate the four RK4 steps
      RHS(r,k,kRK1,null,0.0,bimpact);       // evaluate RHS of ODEs for K1
      step=adapt_step(r,k,kRK1);            // adapt the stepsize
      RHS(r,k,kRK2,kRK1,0.5*step,bimpact);  // evaluate RHS of ODEs for K2
      RHS(r,k,kRK3,kRK2,0.5*step,bimpact);  // evaluate RHS of ODEs for K3
      RHS(r,k,kRK4,kRK3,step,bimpact);      // evaluate RHS of ODEs for K4
      
      // This is the t-coordinate
      r[tco]+=step/6.0*(kRK1[0]+2.0*kRK2[0]+2.0*kRK3[0]+kRK4[0]);
      // This is the r-coordinate
      r[rco]+=step/6.0
	*(kRK1[2]+2.0*kRK2[2]+2.0*kRK3[2]+kRK4[2]);
      // This is the theta-coordinate
      r[theco]+=step/6.0
	*(kRK1[4]+2.0*kRK2[4]+2.0*kRK3[4]+kRK4[4]);
      // This is the phi-coordinate
      r[phico]+=step/6.0
	*(kRK1[1]+2.0*kRK2[1]+2.0*kRK3[1]+kRK4[1]);

      // In order to calculate the t- and phi-components of the
      // 4 momentum of the photon, I will use the Killing vectors,
      // for which I will need the metric at the current position
      metric(r[rco],r[theco],gmunu);
      // This is the t-momentum
      k[tco]=(gmunu[3]+bimpact*gmunu[4])/
	(gmunu[3]*gmunu[0]-gmunu[4]*gmunu[4]);
      // This is the r-momentum
      k[rco]+=step/6.0*(kRK1[3]+2.0*kRK2[3]+2.0*kRK3[3]+kRK4[3]);
      // This is the theta-momentum
      k[theco]+=step/6.0*(kRK1[5]+2.0*kRK2[5]+2.0*kRK3[5]+kRK4[5]);
      // This is the phi-momentum
      k[phico]=(bimpact*gmunu[0]+gmunu[4])/
	(gmunu[3]*gmunu[0]-gmunu[4]*gmunu[4]);
      
      // if error is larger than previous error, reset it
      integral=fmax(integral,
		    fabs(fabs(gmunu[0]*k[tco]*k[tco]/(gmunu[1]*k[rco]*k[rco]
			+gmunu[2]*k[theco]*k[theco]+gmunu[3]*k[phico]*k[phico]
			 +2.0*gmunu[4]*k[tco]*k[phico]))-1.0));
		    
      lambda+=step;                   // increase the value of the affine param

      //print statements for debugging
      //      printf("%e %e %e %e %e\n",r[tco],r[rco]*sin(r[theco])*cos(r[phico]),
      //    r[rco]*sin(r[theco])*sin(r[phico]),r[rco]*cos(r[theco]),integral);
      // printf("%e %e %e %e %e\n",r[tco],r[rco],r[theco],r[phico],integral);
    }

  // stop integration if the ray crossed the equatorial plane. also
  // stop the integraton if the ray hit the region where the metric is not
  // well behaved, i.e., r<2.6M, or if it has been scattered out to 10%
  // further than the distance to the observer*/
  while (cos(r[theco])*cos(r_prev[theco])>0.0 // ray did not cross equator
	 && r[rco]>RMIN_TROUBLE       // ray is outside region of trouble 
	 && r[rco]<1.1*Dobs);         // ray got deflected outwards

  // if the ray has actually crossed the equatorial plane
  if (cos(r[theco])*cos(r_prev[theco])<0.0)
    {                                 // interpolate to find crossing point
      mu1=cos(r_prev[theco]);         // previous cos(theta)
      mu2=cos(r[theco]);              // current cos(theta)
      for (index=0;index<=3;index++)  // for all variables
	{    
	  r[index]=(r[index]-r_prev[index])/(mu2-mu1)*(0.0-mu1)
	    +r_prev[index];
	                              // interpolate momentum 4-vector
	  k[index]=(k[index]-k_prev[index])/(mu2-mu1)*(0.0-mu1)
	    +k_prev[index];
	}
      // in order to conserve enery, recalculate k^t from the condition
      // that the norm of the photon 4-momentum is zero.
      // For this, first evaluate the metric elements
      metric(r[rco],r[theco],gmunu);
      // Then the discriminant of the quadratic
      Delta=-gmunu[0]*(gmunu[1]*k[rco]*k[rco]+
		       gmunu[2]*k[theco]*k[theco]+
			   gmunu[3]*k[phico]*k[phico]);
      // and finally the k^t component
      k[tco]=(gmunu[4]*k[phico]+sqrt(Delta))/gmunu[0];
    }
  
  // after all that we can return!
  return integral;
}

void fillimage(double r_equat[Nx][Ny], double Einfty_E0[Nx][Ny])
/* April 19, 2010 (DP)
   Performs the ray tracing for each point on the image plane and returns two 
   arrays: r_equat[Nx][Ny] that contains the radius in the equatorial plane 
   that appears on the (Nx-th,Ny-th) position on the image plane, and 
   Einfty_E0[Nx][Ny] that contains the ratio between the energy observed
   at infinity and the one emitted at the equatorial plane.
*/
{
  int ix, iy;                     // index along the x- and y-axes
  double stepX, stepY;            // grid spacing along the x- and y-axes
  double X0, Y0;                  // initial position of photon on image plane
  double r_i[4];                  // initial 4-position of photons in the
                                  // frame centered on the black hole
  double k_i[4];                  // initial 4-momentum of photons in the
                                  // 3D frame centered on the black hole
  double bimpact;                 // impact parameter of photon
  double costh0,sinth0;           // auxiliary variables for observer
  double costhi,sinthi;           // aux variables for points on image plane
  double cosphi,sinphi;           // aux variables for points on image plane
  double gmunu[5];                // non-zero metric elements
  double Delta;                   // aux variable for quadratic discriminant
  double error;                   // max error in calculation
  double Enum, Eden;              // numerator and denominator of energy expression
  double velocity[4];             // orbital velocity 4-vector

  costh0=cos(theta0);             // useful to calculate once...
  sinth0=sin(theta0);             // ... the cosine and sine of theta0

  stepX=Xmax/(Nx-1.0);            // step for Nx points along the x-axis
  stepY=Ymax/(Ny-1.0);            // step for Ny points along the y-axis

  for (ix=1;ix<=Nx;ix++)          // for all points on the x-axis
    {
      X0=-0.5*Xmax+stepX*(ix-1.0);// along x-axis, -Xmax/2<=X0<=Xmax/2
      for (iy=1;iy<=Ny;iy++)      // for all points on the y-axis
	{
	  Y0=-0.5*Ymax+stepY*(iy-1.0);//along y-axis -Ymax/2<=Y0<=Ymax/2

	  // Calculate the initial position of each photon in the
	  // coordinate system centered on the black hole
	  r_i[tco]=0.0;                         // this is the t-component
	  r_i[rco]=sqrt(Dobs*Dobs+X0*X0+Y0*Y0);             // r-component
	  costhi=(Dobs*costh0+Y0*sinth0) / r_i[1];   // cos of theta-component
	  sinthi=sqrt(1.0-costhi*costhi);
	  r_i[theco]=acos(costhi);                   // theta-component
	  r_i[phico]=atan( X0 / (Dobs*sinth0-Y0*costh0) );    // phi-component
	  cosphi=cos(r_i[phico]);                    // cos of phi-component
	  sinphi=sin(r_i[phico]);                    // sin of phi-component

	  // Calculate the initial 4-momentum of each photon in the
	  // coordinate system centered on the black hole
	  k_i[rco]=-Dobs/r_i[rco];                   // r-component
   	                                             // theta-component
	  k_i[theco]=(costh0-Dobs*(Dobs*costh0+Y0*sinth0)/(r_i[rco]*r_i[rco]))
	    /sqrt(r_i[rco]*r_i[rco]-
		  (Dobs*costh0+Y0*sinth0)*(Dobs*costh0+Y0*sinth0));
	                                             // phi-component 
	  k_i[phico]=X0*sinth0
	    /(X0*X0+
	      (Dobs*sinth0-Y0*costh0)*(Dobs*sinth0-Y0*costh0));
	  // Set the t-component of the initial 4-momentum of each 
	  // photon so that the norm of the 4-momentum vanishes
	  // For this, first evaluate the metric elements
	  metric(r_i[rco],r_i[theco],gmunu);
	  // Then the discriminant of the quadratic
	  Delta=-gmunu[0]*(gmunu[1]*k_i[rco]*k_i[rco]+
		     gmunu[2]*k_i[theco]*k_i[theco]+
			   gmunu[3]*k_i[phico]*k_i[phico]);
	  k_i[tco]=(gmunu[4]*k_i[phico]+sqrt(Delta))/gmunu[0];

	  // Calculate the impact parameter of the photon as b=Lz/E
	  bimpact=(gmunu[3]*k_i[phico]+gmunu[4]*k_i[tco])
	    /(gmunu[0]*k_i[tco]+gmunu[4]*k_i[phico]);

	  // Calculate the numerator for the expression of ratio of
	  // observed to emitted photon energy
	  Enum=-sqrt(-gmunu[tco])*k_i[tco];

	  // trace the rays back to the accretion disk 
	  error=trace_rays(r_i,k_i,bimpact);

	  // if the ray crossed the accretion disk outside the region
	  // of trouble (don't worry about the outer disk radius for now)
	  if (r_i[rco]>=RMIN_TROUBLE)
	    {
	      // Calculate the orbital velocity at the crossing radius
	      vel_orb(r_i[rco],velocity);

	      // Calculate the metric elements at the crossing radius
	      metric(r_i[rco],M_PI/2.0,gmunu);

	      // Calculate the denominator for the expression of ratio of
	      // observed to emitted photon energy
	      Eden=gmunu[tco]*velocity[tco]*k_i[tco]
		+gmunu[phico]*velocity[phico]*k_i[phico]
		+gmunu[4]*(velocity[tco]*k_i[phico]+
			   velocity[phico]*k_i[tco]);

	      // Take the ratio of numerator to denominator and store it
	      Einfty_E0[ix-1][iy-1]=Enum/Eden;
	      // store the radius of intersection with the equatorial plane
	      r_equat[ix-1][iy-1]=r_i[rco];

	    }
	}
    }

  return;
}
