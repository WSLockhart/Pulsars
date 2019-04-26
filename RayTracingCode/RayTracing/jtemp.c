#include <stdio.h>
#include <math.h>
#include <stdlib.h>
// #include <fcntl.h>
#include <unistd.h>
// #include "definitions.h"      // File with useful definitions and headers
// #include "global.h"           // Global variable headers


double mns;
double rns;
double fns;
double qns;
double zeta;

double current_density(double theta, double phi)
/* July 11, 2018 (WL)
   Calculate the current density J(theta,phi) on the surface of the star.
   This is based on Sam Gralla's analytic results, coded originally in
   Mathematica and converted here.
*/

/* Aug 28, 2018 (WL)
  In gravitational (c=G=1) units, mass is scaled out and everything is measured in just time and length.
  We should be working only in the dimensionless parameter C=M/R, but because the formulae for J are
  written with r-derivatives, this will be a messy conversion. For now we will set M=1 and scale the
  radius r and frequency Omega by mns.

*/

{

  double iota = -zeta;  //inclination of the magnetic dipole axis.
                            //negative because of how the rotation is done
  // double M = 1;  //(see note above)
  // double R = rns/(1.47662*mns); // divide by (G*Msolar/c^2) to make it dimensionless,
  //                               // multiply by (1000km/1m), and divide by mns
  // double Omega = Omegans;       // defined in io.c
  // double M = mns; double R = rns; double Omega = 2*M_PI*fns;

  double M = mns;  //(see note above)
  double R = rns/(1.47662); // divide by (G*Msolar/c^2) to make it dimensionless,
                                // multiply by (1000km/1m), and divide by mns
  double Omega = 2.*M_PI*4.92549e-6*fns;

  double q = qns;  //ratio of quadrupole moment to dipole moment
  double II = 2./5;  //dimensionless moment of inertia (MOI)

  // printf("iota=%f,M=%f,R=%f,Omega=%f,q=%f\n", iota,M,R,Omega,q);
  // double M=0.25; double R=1; double Omega=0.1; double q=0; //test values [C=0.5,epsilon=0.1]
  // double M = 1; double r = 5.8; double Omega = 0.012; //test
  //debug check for dimensionless consistency:
  // M *= (0.5); R *= (0.5); Omega *= 2;

  // optimize code by computing some trig functions beforehand
  double f = 1-(2*M/R);

  double sin_theta = sin(theta);
  double cos_theta = cos(theta);
  double sin_phi = sin(phi);
  double cos_phi = cos(phi);
  double sin_iota = sin(iota);
  double cos_iota = cos(iota);
  double log_f = log(f);


  double MOI = II*M*R*R;
  double OmegaZ = Omega*(2*MOI)/(R*R*R);  //"frame-dragging frequency"

  double mu = 1;
  double mu2 = (9*mu*q*M*M*pow(R,-1)*(3 + 2*log_f - 4*f + f*f)*
     pow(17 - 9*f + 6*log_f*(1 + 3*f) - 9*f*f + pow(f,3),-1))/5.;

  double theta_prime = acos(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta);
  double phi_prime = atan2((sin_theta*sin_phi), (cos_theta*sin_iota + cos_iota*cos_phi*sin_theta));

  double R1 = -(1.5)*((3 + 2*log_f - 4*f + f*f)/(pow(1 - f,3)))/R;
  double alpha = mu*R1*(1+q*cos(theta_prime))*sin(theta_prime)*sin(theta_prime);  //this is alpha ON THE SURFACE (r=R)
  double beta = phi_prime;
  double alpha0 = sqrt(1.5)*mu*Omega*(1 + sin_iota*sin_iota/5.);

  double Jsquared;

  int northcap;  // correction to the paper for southern cap
  if (theta_prime<(M_PI/2)) northcap = 1; else northcap = -1;

  //outside the spacelike region, no 3-current is allowed
  //after adding the quadrupole term, we must also restrict alpha to positive values
  if (alpha>=alpha0 || alpha<=0){
    Jsquared = 0.0;
    // printf("no current detected at (%.1f,%.1f).\n",theta*(180.0/M_PI),phi*(180.0/M_PI));
  }
  else
  {

    double Jt; double Jr; double Jtheta; double Jphi;

    double Lambda = -2*Omega*
      (j0(2*asin(sqrt(alpha/alpha0)))*cos_iota -
      northcap*j1(2*asin(sqrt(alpha/alpha0)))*cos(beta)*sin_iota);  //j0 and j1 are the bessel functions


// massive formulae for the four-current in a rotated coordinate system! //


    Jt = (Omega - OmegaZ)*pow(R,-2)*pow(f,-0.5)*(-((cos_theta*
           pow(pow(sin_phi,2)*pow(sin_theta,2) + pow(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta,2),-1)*sin_phi*
           (cos_theta*sin_iota + cos_iota*cos_phi*sin_theta) -
          pow(pow(sin_phi,2)*pow(sin_theta,2) + pow(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta,2),-1)*sin_phi*sin_theta*
           (cos_iota*cos_phi*cos_theta - sin_iota*sin_theta))*
        ((-3*mu*q*cos_theta*pow(M,-3)*R*R*(3 + 2*log_f - 4*f + f*f)*
             (1 - pow(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta,2))*sin_iota*sin_phi)/16. +
          (3*mu*q*cos_theta*pow(M,-3)*R*R*(3 + 2*log_f - 4*f + f*f)*
             pow(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta,2)*sin_iota*sin_phi)/8. +
          (3*mu*pow(M,-3)*R*R*(3 + 2*log_f - 4*f + f*f)*sin_iota*sin_phi*
             sin_theta*(-(cos_phi*cos_theta*sin_iota) - cos_iota*sin_theta))/8. +
          (3*mu*cos_theta*pow(M,-3)*R*R*(3 + 2*log_f - 4*f + f*f)*sin_iota*
             sin_phi*(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta))/8. +
          (9*mu*q*pow(M,-3)*R*R*(3 + 2*log_f - 4*f + f*f)*sin_iota*sin_phi*
             sin_theta*(-(cos_phi*cos_theta*sin_iota) - cos_iota*sin_theta)*(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta))/8.))\
      + ((-3*mu*q*pow(M,-3)*R*R*(3 + 2*log_f - 4*f + f*f)*
           (1 - pow(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta,2))*(-(cos_phi*cos_theta*sin_iota) - cos_iota*sin_theta))/16. +
        (3*mu*q*pow(M,-3)*R*R*(3 + 2*log_f - 4*f + f*f)*
           pow(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta,2)*(-(cos_phi*cos_theta*sin_iota) - cos_iota*sin_theta))/8. +
        (3*mu*pow(M,-3)*R*R*(3 + 2*log_f - 4*f + f*f)*
           (-(cos_phi*cos_theta*sin_iota) - cos_iota*sin_theta)*(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta))/8.)*
      (2*cos_iota*cos_theta*pow(pow(sin_phi,2)*pow(sin_theta,2) + pow(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta,2),-1)*
         pow(sin_phi,2)*sin_theta + cos_phi*cos_theta*
         pow(pow(sin_phi,2)*pow(sin_theta,2) + pow(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta,2),-1)*
         (cos_theta*sin_iota + cos_iota*cos_phi*sin_theta) +
        cos_phi*pow(pow(sin_phi,2)*pow(sin_theta,2) + pow(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta,2),-1)*sin_theta*
         (cos_iota*cos_phi*cos_theta - sin_iota*sin_theta) -
        cos_iota*pow(pow(sin_phi,2)*pow(sin_theta,2) + pow(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta,2),-2)*pow(sin_phi,2)*
         pow(sin_theta,2)*(2*cos_theta*pow(sin_phi,2)*sin_theta +
           2*(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta)*(cos_iota*cos_phi*cos_theta - sin_iota*sin_theta)) -
        cos_phi*pow(pow(sin_phi,2)*pow(sin_theta,2) + pow(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta,2),-2)*sin_theta*
         (cos_theta*sin_iota + cos_iota*cos_phi*sin_theta)*
         (2*cos_theta*pow(sin_phi,2)*sin_theta + 2*(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta)*
            (cos_iota*cos_phi*cos_theta - sin_iota*sin_theta))) -
     (1/sin_theta)*((-3*mu*q*pow(M,-3)*R*R*(3 + 2*log_f - 4*f + f*f)*
           (1 - pow(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta,2))*sin_iota*sin_phi*sin_theta)/16. +
        (3*mu*q*pow(M,-3)*R*R*(3 + 2*log_f - 4*f + f*f)*
           pow(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta,2)*sin_iota*sin_phi*sin_theta)/8. +
        (3*mu*pow(M,-3)*R*R*(3 + 2*log_f - 4*f + f*f)*sin_iota*sin_phi*
           sin_theta*(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta))/8.)*
      (cos_theta*(cos_theta*pow(pow(sin_phi,2)*pow(sin_theta,2) + pow(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta,2),-1)*sin_phi*
            (cos_theta*sin_iota + cos_iota*cos_phi*sin_theta) -
           pow(pow(sin_phi,2)*pow(sin_theta,2) + pow(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta,2),-1)*sin_phi*sin_theta*
            (cos_iota*cos_phi*cos_theta - sin_iota*sin_theta)) +
        sin_theta*(-(pow(pow(sin_phi,2)*pow(sin_theta,2) + pow(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta,2),-1)*sin_phi*
              sin_theta*(-(cos_theta*sin_iota) - cos_iota*cos_phi*sin_theta)) -
           pow(pow(sin_phi,2)*pow(sin_theta,2) + pow(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta,2),-1)*sin_phi*sin_theta*
            (cos_theta*sin_iota + cos_iota*cos_phi*sin_theta) -
           cos_theta*pow(pow(sin_phi,2)*pow(sin_theta,2) + pow(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta,2),-2)*sin_phi*
            (cos_theta*sin_iota + cos_iota*cos_phi*sin_theta)*
            (2*cos_theta*pow(sin_phi,2)*sin_theta + 2*(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta)*
               (cos_iota*cos_phi*cos_theta - sin_iota*sin_theta)) +
           pow(pow(sin_phi,2)*pow(sin_theta,2) + pow(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta,2),-2)*sin_phi*sin_theta*
            (cos_iota*cos_phi*cos_theta - sin_iota*sin_theta)*
            (2*cos_theta*pow(sin_phi,2)*sin_theta + 2*(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta)*
               (cos_iota*cos_phi*cos_theta - sin_iota*sin_theta)))) +
     (cos_iota*pow(pow(sin_phi,2)*pow(sin_theta,2) + pow(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta,2),-1)*pow(sin_phi,2)*
         pow(sin_theta,2) + cos_phi*pow(pow(sin_phi,2)*pow(sin_theta,2) + pow(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta,2),-1)*
         sin_theta*(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta))*
      (f*(R*R*((-3*mu*pow(M,-3)*R*R*
                 (8*pow(M,2)*pow(R,-4) + 16*M*pow(R,-3) - 8*M*pow(R,-3)*f - 8*pow(M,2)*pow(R,-4)*pow(f,-2) -
                   8*M*pow(R,-3)*pow(f,-1))*(1 - pow(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta,2)))/16. -
              (3*mu*R*pow(M,-3)*(-8*M*pow(R,-2) + 4*M*pow(R,-2)*f + 4*M*pow(R,-2)*pow(f,-1))*
                 (1 - pow(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta,2)))/4. -
              (3*mu*pow(M,-3)*(3 + 2*log_f - 4*f + f*f)*
                 (1 - pow(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta,2)))/8. -
              (9*mu*q*pow(M,-3)*(3 + 2*log_f - 4*f + f*f)*
                 (1 - pow(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta,2))*(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta))/8. -
              (3*mu*q*pow(M,-3)*R*R*(3 + 2*log_f - 4*f + f*f)*
                 (-72*pow(M,2)*pow(R,-4) + 36*M*pow(R,-3) - 72*M*log_f*pow(R,-3) + 24*pow(M,2)*pow(R,-4)*f +
                   72*M*pow(R,-3)*f - 24*pow(M,2)*pow(R,-4)*(1 + 3*f)*pow(f,-2) +
                   144*pow(M,2)*pow(R,-4)*pow(f,-1) - 24*M*pow(R,-3)*(1 + 3*f)*pow(f,-1) -
                   12*M*pow(R,-3)*f*f)*
                 pow(17 - 9*f + 6*log_f*(1 + 3*f) - 9*f*f +
                   pow(f,3),-1)*(1 - pow(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta,2))*
                 (cos_iota*cos_theta - cos_phi*sin_iota*sin_theta))/16. -
              (9*mu*q*R*pow(M,-3)*(3 + 2*log_f - 4*f + f*f)*
                 (-18*M*pow(R,-2) + 36*M*log_f*pow(R,-2) - 36*M*pow(R,-2)*f +
                   12*M*pow(R,-2)*(1 + 3*f)*pow(f,-1) + 6*M*pow(R,-2)*f*f)*
                 pow(17 - 9*f + 6*log_f*(1 + 3*f) - 9*f*f +
                   pow(f,3),-1)*(1 - pow(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta,2))*
                 (cos_iota*cos_theta - cos_phi*sin_iota*sin_theta))/8.) +
           2*R*((-3*mu*pow(M,-3)*R*R*(-8*M*pow(R,-2) + 4*M*pow(R,-2)*f + 4*M*pow(R,-2)*pow(f,-1))*
                 (1 - pow(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta,2)))/16. -
              (3*mu*R*pow(M,-3)*(3 + 2*log_f - 4*f + f*f)*
                 (1 - pow(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta,2)))/8. -
              (9*mu*q*R*pow(M,-3)*(3 + 2*log_f - 4*f + f*f)*
                 (1 - pow(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta,2))*(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta))/16.\
               - (3*mu*q*pow(M,-3)*R*R*(3 + 2*log_f - 4*f + f*f)*
                 (-18*M*pow(R,-2) + 36*M*log_f*pow(R,-2) - 36*M*pow(R,-2)*f +
                   12*M*pow(R,-2)*(1 + 3*f)*pow(f,-1) + 6*M*pow(R,-2)*f*f)*
                 pow(17 - 9*f + 6*log_f*(1 + 3*f) - 9*f*f +
                   pow(f,3),-1)*(1 - pow(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta,2))*
                 (cos_iota*cos_theta - cos_phi*sin_iota*sin_theta))/16.)) +
        (1/sin_theta)*(cos_theta*((-3*mu*q*pow(M,-3)*R*R*(3 + 2*log_f - 4*f + f*f)*
                 (1 - pow(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta,2))*(-(cos_phi*cos_theta*sin_iota) - cos_iota*sin_theta))/
               16. + (3*mu*q*pow(M,-3)*R*R*(3 + 2*log_f - 4*f + f*f)*
                 pow(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta,2)*(-(cos_phi*cos_theta*sin_iota) - cos_iota*sin_theta))/8. +
              (3*mu*pow(M,-3)*R*R*(3 + 2*log_f - 4*f + f*f)*
                 (-(cos_phi*cos_theta*sin_iota) - cos_iota*sin_theta)*(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta))/8.) +
           sin_theta*((3*mu*pow(M,-3)*R*R*(3 + 2*log_f - 4*f + f*f)*
                 pow(-(cos_phi*cos_theta*sin_iota) - cos_iota*sin_theta,2))/8. +
              (9*mu*q*pow(M,-3)*R*R*(3 + 2*log_f - 4*f + f*f)*
                 pow(-(cos_phi*cos_theta*sin_iota) - cos_iota*sin_theta,2)*(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta))/8. -
              (3*mu*q*pow(M,-3)*R*R*(3 + 2*log_f - 4*f + f*f)*
                 (1 - pow(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta,2))*(-(cos_iota*cos_theta) + cos_phi*sin_iota*sin_theta))/
               16. + (3*mu*q*pow(M,-3)*R*R*(3 + 2*log_f - 4*f + f*f)*
                 pow(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta,2)*(-(cos_iota*cos_theta) + cos_phi*sin_iota*sin_theta))/8. +
              (3*mu*pow(M,-3)*R*R*(3 + 2*log_f - 4*f + f*f)*
                 (cos_iota*cos_theta - cos_phi*sin_iota*sin_theta)*(-(cos_iota*cos_theta) + cos_phi*sin_iota*sin_theta))/8.))));


    Jr = Lambda*(1/sin_theta)*pow(R,-1)*pow(f*R*R,-0.5)*
   (-((cos_theta*pow(pow(sin_phi,2)*pow(sin_theta,2) + pow(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta,2),-1)*sin_phi*
           (cos_theta*sin_iota + cos_iota*cos_phi*sin_theta) -
          pow(pow(sin_phi,2)*pow(sin_theta,2) + pow(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta,2),-1)*sin_phi*sin_theta*
           (cos_iota*cos_phi*cos_theta - sin_iota*sin_theta))*
        ((-3*mu*q*pow(M,-3)*R*R*(3 + 2*log_f - 4*f + f*f)*
             (1 - pow(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta,2))*sin_iota*sin_phi*sin_theta)/16. +
          (3*mu*q*pow(M,-3)*R*R*(3 + 2*log_f - 4*f + f*f)*
             pow(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta,2)*sin_iota*sin_phi*sin_theta)/8. +
          (3*mu*pow(M,-3)*R*R*(3 + 2*log_f - 4*f + f*f)*sin_iota*sin_phi*
             sin_theta*(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta))/8.)) +
     (cos_iota*pow(pow(sin_phi,2)*pow(sin_theta,2) + pow(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta,2),-1)*pow(sin_phi,2)*
         pow(sin_theta,2) + cos_phi*pow(pow(sin_phi,2)*pow(sin_theta,2) + pow(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta,2),-1)*
         sin_theta*(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta))*
      ((-3*mu*q*pow(M,-3)*R*R*(3 + 2*log_f - 4*f + f*f)*
           (1 - pow(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta,2))*(-(cos_phi*cos_theta*sin_iota) - cos_iota*sin_theta))/16. +
        (3*mu*q*pow(M,-3)*R*R*(3 + 2*log_f - 4*f + f*f)*
           pow(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta,2)*(-(cos_phi*cos_theta*sin_iota) - cos_iota*sin_theta))/8. +
        (3*mu*pow(M,-3)*R*R*(3 + 2*log_f - 4*f + f*f)*
           (-(cos_phi*cos_theta*sin_iota) - cos_iota*sin_theta)*(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta))/8.));


    Jtheta = Lambda*(1/sin_theta)*pow(R,-1)*(cos_iota*pow(pow(sin_phi,2)*pow(sin_theta,2) + pow(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta,2),-1)*
      pow(sin_phi,2)*pow(sin_theta,2) + cos_phi*pow(pow(sin_phi,2)*pow(sin_theta,2) +
        pow(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta,2),-1)*sin_theta*(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta))*
   ((3*mu*pow(M,-3)*R*R*(-8*M*pow(R,-2) + 4*M*pow(R,-2)*f + 4*M*pow(R,-2)*pow(f,-1))*
        (1 - pow(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta,2)))/16. +
     (3*mu*R*pow(M,-3)*(3 + 2*log_f - 4*f + f*f)*
        (1 - pow(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta,2)))/8. +
     (9*mu*q*R*pow(M,-3)*(3 + 2*log_f - 4*f + f*f)*
        (1 - pow(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta,2))*(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta))/16. +
     (3*mu*q*pow(M,-3)*R*R*(3 + 2*log_f - 4*f + f*f)*
        (-18*M*pow(R,-2) + 36*M*log_f*pow(R,-2) - 36*M*pow(R,-2)*f +
          12*M*pow(R,-2)*(1 + 3*f)*pow(f,-1) + 6*M*pow(R,-2)*f*f)*
        pow(17 - 9*f + 6*log_f*(1 + 3*f) - 9*f*f +
          pow(f,3),-1)*(1 - pow(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta,2))*
        (cos_iota*cos_theta - cos_phi*sin_iota*sin_theta))/16.);


    Jphi = Lambda*pow(R,-1)*(cos_theta*pow(pow(sin_phi,2)*pow(sin_theta,2) + pow(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta,2),-1)*sin_phi*
      (cos_theta*sin_iota + cos_iota*cos_phi*sin_theta) -
     pow(pow(sin_phi,2)*pow(sin_theta,2) + pow(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta,2),-1)*sin_phi*sin_theta*
      (cos_iota*cos_phi*cos_theta - sin_iota*sin_theta))*
   ((-3*mu*pow(M,-3)*R*R*(-8*M*pow(R,-2) + 4*M*pow(R,-2)*f + 4*M*pow(R,-2)*pow(f,-1))*
        (1 - pow(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta,2)))/16. -
     (3*mu*R*pow(M,-3)*(3 + 2*log_f - 4*f + f*f)*
        (1 - pow(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta,2)))/8. -
     (9*mu*q*R*pow(M,-3)*(3 + 2*log_f - 4*f + f*f)*
        (1 - pow(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta,2))*(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta))/16. -
     (3*mu*q*pow(M,-3)*R*R*(3 + 2*log_f - 4*f + f*f)*
        (-18*M*pow(R,-2) + 36*M*log_f*pow(R,-2) - 36*M*pow(R,-2)*f +
          12*M*pow(R,-2)*(1 + 3*f)*pow(f,-1) + 6*M*pow(R,-2)*f*f)*
        pow(17 - 9*f + 6*log_f*(1 + 3*f) - 9*f*f +
          pow(f,3),-1)*(1 - pow(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta,2))*
        (cos_iota*cos_theta - cos_phi*sin_iota*sin_theta))/16.);


    double JmuJmu = -Jt*Jt + Jr*Jr + Jtheta*Jtheta + Jphi*Jphi;

    // Jsquared = JmuJmu;

    //spot type (D): only 3-current hotspot
    if(JmuJmu>0){
      Jsquared = Jr*Jr + Jtheta*Jtheta + Jphi*Jphi;
    } else{
      Jsquared=0;
    }

  }
  return Jsquared;
}


double temp_contours(double q, double inc_angle)
{
  //set global variables
  mns = 1.5;
  rns = 11;
  fns = 300;
  qns = q;
  zeta = inc_angle*(M_PI/180);  //convert from degrees to radians
  //

  int N_th=1024; //number of theta points in mesh
  int N_ph=2048; //number of phi points in mesh
  int th; int ph;
  double stepTH=M_PI/N_th;
  double stepPH=2*M_PI/N_ph;  //step size for iterating over the mesh
  double north_area=0; //count how many points have current in each hemisphere
  double south_area=0; //multiplied by the area element sin_theta
  double north_total=0;
  double south_total=0; //total temperature count

  double theta;
  double phi;
  double JJ;
  double temperature;
  double theta_prime;

  /* Iterate over the whole sphere. Every time we find a positive current,
     identify which hemisphere it's in and add it to the count.
     Then at the end we'll report the average temperature of each hemisphere.
  */
  for (th=0;th<=N_th;th++){
    theta=stepTH*th;       // theta-position in fine grid
    for (ph=0;ph<=N_ph;ph++){
      phi=stepPH*ph;       // phi-position in fine grid

      JJ = current_density(theta,phi);
      if (JJ>0)  //must be > not >= so we don't overcount the spot areas!
      {
        temperature = pow(JJ,(1/8.));
        theta_prime = acos(cos(-zeta)*cos(theta) - cos(phi)*sin(-zeta)*sin(theta));
        if(theta_prime<(M_PI/2)){
          north_total+=temperature*sin(theta);    //figure out which hemisphere to add to
          north_area+=sin(theta);
        } else{
          south_total+=temperature*sin(theta);
          south_area+=sin(theta);
        }
      }

    }
  } //end of theta-phi loops
  double avTempNorth = north_total/north_area;
  double avTempSouth = south_total/south_area;
  // printf("Average temperatures: %f, %f\n", avTempNorth,avTempSouth);

  return (avTempNorth/avTempSouth);
}



// int main(int argc, char **argv)
int main(void)
{
  // double q = strtod(argv[1], NULL);
  // double z = strtod(argv[2], NULL);
  // double JJ = current_density(th*(M_PI/180.0),ph*(M_PI/180.0));

  FILE *output;
  char filename[sizeof "../FINAL_PLOTS/FIG7_Contours/temp_ratios_D.dat"];
  sprintf(filename, "../FINAL_PLOTS/FIG7_Contours/temp_ratios_D.dat");
  output = fopen(filename, "w");

  int N_q=30;
  int N_z=45;  //size of meshes
  double q; double z;
  int qi; int zi;  //dummy variables for iteration
  double stepQ=3.0/N_q;
  double stepZ=90.0/N_z;
  double temperature_ratio;

  double TempRatios[N_q+1][N_z+1];

  for (qi=0;qi<=N_q;qi++){
    q=stepQ*qi;       // theta-position in fine grid
    for (zi=0;zi<=N_z;zi++){
      z=stepZ*zi;       // phi-position in fine grid

      printf("computing temp ratio for q=%1.1f,zeta=%3.0f\r",q,z); fflush(stdout);

      temperature_ratio = temp_contours(q,z);   //enter zeta in degrees
      TempRatios[qi][zi] = temperature_ratio;
      // printf("q=%.2f,z=%.2f, temperature_ratio = %.3f\n", q,z,temperature_ratio);
      fprintf(output,"%e %e %e\n",q,z,temperature_ratio);

    }
  }

  // double temperature_ratio = temp_contours(q,z);   //enter zeta in degrees
  // printf("%f\n", temperature_ratio);

  fclose(output);
  // return 0;
}
