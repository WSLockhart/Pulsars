#include <stdio.h>
#include <math.h>
#include "definitions.h"      // File with useful definitions and headers
#include "global.h"           // Global variable headers



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

  Should the MOI really be 2/5, or do we need corrections from Michi's paper?
*/

{
  double iota = -the_spot;  //inclination of the magnetic dipole axis.
                            //negative because of how the rotation is done
  double q = qns;  //ratio of quadrupole moment to dipole moment

  // double M = 1; double R = 5.8; double Omega = 0.012; //test
  //debug check for dimensionless consistency:
// M *= (0.5); R *= (0.5); Omega *= 2;
  double M = 1;  //(see note above)
  double R = rns/(1.47662*mns); // divide by (G*Msolar/c^2) to make it dimensionless,
                                // multiply by (1000km/1m), and divide by mns
  double Omega = Omegans;       // defined in io.c
  double II = 2./5;  //dimensionless moment of inertia (MOI)

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
  // double mu2 = (9*mu*q*M*M*pow(R,-1)*(3 + 2*log_f - 4*f + f*f)*
  //    pow(17 - 9*f + 6*log_f*(1 + 3*f) - 9*f*f + pow(f,3),-1))/5.;

  double theta_prime = acos(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta);
  double phi_prime = atan2((sin_theta*sin_phi), (cos_theta*sin_iota + cos_iota*cos_phi*sin_theta));

  double R1 = -(1.5)*((3 + 2*log_f - 4*f + f*f)/(pow(1 - f,3)))/R;
  double alpha = mu*R1*(1+q*cos(theta_prime))*sin(theta_prime)*sin(theta_prime);  //this is alpha ON THE SURFACE (r=R)
  double beta = phi_prime;
  double alpha0 = sqrt(1.5)*mu*Omega*(1 + sin_iota*sin_iota/5.);

  double Jsquared;

  int northcap;  // correction to the paper for southern cap
  if (theta_prime<(M_PI/2)) northcap = 1; else northcap = -1;

  //after adding the quadrupole term, we must also restrict alpha to positive values
  if (alpha>=alpha0 || alpha<=0){
    Jsquared = 0.0;
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

    // spot type: full hotspot
    // Jsquared=JmuJmu;

    //spot type (C): only 3-current hotspot
    if(JmuJmu>0){
      Jsquared = Jr*Jr + Jtheta*Jtheta + Jphi*Jphi;
    } else{
      Jsquared=0;
    }

    //spot type (C) for hotspot plots:
    // if(JmuJmu>0){
    //   Jsquared = Jr*Jr + Jtheta*Jtheta + Jphi*Jphi;
    // } else{
    //   Jsquared=-1;
    // }

    // spot type (B): uniform spacelike region
    // if (JmuJmu > 0){
    //   Jsquared = 0.001;
    // } else{
    //   Jsquared = 0;
    // }

    // spot type (A): uniform polar cap
    // Jsquared = 0.001;
}

  return Jsquared;

}


// // int main(int argc, char **argv)
// int main(void)
// {
//   // th = strtod(argv[1], NULL);
//   // ph = strtod(argv[2], NULL);
//   double mag_j = current_density(M_PI/3+0.1,0.1);
//
//   return 0;
// }
