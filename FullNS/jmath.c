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

  - Should the MOI really be 2/5, or do we need corrections from Michi's paper?

  - Double-check the 'northsouth' plus-minus factor in Lambda
*/

{
  double iota = -the_spot;  //inclination of the magnetic dipole axis.
                            //negative because of how the rotation is done

  // double M = 0.25; double r = 1; double Omega = 0.1; //test values [C=0.5,epsilon=0.1]
  // double M = 1; double r = 5.8; double Omega = 0.012; //test

  double M = 1;  //(see note above)
  double r = rns/(1.47662*mns); // divide by (G*Msolar/c^2) to make it dimensionless,
                                // multiply by (1000km/1m), and divide by mns
  double Omega = Omegans;       // defined in io.c

  double mu = 1;
  double II = 2./5;  //dimensionless moment of inertia (MOI)
  double MOI = II*M*r*r;
  double OmegaZ = Omega*(2*MOI)/(r*r*r);  //"frame-dragging frequency"
  double f = 1 - (2*M)/r;

  // optimize code by computing the trig functions beforehand
  double sin_theta = sin(theta);
  double cos_theta = cos(theta);
  double sin_phi = sin(phi);
  double cos_phi = cos(phi);
  double sin_iota = sin(iota);
  double cos_iota = cos(iota);
  double log_f = log(f);

  double theta_prime = acos(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta);
  double phi_prime = atan2((sin_theta*sin_phi), (cos_theta*sin_iota + cos_iota*cos_phi*sin_theta));

  double alpha = (-3*mu*(3 - 4*f + f*f + 2*log_f)*pow(sin(theta_prime),2))/(2.*pow(1 - f,3)*r);
  double beta = phi_prime;
  double alpha0 = sqrt(1.5)*mu*Omega*(1 + sin_iota*sin_iota/5.);

  double Jsquared;

  //outside the spacelike region, no 3-current is allowed
  if (alpha>=alpha0 || theta_prime>(M_PI/2)) //only one spot instead of two
  {
    Jsquared = 0.0;
  }
  else
  {

    double Jt; double Jr; double Jtheta; double Jphi;

    // int northsouth;
    // if (theta_prime<=M_PI) northsouth = -1; else northsouth = 1;
    double Lambda = 2*Omega*(j0(2*asin(sqrt(alpha/alpha0)))*cos_iota -
      j1(2*asin(sqrt(alpha/alpha0)))*cos(beta)*sin_iota);  //j0 and j1 are the bessel functions


// massive formulae for the four-current in a rotated coordinate system! //


    Jt = ((Omega - OmegaZ)*sqrt(f)*(-(((3*mu*r*r*(3 - 4*(f) + f*f +
                     2*log_f)*sin_iota*sin_theta*(-(cos_theta*cos_phi*sin_iota)
                     - cos_iota*sin_theta)*sin_phi)/(8.*M*M*M) +
                (3*mu*r*r*cos_theta*(3 - 4*(f) + f*f + 2*log_f)*
                   sin_iota*(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta)*sin_phi)/
                 (8.*M*M*M))*((cos_theta*(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta)*
                   sin_phi)/(pow(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta,2) +
                   sin_theta*sin_theta*sin_phi*sin_phi) -
                (sin_theta*(cos_iota*cos_theta*cos_phi - sin_iota*sin_theta)*sin_phi)/
                 (pow(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta,2) +
                   sin_theta*sin_theta*sin_phi*sin_phi))) +
           (3*mu*r*r*(3 - 4*(f) + f*f +
                2*log_f)*(-(cos_theta*cos_phi*sin_iota) - cos_iota*sin_theta)*
              (cos_iota*cos_theta - cos_phi*sin_iota*sin_theta)*
              (-((cos_phi*sin_theta*(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta)*
                     (2*(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta)*
                        (cos_iota*cos_theta*cos_phi - sin_iota*sin_theta) +
                       2*cos_theta*sin_theta*sin_phi*sin_phi))/
                   pow(pow(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta,2) +
                     sin_theta*sin_theta*sin_phi*sin_phi,2)) -
                (cos_iota*sin_theta*sin_theta*sin_phi*sin_phi*
                   (2*(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta)*
                      (cos_iota*cos_theta*cos_phi - sin_iota*sin_theta) +
                     2*cos_theta*sin_theta*sin_phi*sin_phi))/
                 pow(pow(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta,2) +
                   sin_theta*sin_theta*sin_phi*sin_phi,2) +
                (cos_theta*cos_phi*(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta))/
                 (pow(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta,2) +
                   sin_theta*sin_theta*sin_phi*sin_phi) +
                (cos_phi*sin_theta*(cos_iota*cos_theta*cos_phi - sin_iota*sin_theta))/
                 (pow(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta,2) +
                   sin_theta*sin_theta*sin_phi*sin_phi) +
                (2*cos_iota*cos_theta*sin_theta*sin_phi*sin_phi)/
                 (pow(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta,2) +
                   sin_theta*sin_theta*sin_phi*sin_phi)))/(8.*M*M*M) +
           ((1/sin_theta)*((3*mu*r*r*(3 - 4*(f) + f*f + 2*log_f)*sin_theta*
                    pow(-(cos_theta*cos_phi*sin_iota) - cos_iota*sin_theta,2))/(8.*M*M*M)\
                  + (3*mu*r*r*cos_theta*(3 - 4*(f) + f*f + 2*log_f)*
                    (-(cos_theta*cos_phi*sin_iota) - cos_iota*sin_theta)*
                    (cos_iota*cos_theta - cos_phi*sin_iota*sin_theta))/(8.*M*M*M) +
                 (3*mu*r*r*(3 - 4*(f) + f*f + 2*log_f)*sin_theta*
                    (cos_iota*cos_theta - cos_phi*sin_iota*sin_theta)*
                    (-(cos_iota*cos_theta) + cos_phi*sin_iota*sin_theta))/(8.*M*M*M)) +
              (f)*(r*r*((-3*mu*((-8*M)/r*r + (4*M)/((f)*r*r) +
               (4*M*(f))/r*r)*r*(1 - pow(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta,2)))/
                     (4.*M*M*M) - (3*mu*((8*M*M)/r*r*r*r -
                         (8*M*M)/(f*f*r*r*r*r) + (16*M)/r*r*r - (8*M)/((f)*r*r*r) -
                         (8*M*(f))/r*r*r)*r*r*
                       (1 - pow(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta,2)))/
                     (16.*M*M*M) - (3*mu*(3 - 4*(f) + f*f + 2*log_f)*
                       (1 - pow(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta,2)))/
                     (8.*M*M*M)) + 2*r*((-3*mu*((-8*M)/r*r + (4*M)/((f)*r*r) +
                         (4*M*(f))/r*r)*r*r*
                       (1 - pow(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta,2)))/
                     (16.*M*M*M) - (3*mu*r*(3 - 4*(f) + f*f +
                         2*log_f)*(1 - pow(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta,2)))/
                     (8.*M*M*M))))*
            ((cos_phi*sin_theta*(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta))/
               (pow(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta,2) +
                 sin_theta*sin_theta*sin_phi*sin_phi) +
              (cos_iota*sin_theta*sin_theta*sin_phi*sin_phi)/
               (pow(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta,2) +
                 sin_theta*sin_theta*sin_phi*sin_phi)) -
           (3*mu*r*r*(3 - 4*(f) + f*f +
                2*log_f)*sin_iota*(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta)*
              sin_phi*(sin_theta*(-((cos_theta*(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta)*sin_phi*
                        (2*(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta)*
                           (cos_iota*cos_theta*cos_phi - sin_iota*sin_theta) +
                          2*cos_theta*sin_theta*sin_phi*sin_phi))/
                      pow(pow(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta,2) +
                        sin_theta*sin_theta*sin_phi*sin_phi,2)) +
                   (sin_theta*(cos_iota*cos_theta*cos_phi - sin_iota*sin_theta)*sin_phi*
                      (2*(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta)*
                         (cos_iota*cos_theta*cos_phi - sin_iota*sin_theta) +
                        2*cos_theta*sin_theta*sin_phi*sin_phi))/
                    pow(pow(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta,2) +
                      sin_theta*sin_theta*sin_phi*sin_phi,2) -
                   (sin_theta*(-(cos_theta*sin_iota) - cos_iota*cos_phi*sin_theta)*sin_phi)/
                    (pow(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta,2) +
                      sin_theta*sin_theta*sin_phi*sin_phi) -
                   (sin_theta*(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta)*sin_phi)/
                    (pow(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta,2) +
                      sin_theta*sin_theta*sin_phi*sin_phi)) +
                cos_theta*((cos_theta*(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta)*sin_phi)/
                    (pow(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta,2) +
                      sin_theta*sin_theta*sin_phi*sin_phi) -
                   (sin_theta*(cos_iota*cos_theta*cos_phi - sin_iota*sin_theta)*sin_phi)/
                    (pow(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta,2) +
                      sin_theta*sin_theta*sin_phi*sin_phi))))/(8.*M*M*M)))/(r*(-2*M + r));

    // //outside the spacelike region, no 3-current is allowed
    // if (alpha>alpha0){
    //   Jr = 0.0; Jtheta = 0.0; Jphi = 0.0;
    // }
    // else{

    Jr = (Lambda*(1/sin_theta)*((-3*mu*r*r*(3 - 4*(f) + f*f +
                2*log_f)*sin_iota*sin_theta*
              (cos_iota*cos_theta - cos_phi*sin_iota*sin_theta)*sin_phi*
              ((cos_theta*(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta)*sin_phi)/
                 (pow(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta,2) +
                   sin_theta*sin_theta*sin_phi*sin_phi) -
                (sin_theta*(cos_iota*cos_theta*cos_phi - sin_iota*sin_theta)*sin_phi)/
                 (pow(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta,2) +
                   sin_theta*sin_theta*sin_phi*sin_phi)))/(8.*M*M*M) +
           (3*mu*r*r*(3 - 4*(f) + f*f +
                2*log_f)*(-(cos_theta*cos_phi*sin_iota) - cos_iota*sin_theta)*
              (cos_iota*cos_theta - cos_phi*sin_iota*sin_theta)*
              ((cos_phi*sin_theta*(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta))/
                 (pow(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta,2) +
                   sin_theta*sin_theta*sin_phi*sin_phi) +
                (cos_iota*sin_theta*sin_theta*sin_phi*sin_phi)/
                 (pow(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta,2) +
                   sin_theta*sin_theta*sin_phi*sin_phi)))/(8.*M*M*M)))/
       (r*sqrt(r*(-2*M + r)));


    Jtheta = -((Lambda*(1/sin_theta)*((-3*mu*((-8*M)/r*r + (4*M)/((f)*r*r) +
                  (4*M*(f))/r*r)*r*r*
                (1 - pow(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta,2)))/(16.*M*M*M)
              - (3*mu*r*(3 - 4*(f) + f*f + 2*log_f)*
                (1 - pow(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta,2)))/(8.*M*M*M))*
           ((cos_phi*sin_theta*(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta))/
              (pow(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta,2) +
                sin_theta*sin_theta*sin_phi*sin_phi) +
             (cos_iota*sin_theta*sin_theta*sin_phi*sin_phi)/
              (pow(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta,2) +
                sin_theta*sin_theta*sin_phi*sin_phi)))/r);


    Jphi = (Lambda*((-3*mu*((-8*M)/r*r + (4*M)/((f)*r*r) + (4*M*(f))/r*r)*r*r*
              (1 - pow(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta,2)))/(16.*M*M*M) -
           (3*mu*r*(3 - 4*(f) + f*f + 2*log_f)*
              (1 - pow(cos_iota*cos_theta - cos_phi*sin_iota*sin_theta,2)))/(8.*M*M*M))*
         ((cos_theta*(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta)*sin_phi)/
            (pow(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta,2) +
              sin_theta*sin_theta*sin_phi*sin_phi) -
           (sin_theta*(cos_iota*cos_theta*cos_phi - sin_iota*sin_theta)*sin_phi)/
            (pow(cos_theta*sin_iota + cos_iota*cos_phi*sin_theta,2) +
              sin_theta*sin_theta*sin_phi*sin_phi)))/r;

    //We're only interested in the 3-current for our model of temperature
    Jsquared = Jr*Jr + Jtheta*Jtheta + Jphi*Jphi;
    // Jsquared = -Jt*Jt + Jr*Jr + Jtheta*Jtheta + Jphi*Jphi;
}

  return Jsquared;
}


// int main(void)
// {
//   double j[4];
//   double mag_j = current_density(1.,0.1,0.1,j);  //call it with (r,theta,phi,j[])
//   printf("current density: \n"); printf("Jt = %e, Jr = %e, Jtheta = %e, Jphi = %e\n",j[0],j[1],j[2],j[3]);
//   // double normj = sqrt(-j[0]*j[0] + j[1]*j[1] + j[2]*j[2] + j[3]*j[3]);
//   printf("J^2 = %e\n",mag_j);
//   return 0;
// }
