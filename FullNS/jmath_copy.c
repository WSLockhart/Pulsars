#include <stdio.h>
#include <math.h>
#include "definitions.h"      // File with useful definitions and headers
#include "global.h"           // Global variable headers


double current_density(double r, double theta, double phi)
/* July 11, 2018 (WL)
   Calculate the current density J(theta,phi) on the surface of the star.
   This is based on Sam Gralla's analytic results, coded originally in
   Mathematica and transferred here.
*/
{
/* Aug 15, 2018 (WL)
Still need to sync these parameters with the global inputs, and
also find the way to code the moment of inertia consistently (see Michi's paper)
*/
  double iota = -the_spot;  //inclination of the magnetic dipole axis.
                            //negative because of how the rotation is done
  // double M = mns; // solar mass is the fundamental unit
  // double Rstar = rns / 1.47662; // divide by G*Msolar/c^2 to make it dimensionless
  // double Omega = Omegans;  //global variable with units of 1/M
  double M = 0.25; // solar mass is the fundamental unit
  double Rstar = 1; // divide by G*Msolar/c^2 to make it dimensionless
  double Omega = 0.1;  //global variable with units of 1/M

  double mu = 1;
  double II = 2./5;
  double MOI = II*M*Rstar*Rstar;
  double OmegaZ = Omega*(2*MOI)/(r*r*r);  //"frame-dragging frequency"
  double f = 1 - (2*M)/r;

  double theta_prime = acos(cos(iota)*cos(theta) - cos(phi)*sin(iota)*sin(theta));
  double phi_prime = atan2((sin(theta)*sin(phi)), (cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta)));
  // if (phi > M_PI) phi_prime+=(2*M_PI);

  double alpha = (-3*mu*(3 - 4*f + pow(f,2) + 2*log(f))*pow(sin(theta_prime),2))/(2.*pow(1 - f,3)*r);
  double beta = phi_prime;
  double alpha0 = sqrt(1.5)*mu*Omega*(1 + pow(sin(iota),2)/5.);

  double Jsquared; double Jt; double Jr; double Jtheta; double Jphi;

  //outside the spacelike region, no 3-current is allowed
  if (alpha>=alpha0){
    Jsquared = 0.0;
  }
  else
  {

    int northsouth;
    if (theta_prime<=M_PI) northsouth = -1; else northsouth = 1;
    double Lambda = northsouth*2*Omega*(j0(2*asin(sqrt(alpha/alpha0)))*cos(iota) -
      j1(2*asin(sqrt(alpha/alpha0)))*cos(beta)*sin(iota));  //j0 and j1 are the bessel functions


// massive formulae for the four-current in a rotated coordinate system! //

    Jt = ((Omega - OmegaZ)*sqrt(1 - (2*M)/r)*
         (-(((3*mu*pow(r,2)*(3 - 4*(1 - (2*M)/r) + pow(1 - (2*M)/r,2) +
                     2*log(1 - (2*M)/r))*sin(iota)*sin(theta)*
                   (-(cos(theta)*cos(phi)*sin(iota)) - cos(iota)*sin(theta))*sin(phi))/(8.*pow(M,3)) +
                (3*mu*pow(r,2)*cos(theta)*
                   (3 - 4*(1 - (2*M)/r) + pow(1 - (2*M)/r,2) + 2*log(1 - (2*M)/r))*
                   sin(iota)*(cos(iota)*cos(theta) - cos(phi)*sin(iota)*sin(theta))*sin(phi))/
                 (8.*pow(M,3)))*((cos(theta)*(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta))*
                   sin(phi))/
                 (pow(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta),2) +
                   pow(sin(theta),2)*pow(sin(phi),2)) -
                (sin(theta)*(cos(iota)*cos(theta)*cos(phi) - sin(iota)*sin(theta))*sin(phi))/
                 (pow(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta),2) +
                   pow(sin(theta),2)*pow(sin(phi),2)))) +
           (3*mu*pow(r,2)*(3 - 4*(1 - (2*M)/r) + pow(1 - (2*M)/r,2) +
                2*log(1 - (2*M)/r))*(-(cos(theta)*cos(phi)*sin(iota)) - cos(iota)*sin(theta))*
              (cos(iota)*cos(theta) - cos(phi)*sin(iota)*sin(theta))*
              (-((cos(phi)*sin(theta)*(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta))*
                     (2*(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta))*
                        (cos(iota)*cos(theta)*cos(phi) - sin(iota)*sin(theta)) +
                       2*cos(theta)*sin(theta)*pow(sin(phi),2)))/
                   pow(pow(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta),2) +
                     pow(sin(theta),2)*pow(sin(phi),2),2)) -
                (cos(iota)*pow(sin(theta),2)*pow(sin(phi),2)*
                   (2*(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta))*
                      (cos(iota)*cos(theta)*cos(phi) - sin(iota)*sin(theta)) +
                     2*cos(theta)*sin(theta)*pow(sin(phi),2)))/
                 pow(pow(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta),2) +
                   pow(sin(theta),2)*pow(sin(phi),2),2) +
                (cos(theta)*cos(phi)*(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta)))/
                 (pow(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta),2) +
                   pow(sin(theta),2)*pow(sin(phi),2)) +
                (cos(phi)*sin(theta)*(cos(iota)*cos(theta)*cos(phi) - sin(iota)*sin(theta)))/
                 (pow(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta),2) +
                   pow(sin(theta),2)*pow(sin(phi),2)) +
                (2*cos(iota)*cos(theta)*sin(theta)*pow(sin(phi),2))/
                 (pow(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta),2) +
                   pow(sin(theta),2)*pow(sin(phi),2))))/(8.*pow(M,3)) +
           ((1/sin(theta))*((3*mu*pow(r,2)*(3 - 4*(1 - (2*M)/r) + pow(1 - (2*M)/r,2) +
                      2*log(1 - (2*M)/r))*sin(theta)*
                    pow(-(cos(theta)*cos(phi)*sin(iota)) - cos(iota)*sin(theta),2))/(8.*pow(M,3))\
                  + (3*mu*pow(r,2)*cos(theta)*
                    (3 - 4*(1 - (2*M)/r) + pow(1 - (2*M)/r,2) + 2*log(1 - (2*M)/r))*
                    (-(cos(theta)*cos(phi)*sin(iota)) - cos(iota)*sin(theta))*
                    (cos(iota)*cos(theta) - cos(phi)*sin(iota)*sin(theta)))/(8.*pow(M,3)) +
                 (3*mu*pow(r,2)*(3 - 4*(1 - (2*M)/r) + pow(1 - (2*M)/r,2) +
                      2*log(1 - (2*M)/r))*sin(theta)*
                    (cos(iota)*cos(theta) - cos(phi)*sin(iota)*sin(theta))*
                    (-(cos(iota)*cos(theta)) + cos(phi)*sin(iota)*sin(theta)))/(8.*pow(M,3))) +
              (1 - (2*M)/r)*(pow(r,2)*
                  ((-3*mu*((-8*M)/pow(r,2) + (4*M)/((1 - (2*M)/r)*pow(r,2)) +
                         (4*M*(1 - (2*M)/r))/pow(r,2))*r*
                       (1 - pow(cos(iota)*cos(theta) - cos(phi)*sin(iota)*sin(theta),2)))/
                     (4.*pow(M,3)) - (3*mu*
                       ((8*pow(M,2))/pow(r,4) -
                         (8*pow(M,2))/(pow(1 - (2*M)/r,2)*pow(r,4)) +
                         (16*M)/pow(r,3) - (8*M)/((1 - (2*M)/r)*pow(r,3)) -
                         (8*M*(1 - (2*M)/r))/pow(r,3))*pow(r,2)*
                       (1 - pow(cos(iota)*cos(theta) - cos(phi)*sin(iota)*sin(theta),2)))/
                     (16.*pow(M,3)) -
                    (3*mu*(3 - 4*(1 - (2*M)/r) + pow(1 - (2*M)/r,2) + 2*log(1 - (2*M)/r))*
                       (1 - pow(cos(iota)*cos(theta) - cos(phi)*sin(iota)*sin(theta),2)))/
                     (8.*pow(M,3))) +
                 2*r*((-3*mu*((-8*M)/pow(r,2) + (4*M)/((1 - (2*M)/r)*pow(r,2)) +
                         (4*M*(1 - (2*M)/r))/pow(r,2))*pow(r,2)*
                       (1 - pow(cos(iota)*cos(theta) - cos(phi)*sin(iota)*sin(theta),2)))/
                     (16.*pow(M,3)) -
                    (3*mu*r*(3 - 4*(1 - (2*M)/r) + pow(1 - (2*M)/r,2) +
                         2*log(1 - (2*M)/r))*
                       (1 - pow(cos(iota)*cos(theta) - cos(phi)*sin(iota)*sin(theta),2)))/
                     (8.*pow(M,3)))))*
            ((cos(phi)*sin(theta)*(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta)))/
               (pow(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta),2) +
                 pow(sin(theta),2)*pow(sin(phi),2)) +
              (cos(iota)*pow(sin(theta),2)*pow(sin(phi),2))/
               (pow(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta),2) +
                 pow(sin(theta),2)*pow(sin(phi),2))) -
           (3*mu*pow(r,2)*(3 - 4*(1 - (2*M)/r) + pow(1 - (2*M)/r,2) +
                2*log(1 - (2*M)/r))*sin(iota)*(cos(iota)*cos(theta) - cos(phi)*sin(iota)*sin(theta))*
              sin(phi)*(sin(theta)*(-((cos(theta)*(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta))*sin(phi)*
                        (2*(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta))*
                           (cos(iota)*cos(theta)*cos(phi) - sin(iota)*sin(theta)) +
                          2*cos(theta)*sin(theta)*pow(sin(phi),2)))/
                      pow(pow(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta),2) +
                        pow(sin(theta),2)*pow(sin(phi),2),2)) +
                   (sin(theta)*(cos(iota)*cos(theta)*cos(phi) - sin(iota)*sin(theta))*sin(phi)*
                      (2*(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta))*
                         (cos(iota)*cos(theta)*cos(phi) - sin(iota)*sin(theta)) +
                        2*cos(theta)*sin(theta)*pow(sin(phi),2)))/
                    pow(pow(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta),2) +
                      pow(sin(theta),2)*pow(sin(phi),2),2) -
                   (sin(theta)*(-(cos(theta)*sin(iota)) - cos(iota)*cos(phi)*sin(theta))*sin(phi))/
                    (pow(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta),2) +
                      pow(sin(theta),2)*pow(sin(phi),2)) -
                   (sin(theta)*(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta))*sin(phi))/
                    (pow(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta),2) +
                      pow(sin(theta),2)*pow(sin(phi),2))) +
                cos(theta)*((cos(theta)*(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta))*sin(phi))/
                    (pow(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta),2) +
                      pow(sin(theta),2)*pow(sin(phi),2)) -
                   (sin(theta)*(cos(iota)*cos(theta)*cos(phi) - sin(iota)*sin(theta))*sin(phi))/
                    (pow(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta),2) +
                      pow(sin(theta),2)*pow(sin(phi),2)))))/(8.*pow(M,3))))/(r*(-2*M + r));


    // //outside the spacelike region, no 3-current is allowed
    // if (alpha>alpha0){
    //   Jr = 0.0; Jtheta = 0.0; Jphi = 0.0;
    // }
    // else{

    Jr = (Lambda*(1/sin(theta))*((-3*mu*pow(r,2)*(3 - 4*(1 - (2*M)/r) + pow(1 - (2*M)/r,2) +
                2*log(1 - (2*M)/r))*sin(iota)*sin(theta)*
              (cos(iota)*cos(theta) - cos(phi)*sin(iota)*sin(theta))*sin(phi)*
              ((cos(theta)*(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta))*sin(phi))/
                 (pow(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta),2) +
                   pow(sin(theta),2)*pow(sin(phi),2)) -
                (sin(theta)*(cos(iota)*cos(theta)*cos(phi) - sin(iota)*sin(theta))*sin(phi))/
                 (pow(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta),2) +
                   pow(sin(theta),2)*pow(sin(phi),2))))/(8.*pow(M,3)) +
           (3*mu*pow(r,2)*(3 - 4*(1 - (2*M)/r) + pow(1 - (2*M)/r,2) +
                2*log(1 - (2*M)/r))*(-(cos(theta)*cos(phi)*sin(iota)) - cos(iota)*sin(theta))*
              (cos(iota)*cos(theta) - cos(phi)*sin(iota)*sin(theta))*
              ((cos(phi)*sin(theta)*(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta)))/
                 (pow(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta),2) +
                   pow(sin(theta),2)*pow(sin(phi),2)) +
                (cos(iota)*pow(sin(theta),2)*pow(sin(phi),2))/
                 (pow(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta),2) +
                   pow(sin(theta),2)*pow(sin(phi),2))))/(8.*pow(M,3))))/
       (r*sqrt(r*(-2*M + r)));


    Jtheta = -((Lambda*(1/sin(theta))*((-3*mu*((-8*M)/pow(r,2) + (4*M)/((1 - (2*M)/r)*pow(r,2)) +
                  (4*M*(1 - (2*M)/r))/pow(r,2))*pow(r,2)*
                (1 - pow(cos(iota)*cos(theta) - cos(phi)*sin(iota)*sin(theta),2)))/(16.*pow(M,3))\
              - (3*mu*r*(3 - 4*(1 - (2*M)/r) + pow(1 - (2*M)/r,2) + 2*log(1 - (2*M)/r))*
                (1 - pow(cos(iota)*cos(theta) - cos(phi)*sin(iota)*sin(theta),2)))/(8.*pow(M,3)))*
           ((cos(phi)*sin(theta)*(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta)))/
              (pow(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta),2) +
                pow(sin(theta),2)*pow(sin(phi),2)) +
             (cos(iota)*pow(sin(theta),2)*pow(sin(phi),2))/
              (pow(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta),2) +
                pow(sin(theta),2)*pow(sin(phi),2))))/r);


    Jphi = (Lambda*((-3*mu*((-8*M)/pow(r,2) + (4*M)/((1 - (2*M)/r)*pow(r,2)) +
                (4*M*(1 - (2*M)/r))/pow(r,2))*pow(r,2)*
              (1 - pow(cos(iota)*cos(theta) - cos(phi)*sin(iota)*sin(theta),2)))/(16.*pow(M,3)) -
           (3*mu*r*(3 - 4*(1 - (2*M)/r) + pow(1 - (2*M)/r,2) + 2*log(1 - (2*M)/r))*
              (1 - pow(cos(iota)*cos(theta) - cos(phi)*sin(iota)*sin(theta),2)))/(8.*pow(M,3)))*
         ((cos(theta)*(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta))*sin(phi))/
            (pow(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta),2) +
              pow(sin(theta),2)*pow(sin(phi),2)) -
           (sin(theta)*(cos(iota)*cos(theta)*cos(phi) - sin(iota)*sin(theta))*sin(phi))/
            (pow(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta),2) +
              pow(sin(theta),2)*pow(sin(phi),2))))/r;


    // current[0]=Jt; current[1]=Jr; current[2]=Jtheta; current[3]=Jphi;
    Jsquared = -Jt*Jt + Jr*Jr + Jtheta*Jtheta + Jphi*Jphi;
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
