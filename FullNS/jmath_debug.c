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
  // double iota = -the_spot;  //inclination of the magnetic dipole axis.
                            //negative because of how the rotation is done
  double iota = -M_PI/6.;  //inclination of the magnetic dipole axis.

  // double M = 0.25; double r = 1; double Omega = 0.1; //test values [C=0.5,epsilon=0.1]
  double M = 1; double r = 5.8; double Omega = 0.012; //test
  // M *= (0.5); r *= (0.5); Omega *= 2;
  // double M = 1;  //(see note above)
  // double r = rns/(1.47662*mns); // divide by (G*Msolar/c^2) to make it dimensionless,
  //                               // multiply by (1000km/1m), and divide by mns
  // double Omega = Omegans;       // defined in io.c
  double II = 2./5;  //dimensionless moment of inertia (MOI)

  //debug check for dimensionless consistency:
  // r *= (1./2); M *= (1./2); Omega *= 2;

  double mu = 1;
  double MOI = II*M*r*r;
  double OmegaZ = Omega*(2*MOI)/(r*r*r);  //"frame-dragging frequency"
  double C = (2*M)/r;
  double f = 1-C;

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

  double alpha = -(1.5)*((3 + 2*log_f - 4*f + f*f)/(pow(1 - f,3)))*(mu/r)*sin(theta_prime)*sin(theta_prime);
  double beta = phi_prime;
  double alpha0 = sqrt(1.5)*mu*Omega*(1 + sin_iota*sin_iota/5.);

  double Jsquared;

  //outside the spacelike region, no 3-current is allowed
  // if (alpha>=alpha0 || theta_prime>(M_PI/2)){ //only one spot instead of two
  if (alpha>=alpha0){
    Jsquared = 0.0;
    printf("no current detected at (%.1f,%.1f).\n",theta,phi);
  }
  else
  {

    double Jt; double Jr; double Jtheta; double Jphi;

    int plusminus;  // correction to the paper for southern cap
    if (theta<(M_PI/2)) plusminus = 1; else plusminus = -1;
    double Lambda = 2*Omega*(j0(2*asin(sqrt(alpha/alpha0)))*cos_iota -
      plusminus*j1(2*asin(sqrt(alpha/alpha0)))*cos(beta)*sin_iota);  //j0 and j1 are the bessel functions


// massive formulae for the four-current in a rotated coordinate system! //


    Jt = (Omega - OmegaZ)*pow(r,-2)*pow(1 - 2*M*pow(r,-1),-0.5)*(-((cos(theta)*
           pow(pow(sin(phi),2)*pow(sin(theta),2) + pow(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta),2),-1)*sin(phi)*
           (cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta)) -
          pow(pow(sin(phi),2)*pow(sin(theta),2) + pow(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta),2),-1)*sin(phi)*sin(theta)*
           (cos(iota)*cos(phi)*cos(theta) - sin(iota)*sin(theta)))*
        ((3*mu*pow(M,-3)*pow(r,2)*(3 + 2*log(1 - 2*M*pow(r,-1)) - 4*(1 - 2*M*pow(r,-1)) + pow(1 - 2*M*pow(r,-1),2))*sin(iota)*sin(phi)*
             sin(theta)*(-(cos(phi)*cos(theta)*sin(iota)) - cos(iota)*sin(theta)))/8. +
          (3*mu*cos(theta)*pow(M,-3)*pow(r,2)*(3 + 2*log(1 - 2*M*pow(r,-1)) - 4*(1 - 2*M*pow(r,-1)) + pow(1 - 2*M*pow(r,-1),2))*sin(iota)*
             sin(phi)*(cos(iota)*cos(theta) - cos(phi)*sin(iota)*sin(theta)))/8.)) +
     (3*mu*pow(M,-3)*pow(r,2)*(3 + 2*log(1 - 2*M*pow(r,-1)) - 4*(1 - 2*M*pow(r,-1)) + pow(1 - 2*M*pow(r,-1),2))*
        (-(cos(phi)*cos(theta)*sin(iota)) - cos(iota)*sin(theta))*(cos(iota)*cos(theta) - cos(phi)*sin(iota)*sin(theta))*
        (2*cos(iota)*cos(theta)*pow(pow(sin(phi),2)*pow(sin(theta),2) + pow(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta),2),-1)*
           pow(sin(phi),2)*sin(theta) + cos(phi)*cos(theta)*
           pow(pow(sin(phi),2)*pow(sin(theta),2) + pow(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta),2),-1)*
           (cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta)) +
          cos(phi)*pow(pow(sin(phi),2)*pow(sin(theta),2) + pow(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta),2),-1)*sin(theta)*
           (cos(iota)*cos(phi)*cos(theta) - sin(iota)*sin(theta)) -
          cos(iota)*pow(pow(sin(phi),2)*pow(sin(theta),2) + pow(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta),2),-2)*pow(sin(phi),2)*
           pow(sin(theta),2)*(2*cos(theta)*pow(sin(phi),2)*sin(theta) +
             2*(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta))*(cos(iota)*cos(phi)*cos(theta) - sin(iota)*sin(theta))) -
          cos(phi)*pow(pow(sin(phi),2)*pow(sin(theta),2) + pow(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta),2),-2)*sin(theta)*
           (cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta))*
           (2*cos(theta)*pow(sin(phi),2)*sin(theta) + 2*(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta))*
              (cos(iota)*cos(phi)*cos(theta) - sin(iota)*sin(theta)))))/8. +
     (cos(iota)*pow(pow(sin(phi),2)*pow(sin(theta),2) + pow(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta),2),-1)*pow(sin(phi),2)*
         pow(sin(theta),2) + cos(phi)*pow(pow(sin(phi),2)*pow(sin(theta),2) + pow(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta),2),-1)*
         sin(theta)*(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta)))*
      ((1 - 2*M*pow(r,-1))*(pow(r,2)*((-3*mu*pow(M,-3)*pow(r,2)*
                 (8*pow(M,2)*pow(r,-4) + 16*M*pow(r,-3) - 8*M*pow(r,-3)*(1 - 2*M*pow(r,-1)) -
                   8*pow(M,2)*pow(r,-4)*pow(1 - 2*M*pow(r,-1),-2) - 8*M*pow(r,-3)*pow(1 - 2*M*pow(r,-1),-1))*
                 (1 - pow(cos(iota)*cos(theta) - cos(phi)*sin(iota)*sin(theta),2)))/16. -
              (3*mu*r*pow(M,-3)*(-8*M*pow(r,-2) + 4*M*pow(r,-2)*(1 - 2*M*pow(r,-1)) + 4*M*pow(r,-2)*pow(1 - 2*M*pow(r,-1),-1))*
                 (1 - pow(cos(iota)*cos(theta) - cos(phi)*sin(iota)*sin(theta),2)))/4. -
              (3*mu*pow(M,-3)*(3 + 2*log(1 - 2*M*pow(r,-1)) - 4*(1 - 2*M*pow(r,-1)) + pow(1 - 2*M*pow(r,-1),2))*
                 (1 - pow(cos(iota)*cos(theta) - cos(phi)*sin(iota)*sin(theta),2)))/8.) +
           2*r*((-3*mu*pow(M,-3)*pow(r,2)*(-8*M*pow(r,-2) + 4*M*pow(r,-2)*(1 - 2*M*pow(r,-1)) + 4*M*pow(r,-2)*pow(1 - 2*M*pow(r,-1),-1))*
                 (1 - pow(cos(iota)*cos(theta) - cos(phi)*sin(iota)*sin(theta),2)))/16. -
              (3*mu*r*pow(M,-3)*(3 + 2*log(1 - 2*M*pow(r,-1)) - 4*(1 - 2*M*pow(r,-1)) + pow(1 - 2*M*pow(r,-1),2))*
                 (1 - pow(cos(iota)*cos(theta) - cos(phi)*sin(iota)*sin(theta),2)))/8.)) +
        (1./sin(theta))*((3*mu*pow(M,-3)*pow(r,2)*(3 + 2*log(1 - 2*M*pow(r,-1)) - 4*(1 - 2*M*pow(r,-1)) + pow(1 - 2*M*pow(r,-1),2))*
              pow(-(cos(phi)*cos(theta)*sin(iota)) - cos(iota)*sin(theta),2)*sin(theta))/8. +
           (3*mu*cos(theta)*pow(M,-3)*pow(r,2)*(3 + 2*log(1 - 2*M*pow(r,-1)) - 4*(1 - 2*M*pow(r,-1)) + pow(1 - 2*M*pow(r,-1),2))*
              (-(cos(phi)*cos(theta)*sin(iota)) - cos(iota)*sin(theta))*(cos(iota)*cos(theta) - cos(phi)*sin(iota)*sin(theta)))/8. +
           (3*mu*pow(M,-3)*pow(r,2)*(3 + 2*log(1 - 2*M*pow(r,-1)) - 4*(1 - 2*M*pow(r,-1)) + pow(1 - 2*M*pow(r,-1),2))*sin(theta)*
              (cos(iota)*cos(theta) - cos(phi)*sin(iota)*sin(theta))*(-(cos(iota)*cos(theta)) + cos(phi)*sin(iota)*sin(theta)))/8.)) -
     (3*mu*pow(M,-3)*pow(r,2)*(3 + 2*log(1 - 2*M*pow(r,-1)) - 4*(1 - 2*M*pow(r,-1)) + pow(1 - 2*M*pow(r,-1),2))*sin(iota)*sin(phi)*
        (cos(iota)*cos(theta) - cos(phi)*sin(iota)*sin(theta))*
        (cos(theta)*(cos(theta)*pow(pow(sin(phi),2)*pow(sin(theta),2) + pow(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta),2),-1)*
              sin(phi)*(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta)) -
             pow(pow(sin(phi),2)*pow(sin(theta),2) + pow(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta),2),-1)*sin(phi)*sin(theta)*
              (cos(iota)*cos(phi)*cos(theta) - sin(iota)*sin(theta))) +
          sin(theta)*(-(pow(pow(sin(phi),2)*pow(sin(theta),2) + pow(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta),2),-1)*sin(phi)*
                sin(theta)*(-(cos(theta)*sin(iota)) - cos(iota)*cos(phi)*sin(theta))) -
             pow(pow(sin(phi),2)*pow(sin(theta),2) + pow(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta),2),-1)*sin(phi)*sin(theta)*
              (cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta)) -
             cos(theta)*pow(pow(sin(phi),2)*pow(sin(theta),2) + pow(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta),2),-2)*sin(phi)*
              (cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta))*
              (2*cos(theta)*pow(sin(phi),2)*sin(theta) +
                2*(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta))*(cos(iota)*cos(phi)*cos(theta) - sin(iota)*sin(theta))) +
             pow(pow(sin(phi),2)*pow(sin(theta),2) + pow(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta),2),-2)*sin(phi)*sin(theta)*
              (cos(iota)*cos(phi)*cos(theta) - sin(iota)*sin(theta))*
              (2*cos(theta)*pow(sin(phi),2)*sin(theta) +
                2*(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta))*(cos(iota)*cos(phi)*cos(theta) - sin(iota)*sin(theta))))))/8.);


    Jr = Lambda*(1./sin(theta))*pow(r,-1)*pow((1 - 2*M*pow(r,-1))*pow(r,2),-0.5)*
   ((3*mu*pow(M,-3)*pow(r,2)*(3 + 2*log(1 - 2*M*pow(r,-1)) - 4*(1 - 2*M*pow(r,-1)) + pow(1 - 2*M*pow(r,-1),2))*
        (-(cos(phi)*cos(theta)*sin(iota)) - cos(iota)*sin(theta))*(cos(iota)*cos(theta) - cos(phi)*sin(iota)*sin(theta))*
        (cos(iota)*pow(pow(sin(phi),2)*pow(sin(theta),2) + pow(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta),2),-1)*pow(sin(phi),2)*
           pow(sin(theta),2) + cos(phi)*pow(pow(sin(phi),2)*pow(sin(theta),2) + pow(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta),2),
            -1)*sin(theta)*(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta))))/8. -
     (3*mu*pow(M,-3)*pow(r,2)*(3 + 2*log(1 - 2*M*pow(r,-1)) - 4*(1 - 2*M*pow(r,-1)) + pow(1 - 2*M*pow(r,-1),2))*sin(iota)*sin(phi)*sin(theta)*
        (cos(iota)*cos(theta) - cos(phi)*sin(iota)*sin(theta))*
        (cos(theta)*pow(pow(sin(phi),2)*pow(sin(theta),2) + pow(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta),2),-1)*sin(phi)*
           (cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta)) -
          pow(pow(sin(phi),2)*pow(sin(theta),2) + pow(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta),2),-1)*sin(phi)*sin(theta)*
           (cos(iota)*cos(phi)*cos(theta) - sin(iota)*sin(theta))))/8.);


    Jtheta = Lambda*(1./sin(theta))*pow(r,-1)*((3*mu*pow(M,-3)*pow(r,2)*(-8*M*pow(r,-2) + 4*M*pow(r,-2)*(1 - 2*M*pow(r,-1)) +
          4*M*pow(r,-2)*pow(1 - 2*M*pow(r,-1),-1))*(1 - pow(cos(iota)*cos(theta) - cos(phi)*sin(iota)*sin(theta),2)))/16. +
     (3*mu*r*pow(M,-3)*(3 + 2*log(1 - 2*M*pow(r,-1)) - 4*(1 - 2*M*pow(r,-1)) + pow(1 - 2*M*pow(r,-1),2))*
        (1 - pow(cos(iota)*cos(theta) - cos(phi)*sin(iota)*sin(theta),2)))/8.)*
   (cos(iota)*pow(pow(sin(phi),2)*pow(sin(theta),2) + pow(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta),2),-1)*pow(sin(phi),2)*
      pow(sin(theta),2) + cos(phi)*pow(pow(sin(phi),2)*pow(sin(theta),2) + pow(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta),2),-1)*
      sin(theta)*(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta)));


    Jphi = Lambda*pow(r,-1)*((-3*mu*pow(M,-3)*pow(r,2)*(-8*M*pow(r,-2) + 4*M*pow(r,-2)*(1 - 2*M*pow(r,-1)) + 4*M*pow(r,-2)*pow(1 - 2*M*pow(r,-1),-1))*
        (1 - pow(cos(iota)*cos(theta) - cos(phi)*sin(iota)*sin(theta),2)))/16. -
     (3*mu*r*pow(M,-3)*(3 + 2*log(1 - 2*M*pow(r,-1)) - 4*(1 - 2*M*pow(r,-1)) + pow(1 - 2*M*pow(r,-1),2))*
        (1 - pow(cos(iota)*cos(theta) - cos(phi)*sin(iota)*sin(theta),2)))/8.)*
   (cos(theta)*pow(pow(sin(phi),2)*pow(sin(theta),2) + pow(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta),2),-1)*sin(phi)*
      (cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta)) -
     pow(pow(sin(phi),2)*pow(sin(theta),2) + pow(cos(theta)*sin(iota) + cos(iota)*cos(phi)*sin(theta),2),-1)*sin(phi)*sin(theta)*
      (cos(iota)*cos(phi)*cos(theta) - sin(iota)*sin(theta)));

    //Jsquared = Jr*Jr + Jtheta*Jtheta + Jphi*Jphi;
    Jsquared = -Jt*Jt + Jr*Jr + Jtheta*Jtheta + Jphi*Jphi;

    printf("iota = %.2f. At (%.2f,%.2f) Current:\n Jt = %e\n Jr = %e\n Jtheta = %e\n Jphi = %e\n" \
      ,iota,theta,phi, Jt,Jr,Jtheta,Jphi);

}

  return Jsquared;
}


// int main(int argc, char **argv)
int main(void)
{
  // th = strtod(argv[1], NULL);
  // ph = strtod(argv[2], NULL);
  double mag_j = current_density(M_PI-((M_PI/6.)+0.2),M_PI+0.2);

  return 0;
}
