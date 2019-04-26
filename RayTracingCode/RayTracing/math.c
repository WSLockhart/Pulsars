/* 
************************************************************
             File: math.c
   Functions for mathematical manipulations
************************************************************
*/

#include "definitions.h"      // File with useful definitions and headers
#include "global.h"           // Global variable headers

double rootbis(double xmin, double xmax, 
	       double (*foo)(double))
/* March 25, 2010 (DP)
   Root finder with bisection. The independent variable
   is in the range [xmin,xmax] in which a solution is known to 
   exist. The function is provided externally as (*foo)(double).
   The iteration stops if either the midpoint of the current 
   interval is the exact solution, or if the fractional width of
   the interval is smaller than the variable EPS.
   The function returns the root. There is a hidden assumption
   in the "while" statement that the root is never at xmid=0.0
*/
{
  double xmid;
  double EPS=1.e-6;               // accuracy

  if ((*foo)(xmin)*(*foo)(xmax)>0)// if the function does not change sign
    {                             // at the initial interval, it's an error
      printf("!No solution or more than one solutions exist");
      printf(" in the bisection interval\n");  
      return 0.0;                 // error code
    }
  do                              // iterate
    {
      xmid=0.5*(xmin+xmax);       // calculate mid point
                                  // if the solution is in the first half  
      if ((*foo)(xmin)*(*foo)(xmid)<0.0)
	{
	  xmax=xmid;              // bring down the upper end
	}
      else
	{
	  xmin=xmid;              // otherwise bring up the lower end
	}
    }                             // as long as the interval is wide and
                                  // haven't hit the root
  while (fabs(xmax-xmin)>EPS*xmid && (*foo)(xmid)!=0.0);

  return xmid;
  
}

