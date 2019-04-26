/*
************************************************************
             File: integrals
   Functions that set the grid of energies and perform the
   various integrals to convert the result of ray-tracing
   into a transfer function.
************************************************************
*/

#include "definitions.h"      // File with useful definitions and headers
#include "global.h"           // Global variable headers

/*
****************************************************************
   Protypes of external functions used in this file
**************************************************************** */

void set_energy_grid(double Energy[])
/* April 21, 2010 (DP)
   Sets a linear grid over photon energies and stores it in the array
   Energy[Nen]. The minimum and maximum values of energy, Emin
   and Emax, as well as the number of grid points Nen can be found
   in the global variables.
*/
{
  int index;
  for (index=1;index<=Nen;index++)
    {
      Energy[index-1]=Emin+(Emax-Emin)/(Nen-1.0)*(index-1.0);
      // Energy[index-1]=pow(10.,log10(Emin)+(log10(Emax)-log10(Emin))/(Nen-1.0)*(index-1.0));
    }

  return;
}
