
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**
**  FILE:      c_shwave.c
**
**  FUNCTION:  double c_shwave()
**
**  PURPOSE:   This code was extracted from the petfunc function from the
**             subwatr.f Soilwater Model source code file.
**             Calculate the short wave radiation outside the atmosphere using
**             Pennmans equation (1948)
**
**  INPUTS:
**    jday      - current Julian day (1-366)
**    month     - current month (1-12)
**    rlatitude - latitude of the site (in radians)
**
**  GLOBAL VARIABLES:
**    None
**
**  LOCAL VARIABLES:
**    ahou            - ?
**    declin          - declination (radians)
**    par1, par2      - parameters in computation of ahou
**    solrad          - solar radiation (ly/day)
**    transcof(month) - transmission coefficient for the month
**
**  OUTPUTS:
**    c_shwave - short wave solar radiation outside the atmosphere
**
**  CALLED BY:
**    snowCent()
**    snowmodel()
**
**  CALLS:
**    None
**
*****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "swconst.h"

    double c_shwave(int month, double rlatitude, int jday)
    {

      double  ahou;
      double  declin;
      double  par1;
      double  par2;
      double  solrad;
      double temp;
      static double transcof[] = {0.0,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8};

      /* Calculate the short wave solar radiation on a clear day using a */
      /* equation presented by Sellers(1965) */

      declin=0.401426 *sin(6.283185*(jday-77.0)/365.0);
      temp = 1.0-pow(-tan(rlatitude)*tan(declin),2);
      if (temp < 0.0) {
        temp = 0.0;
      }
      par1=sqrt(temp);
      par2=((-tan(rlatitude)*tan(declin)));

      ahou=atan2(par1,par2);
      if(ahou < 0.0) {
        ahou=0.0;
      }

      solrad=917.0*transcof[month]*(ahou*sin(rlatitude)*
             sin(declin)+cos(rlatitude)*cos(declin)*sin(ahou));

      /* Determine the short wave radiation outside the atmosphere */
      return(solrad/transcof[month]);
    }
