
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      petrad.c
**
**  FUNCTION:  double petrad()
**
**  PURPOSE:   Calculate the potential evapotranspiration rate using
**             Penmans equation (1948).
**
**  REWRITE:  Melannie Hartman  9/21/93 - 10/6/93 
**
**  NOTES:    Mathematical function substitution:
**              atan2(y,x) is unchanged - arctan(y/x) in the range [-PI, PI]
**
**  HISTORY:
**    4/30/92  (SLC)
**
**  INPUTS: 
**    albedo    - fraction of light reflected by snow
**    avgtemp   - average air temperature for the day (deg C)
**    cldcov    - average cloud cover for the month (%, 1..100)
**    jday      - current julian day (1..366)
**    month     - current month of the year (1..12)
**    reflec    - fraction of light reflected by vegetation
**    rhumid    - average relative humidity for the day (% 1..100)
**    rlatitude - latitude of the site (in radians)
**    snowpack  - current snowpack (equiv. cm H2O)
**    solrad    - total incoming shortwave radiation (langleys/day)
**    windsp    - average daily windspeed at 2 meters (mph)
**
**  GLOBAL VARIABLES:
**    None
**
**  LOCAL VARIABLES:
**    ahou        - intermediate variable for calculations
**    arads       - intermediate variable for calculations
**    clrsky      - relative amount of clear sky (decimal fraction)
**    declin      - declination (radians)
**    fhumid      - vapor pressure (mm Hg)
**    ftemp       - intermediate variable for calculations
**    kelvin      - avgtemp in kelvin degrees
**    KELVIN_CONV - conversion factor for converting degrees C to degrees K
**    netreflec   - net reflectivity (use either reflec or albedo)
**    par1, par2  - parameters in computation of pet
**    shwave      - short wave solar radiation (mm H2O)
**
**  OUTPUTS:
**    pet - potential evapotranspiration rate (cm/day)
**
**  CALLED BY:
**    calcpet()
**
**  CALLS:
**    svapor() - calculate the saturation vapor pressure of water for air
**               temperature
**
*****************************************************************************/

#include <math.h>
#include "soilwater.h"

#define KELVIN_CONV 273

    double petrad(int jday, int month, double cldcov, double windsp,
                 double rhumid, double rlatitude, double avgtemp, double reflec,
                 double solrad, double albedo, double snowpack)
    {
      double pet;      /* return value */
      double declin, ahou;
      double shwave,kelvin,arads;
      double clrsky,fhumid,ftemp;
      double par1, par2;
      double netreflec;

      /* calculate the short wave solar radiation */
      /* on a clear day using a equation presented by sellers(1965) */

      /* if there is snow on the ground, use snow albedo as */
      /* reflectivity, else use vegetation albedo.  Suggested */
      /* by Bill Parton 9/15/94 */

      netreflec = reflec;
      if (jday > 1) {
        if (snowpack > 0) {
          netreflec = albedo;
        }
      }

      declin = 0.401426*sin(2*PI*(jday-77)/365);

      /* Original FORTRAN statement for par1:             */
      /*   par1=sqrt(1.-(-tan(rlatitude)*tan(declin))**2) */

      par1 = (1-pow((-tan(rlatitude)*tan(declin)),2));
      par1 = sqrt(par1);
      par2 = (-tan(rlatitude)*tan(declin));

      ahou = atan2(par1,par2);

      if(ahou < 0) {
        ahou = 0.0;
      }

      /* shwave = mm of water equivalent for total incoming radiation */
      /* 0.0168 converts langleys/day to mm H2O (-mdh 9/94) */
      /* 10*(1/597) = 0.0168: mult by 10 to get mm, divide by latent heat */
      /* of vaporization (cal/cm3) at 0 deg C */

      shwave = solrad*0.0168;

      kelvin = avgtemp + KELVIN_CONV;
      arads = svapor(avgtemp) * (3010.21)/(kelvin*kelvin);

      clrsky = 1 - cldcov/100;
      fhumid = rhumid * svapor(avgtemp)/100;

      ftemp = (avgtemp + KELVIN_CONV)*0.01;
      ftemp = pow(ftemp, 4.0) * (0.201);

      par1 = (0.35)*(svapor(avgtemp)-fhumid) * (1+0.0098*(windsp*24));

      par2=shwave*(1-netreflec)-ftemp*(0.56-0.092*
           sqrt(fhumid))*(0.10+0.90*clrsky);

      pet = ((arads*par2 + 0.27*par1)/(arads + 0.27))/10;

      if(pet < 0.03) {
        pet = 0.03;
      }

      return(pet);
    }
