
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      calcdefac.c
**
**  FUNCTION:  void calcdefac()
**
**  PURPOSE:   This routine calculates the decomposition factor based on 
**             water and temperature, the temperature effect on decomposition,
**             the water effect on decomposition, and the average water filled
**             pore space in the top 15 centimeters of the soil profile.
**
**  INPUTS:
**    idef    - flag used to determine which of three water curves options to
**              use for calculating the water effect on decompostion, read
**              from fix.100 file
**                idef = 1, use relative water content
**                idef = 2, use ratio of precipitation to potential
**                          evapotranspiration
**                idef = 3, use water filled pore space
**    ppt     - daily precip (cm)
**    rprpet  - ratio of precipitation to PET
**    snow    - snow cover (cm SWE)
**    teff[]  - coefficients for temperature function read from fix.100 file
**    texture - a constant to designate coarses/fineness of the soil
**              (i.e. COARSE, MEDIUM, FINE, VERYFINE - see n2o_model.h)
**
**  GLOBAL VARIABLES:
**    COARSE          - designates a coarse, sandy soil, texture (1)
**    FINE            - designates a fine soil texture (3)
**    MEDIUM          - designates a medium, loamy soil, texture (2)
**    VERYFINE        - designates a very fine, volcanic soil, texture (4)
**
**  EXTERNAL VARIABLES:
**    layers           - soil water soil layer structure
**    layers->wfps[]   - water-filled pore space by layer (fraction 0.0-1.0)
**                       (fraction of a porespace that is filled with water)
**    layers->width[]  - the thickness of soil water model layers (cm)
**    soil             - soil temperature structure
**    soil->soiltavg[] - average soil temperature of layer (degrees C)
**
**  LOCAL VARIABLES:
**    a, b, c, d - intermediate variable for calculations
**    agwfunc    - water effect on surface decomposition
**    avg_rel_wc - average relative soil water content
**    base1      - intermediate base variable for calculations
**    base2      - intermediate base variable for calculations
**    e1, e2     - intermediate exponent variables for calculations
**    ilyr       - index
**    krainwfunc - increase of wfunc due to moisture and rain >= 1.0
**    rel_wc[]   - relative soil water content by layer
**    avgstemp - weighted average of the average soil temperature in the
**               second and third soil layer used when calculating tfunc
**               (degrees C)
**
**  OUTPUTS:
**    agdefac  - decomposition factor based on water and temperature for
**               surface decomposition
**    avgwfps  - average wfps in top 15 cm (0-1)
**    bgdefac  - decomposition factor based on water and temperature for
**               soil decomposition
**    bgwfunc  - water effect on soil decomposition
**    tfunc    - temperature effect on decomposition
**
**  CALLED BY:
**    dailymoist()
**
**  CALLS:
**    tcalc         - computes the effect of temperature on decomposition
**    wfunc_pulse() - increase in wfunc due to moisture and rain
**    wfps()        - compute the daily water-filled pore space by layer
**
*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "soilwater.h"
#include "n2o_model.h"

/* Prototype to allow C to call a Fortran function */
extern void tcalc(double *, double [], double *);

    void calcdefac(int *texture, double *tfunc, double *agwfunc, double *bgwfunc, 
                   double *agdefac, double *bgdefac, double *avgwfps, double teff[],
                   double *rprpet, int *idef, double *ppt, double *snow, int *experiment)
    {
      int   ilyr;
      double a, b, c, d;
      double base1, base2;
      double e1, e2;
      double rel_wc[4], avg_rel_wc;
      double krainwfunc;
      double avgstemp;

      extern LAYERPAR_SPT layers;
      extern SOIL_SPT soil;

      if (*experiment == 1) 
      {
          *agdefac = 1.0;
          *bgdefac = 1.0;
          return;
      }
     
      /* avgwfps is passed to the trace_gas_model, calculate this value */
      /* every time calcdefac is called, cak - 09/22/03 */
      wfps(layers);
      *avgwfps = (layers->wfps[0]*layers->width[0] +
                 layers->wfps[1]*layers->width[1] +
                 layers->wfps[2]*layers->width[2]) /
                 (layers->width[0] + layers->width[1] + layers->width[2]);

      switch (*idef) {
        case 1:
          /* Compute water effect for surface decomposition using the */
          /* top soil layer, cak - 04/01/04 */
          rel_wc[0] = ((double)layers->swc[0]/(layers->width[0]) -
                       layers->swclimit[0]) /
                       (layers->fieldc[0] - layers->swclimit[0]);
          if (rel_wc[0] > 1.0) {
            *agwfunc = 1.0;
          } else {
            if (rel_wc[0] < 0.0) {
              rel_wc[0] = 0.0;
            }
            *agwfunc = 1.0/(1.0 + 30.0 * (double)exp(-9.0 * rel_wc[0]));
          }
          /* Compute water effect for soil decomposition using a weighted */
          /* averaging of the 2nd, and 3rd soil layers, cak - 08/01/04 */
          for (ilyr = 1; ilyr < 3; ilyr ++) {
            rel_wc[ilyr] = ((double)layers->swc[ilyr]/(layers->width[ilyr]) -
                            layers->swclimit[ilyr]) /
                            (layers->fieldc[ilyr] - layers->swclimit[ilyr]);
            if (rel_wc[ilyr] < 0.0) {
              rel_wc[ilyr] = 0.0;
            } else if (rel_wc[ilyr] > 1.0) {
              rel_wc[ilyr] = 1.0;
            }
            rel_wc[ilyr] *= layers->width[ilyr];
          }
          avg_rel_wc = (rel_wc[1] + rel_wc[2]) / 
                       (layers->width[1] + layers->width[2]);
          if (avg_rel_wc > 1.0) {
            *bgwfunc = 1.0;
          } else {
            *bgwfunc = 1.0/(1.0 + 30.0 * exp(-9.0 * avg_rel_wc));
          }
          break;
        case 2:
          if (*rprpet > 9) {
            *agwfunc = 1.0;
          } else {
            *agwfunc = 1.0/(1.0 + 30.0 * exp(-8.5 * *rprpet));
          }
          *bgwfunc = *agwfunc;
          break;
        case 3:
          /* Note:  TEXTURE is set in the initlyrs subroutine using a weighted */
          /*        average of the sand content in the top 3 soil layers to */
          /*        match this calculation for avgwfps.  If the calculation */
          /*        for avgwfps is changed make an approptiate change to the */
          /*        TEXTURE calculation in the initlyrs subroutine. */
          /* CAK - 05/31/01) */
          switch (*texture) {
            case COARSE:
              a = 0.55;
              b = 1.70;
              c = -0.007;
              d = 3.22;
              break;
            case MEDIUM:
              a = 0.60;
              b = 1.27;
              c = 0.0012;
              d = 2.84;
              break;
            case FINE:
            case VERYFINE:
              a = 0.60;
              b = 1.27;
              c = 0.0012;
              d = 2.84;
              break;
            default:
              printf("Error in calcdefac, unknown texture = %1d\n", *texture);
              exit(1);
          }
          /* Compute the soil moisture effect on decomposition */
          base1 =((*avgwfps-b) / (a-b));
          base2 =((*avgwfps-c) / (a-c));
          e1 = d * ((b-a)/(a-c));
          e2 = d;
          *agwfunc = (pow(base1, e1)) * (pow(base2,e2));
          *bgwfunc = *agwfunc;
          break;
        default:
          printf("Error in calcdefac, unknown idef value = %1d\n", *idef);
          exit(1);
      }

      /* Compute the soil temperature effect on decomposition */
      /* Compute the temperature effect on decomposition */
      /* using an arctangent function.  CAK - 03/16/01   */
      /* Use a weighted average of the average soil temperature */
      /* in the second and third soil layer when calculating tfunc, */
      /* cak - 07/30/04 */
      avgstemp = (soil->soiltavg[1] * layers->width[1] + 
                  soil->soiltavg[2] * layers->width[2]) /
                 (layers->width[1] + layers->width[2]);
      tcalc(&avgstemp, teff, tfunc);

      /* On days when it rains the surface decomposition is elevated, */
      /* cak - 04/01/04 */
      if (*ppt > 0.1) {
        *agdefac = max(0.000001, *tfunc * 3.0);
      } else {
        *agdefac = max(0.000001, *tfunc * *agwfunc);
      }
      /* Pulse multiplier for enhanced soil decomposition following */
      /* drying and re-wetting of the soil, cak - 08/25/04 */
      if ( *experiment == 0) 
      {
        krainwfunc = wfunc_pulse(ppt, snow);
      }
      else
      {
        krainwfunc = 1.0;
      }

      *bgdefac = max(0.000001, *tfunc * *bgwfunc*krainwfunc);

      return;
    }
