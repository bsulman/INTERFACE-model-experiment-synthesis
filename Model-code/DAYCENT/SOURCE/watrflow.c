
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:  watrflow.c
**
**  FUNCTION:  void watrflow()
**
**  PURPOSE: Water-flow submodel.  This submodel is a rewrite of a
**           model originally written by William Parton.  It simulates
**           the flow of water through the plant canopy and soil.
**           See "Abiotic Section of ELM", as a reference.
** 
**  AUTHOR:  Susan Chaffee    4/30/92 
** 
**  REWRITE:  Melannie Hartman  9/20/93 - 8/21/96
**
**  HISTORY:
**
**  INPUTS:
**    avgtemp    - average air temperature for the day (deg C)
**    nlayer     - number of layers in Century soil profile
**    pet        - daily potential evapotranspiration rates (cm H2O)
**    pptactual  - the current day's precipitation (cm H2O)
**    tempmax    - maximum air temperature for the day (deg C)
**    tempmin    - minimum air temperature for the day (deg C)
**
**  GLOBAL VARIABLES:
**    CENTMAXLYR - maximum number of Century soil layers (10)
**    MAXLYR     - maximum number of soil water model layers (21)
**
**  EXTERNAL VARIABLES:
**    layers                - soil water soil layer structure
**    layers->depth[]       - the distance from the surface to the middle of
**                            the soil layer (cm)
**    layers->fieldc[]      - volumetric water content at field capacity for
**                            layer (cm H2O/cm of soil)
**    layers->lbnd[]        - the index of the lower soil water model layer
**                            which corresponds to given soil layer in Century
**    layers->numlyrs       - total number of layers in the soil water model
**                            soil profile
**    layers->swc[]         - soil water content by layer (cm H2O)
**    layers->swcfc[]       - volumetric soil water content at field capacity
**                            for layer (cm H2O/cm of soil)
**    layers->width[]       - the thickness of soil water model layers (cm)
**
**    soil                  - soil temperature structure
**    soil->soiltavg[]      - average soil temperature of layer (degrees C)
**    soil->soiltmax[]      - maximum soil temperature by layer (degrees C)
**    soil->soiltmin[]      - minimum soil temperature by layer (degrees C)
**    soil->stmtemp[]       - soil surface temperature (degrees Celsius)
**
**  LOCAL VARIABLES:
**    pptsoil    - the amount of water entering the soil including rain and melt (cm)
**    swctemp[]  - temporary storage location for layers->swc[] array values
**
**  OUTPUTS:
**    amovdly    - the amount of water moving through the bottom of a soil
**                 layer (cm H2O) (positive is downward, negative is upward)
**    asmos[]    - soil water content by layer (cm H2O)
**    avh2o[]    - water available for plant growth (avh2o[0]), plant survival
**                 (avh2o[1]), and in the first two Century soil layers 
**                 (avh2o[2])
**    rwcf[]     - relative water content by layer
**    wfluxout[] - the amount of water moving through the bottom of a soil
**                 layer (cm H2O) (positive is downward, negative is upward)
**
**  CALLED BY:
**    dailymoist()
**
**  CALLS:
**    setasmos()     - set asmos, avh2o and rfwc (Century variables) from swc,
**                     the soil water content in the daily soil water model
**
*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include "soilwater.h"

void watrflow(int *nlayer, 
    double *avgtemp, double *tempmin, double *tempmax, double *pptactual, double *pet, 
    double wfluxout[], double amovdly[CENTMAXLYR], double rwcf[CENTMAXLYR], 
    double asmos[CENTMAXLYR], double avh2o[3], double *rprpet)
    {
        extern LAYERPAR_SPT layers;
        extern SOIL_SPT soil;

        int ilyr;
        double pptsoil;
        double swctemp[MAXLYR];

        /* Initilize function arguments that are updated by this function */

        /* Soil temperature variables */
        for(ilyr=0; ilyr < layers->numlyrs; ilyr++) 
        {
            /*soil->soiltavg[ilyr] = *avgtemp;*/
            /*soil->soiltmax[ilyr] = *tempmax;*/
            /*soil->soiltmin[ilyr] = max(*tempmin, -1.0);*/
            soil->soiltavg[ilyr] = 30.0;
            soil->soiltmax[ilyr] = 30.0;
            soil->soiltmin[ilyr] = 30.0;
        }

        for(ilyr=0; ilyr < MAXSTLYR; ilyr++)
        {
            soil->stmtemp[ilyr] = 30;
        }

        /* soil moisture and soil water flux variables */

	/* Century soil layers - this are updated below */
        for(ilyr=0; ilyr <= layers->numlyrs; ilyr++) 
        {
            swctemp[ilyr] = layers->swc[ilyr];
        }
        setasmos(asmos, nlayer, swctemp, &layers->numlyrs, avh2o, rwcf);

        /* daily soil water layers */
        for(ilyr=0; ilyr < layers->numlyrs; ilyr++) 
        {
            layers->swc[ilyr] = layers->fieldc[ilyr] * layers->width[ilyr];
            wfluxout[ilyr] = 0.0;
        }
        for(ilyr=0; ilyr < *nlayer; ilyr++) 
        {
            amovdly[ilyr] = 0.0;
        }

        pptsoil = *pptactual;
        *rprpet = (pptsoil + avh2o[2])/ *pet;

        return;
    }
