
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      initsite_tg.c
**
**  FUNCTION:  void initsite()
**
**  PURPOSE:   Read in site specific parameters (from sitepar.in)
**
**  AUTHOR:    Susan Chaffee    March 10, 1992
**
**  REWRITE:   Melannie Hartman  9/9/93 - 9/23/93
**
**  HISTORY:
**    8/13/92 (SLC) - Change the way the parameter for minimum soil water
**                    content is used.  No longer a function of wilting point,
**                    now it simply gives the percent of water at each layer.
**    5/08/08 (CAK) - Add checks for end of file, if the end of the file is
**                    reached before all of the variables have been read and
**                    initialized this will result in a fatal error.
**
**  INPUTS:
**    flags          - structure containing debugging flags
**    flags->debug   - flag to set debugging mode, 0 = off, 1 = on
**    flags->verbose - flag to set verbose debugging mode, 0 = off, 1 = on
**    sitename       - data file name containing site parameters
**    sitepar        - site specific parameters structure for soil water model
**    layers         - soil water soil layer structure
**
**  GLOBAL VARIABLES:
**    INPTSTRLEN - maximum length of input file line (120)
**    NTDEPTHS   - maximum number of soil regions (4)
**
**  LOCAL VARIABLES:
**    errmsg[]   - string containing error message
**    fp_in      - pointer to input file
**    imo        - current month (1..12)
**    inptline[] - line read from input file
**    m          - month
**
**  OUTPUTS:
**    Read from file sitename:
**    bioabsorp              - litter biomass at full absorption of radiation
**                             (grams biomass)
**    maxphoto               - maximum carbon loss due to photodecomposition
**                             (ug C/KJ srad)
**    sitepar->albedo        - fraction of light reflected by snow
**    sitepar->aspect        - site aspect (degrees)
**    sitepar->cldcov[]      - average cloud cover for the month (%, 1..100)
**    sitepar->dmp           - damping factor for calculating soil temperature
**                             by layer
**    sitepar->dmpflux       - damping factor for soil water flux (in h2oflux)
**    sitepar->drainlag      - number of days that soil drainage should lag
**                             behind rain/irrigation/melt events (1-5)
**    sitepar->ehoriz        - site east horizon, degrees
**    sitepar->elevation     - site elevation (meters)
**    sitepar->fswcinit      - initial soil water content, fraction of field
**                             capacity (0.0 - 1.0)
**    sitepar->hours_rain    - the duration of the rainfall/snowmelt event
**                             (hours)
**    sitepar->hpotdeep      - hydraulic water potential of deep storage
**                             layer, the more negative the number the dryer
**                             the soil layer (units?)
**    sitepar->jdayEnd       - the Julian day to end the turning off of the
**                             restriction of the CO2 effect on
**                             denitrification
**    sitepar->jdayStart     - the Julian day to start the turning off of the
**                             restriction of the CO2 effect on
**                             denitrification
**    sitepar->ksatdeep      - saturated hydraulic conductivity of deep
**                             storage layer (cm/sec)
**    sitepar->N2Oadjust     - nitrification N2O adjustment factor (0.0-1.0)
**    sitepar->Ncoeff        - minimum water/temperature limitation
**                             coefficient for nitrification
**    sitepar->reflec        - fraction of light reflected by vegetation
**    sitepar->slope         - site slope (degrees)
**    sitepar->sublimscale   - multiplier to scale sublimation
**    sitepar->tbotmn        - minimum emperature for bottom soil layer
**                             (degrees C)
**    sitepar->tbotmx        - maximum emperature for bottom soil layer
**                             (degrees C)
**    sitepar->texture       - texture classification for trace gas model
**                             (1 = coarse, 2 = medium, 3 = fine)
**    sitepar->timlag        - days from Jan 1 to coolest temp at bottom of
**                             soil (days)
**    sitepar->usexdrvrs     - 0 = use air temperature to drive PET rates
**                             1 = use extra weather drivers (solrad, rhumid,
**                                 windsp) for PET calculation
**                             2 = use extra weather drivers (srad, vpd)
**                                 for photosynthesis calculation
**                             3 = use extra drivers for both PET and
**                                 photosynthesis calculations
**    sitepar->watertable[]  - 1 = simulate water table, 0 = no water table
**    sitepar->whoriz        - site west horizon, degrees
**    sradadj[]              - solar radiation adjustment for cloud cover and
**                             transmission coeffient
**    tminintercept          - intercept used to adjust minimum temperature
**                             for calculating VPD at dewpoint
**    tminslope              - slope used to adjust minimum temperature
**                             for calculating VPD at dewpoint
**
**  EXAMPLE INPUT FILE:
**  0        / 0 = don't use extra weather drivers for PET 
**             1 = use extra weather drivers for PET (solrad, rhumid, windsp)
**             2 = use extra weather drivers (solrad, vpd) for photosynthesis
**             3 = use extra drivers for both PET and photosynthesis
**  1.0      / sublimscale
**  0.18     / reflec - vegetation reflectivity/albedo (frac)
**  0.65     / albedo - snow albedo (frac)
**  0.90     / fswcinit - initial swc, fraction of field capacity
**  0.000001 / dmpflux - in h2oflux routine (0.000001 = original value)
**  4        / hours_rain - duration of each rain event
**  0        / # of days between rainfall event and drainage of soil (-1=computed)
**  1  0     / watertable[month] - 0 = no water table, 1 = water table
**  2  0
**  3  0
**  4  0
**  5  0
**  6  0
**  7  0
**  8  0
**  9  0
**  10 0
**  11 0
**  12 0
**  -200     / hpotdeep - hydraulic water potential of deep storage layer
**             (units?)
**  0.0002   / ksatdeep - saturated hydraulic conductivity of deep storage
**             layer (cm/sec)
**  1  58    / cldcov[month] - cloud cover (%)
**  2  58
**  3  58
**  4  58
**  5  58
**  6  58
**  7  58
**  8  58
**  9  58
**  10 58
**  11 58
**  12 58
**  5.0 16.4 / min and max temperature for bottom soil layer (degrees C)
**  0.003    / damping factor for calculating soil temperature by layer
**  30.0     / timlag, days from Jan 1 to coolest temp at bottom of soil (days)
**  0.03     / min water/temperature limitation coefficient for nitrify
**  0 0      / turn off resp restraint on denit between these days of the year
**  0.8      / nitrification N2O adjustment factor (0.0-1.0)
**  1525     / elevation, meters
**  0.0      / site slope, degrees
**  0.0      / site aspect, degrees
**  0.0      / site east horizon, degrees
**  0.0      / site west horizon, degrees
**  1   0.63 / sradadj[month] srad adjust for cloud cover & transmission coeff
**  2   0.63 
**  3   0.63 
**  4   0.63 
**  5   0.63 
**  6   0.63 
**  7   0.63 
**  8   0.63 
**  9   0.63 
**  10  0.63 
**  11  0.63
**  12  0.63
**  1.0      / slope for adjusting minimum temperature for VPD dewpoint calc
**  0.0      / intercept for adjusting minimum temperature for VPD dewpoint calc
**  1.16     / maximum carbon loss due to photodecomposition (ug C/KJ srad)
**  200.0    / litter biomass for full absorption of solar radiation (g biomass)
**
**  CALLED BY:
**    initsw()
**
**  CALLS:
**    None
**
*****************************************************************************/

#include "soilwater.h"
#include <stdlib.h>
#include <stdio.h>

    void initsite(char *sitename, SITEPAR_SPT sitepar, LAYERPAR_SPT layers,
                  FLAG_SPT flags, double sradadj[NMONTH], double *tminslope,
                  double *tminintercept, double *maxphoto, double *bioabsorp)
    {

      int  imo, m;
      char inptline[INPTSTRLEN];
      char errmsg[INPTSTRLEN];
      FILE *fp_in;

      if (flags->debug) {
        printf("Entering function initsite\n");
      }

      if ((fp_in = fopen(sitename, "r")) == NULL) {
        sprintf(errmsg, "Cannot open file %s\n", sitename);
        perror(errmsg);
        exit(1);
      }

      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%d", &sitepar->usexdrvrs);
        if (flags->debug) {
          printf("usexdrvrs: %1d\n", sitepar->usexdrvrs);
        }
      } else {
        printf("ERROR:  Missing extra weather drivers value\n");
        printf("in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%lf", &sitepar->sublimscale);
        if (flags->debug) {
          printf("sublimscale: %lf\n", sitepar->sublimscale);
        }
      } else {
        printf("ERROR:  Missing sublimation scalar value\n");
        printf("in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%lf", &sitepar->reflec);
        if (flags->debug) {
          printf("reflec: %lf\n", sitepar->reflec);
        }
      } else {
        printf("ERROR:  Missing fraction of light reflected by vegetation\n");
        printf("value in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%lf", &sitepar->albedo);
        if (flags->debug) {
          printf("albedo: %lf\n", sitepar->albedo);
        }
      } else {
        printf("ERROR:  Missing fraction of light reflected by\n");
        printf("snow value in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%lf", &sitepar->fswcinit);
        if (flags->debug) {
          printf("fswcinit: %lf\n", sitepar->fswcinit);
        }
      } else {
        printf("ERROR:  Missing initial soil water content value\n");
        printf("in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%lf", &sitepar->dmpflux);
        if (flags->debug) {
          printf("dmpflux: %lf\n", sitepar->dmpflux);
        }
      } else {
        printf("ERROR:  Missing damping factor for soil water flux value\n");
        printf("in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%lf", &sitepar->hours_rain);
        if (flags->debug) {
          printf("hours_rain: %lf\n", sitepar->hours_rain);
        }
      } else {
        printf("ERROR:  Missing duration of the rainfall/snowmelt event\n");
        printf("value in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      /* Allow user to set number of days between rainfall event and */
      /* drainage of soil profile.  If a value of -1 is entered set the */
      /* number of days to drainage based in the soil texture.  Constrain */
      /* the number of days to drainage to be <=5 to prevent numerical */
      /* instabilities in the h2oflux subroutine.  cak - 02/13/04 */
      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%d", &sitepar->drainlag);
        if (sitepar->drainlag < 0) {
          sitepar->drainlag = sitepar->texture - 1;
        }
        if (sitepar->drainlag > 5) {
          printf("lag period for drainage too long, setting to max value\n");
          sitepar->drainlag = 5;
        }
        if (flags->debug) {
          printf("drainlag: %d\n", sitepar->drainlag);
        }
      } else {
        printf("ERROR:  Missing number of days that soil drainage should\n");
        printf("lag behind rain/irrigation/melt events value\n");
        printf("in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      for(imo=1; imo<=12; imo++) {
        if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
          sscanf(inptline, "%d %d", &m, &sitepar->watertable[imo]);
          if (flags->debug) {
            printf("watertable: %d  %d\n", imo, sitepar->watertable[imo]);
          }
        } else {
          printf("ERROR:  Missing simulate water table flag\n");
          printf("value in sitepar.in file.  Ending simulation!\n");
          exit(1);
        }
      }

      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%lf", &sitepar->hpotdeep);
        if (flags->debug) {
          printf("hpotdeep: %lf\n", sitepar->hpotdeep);
        }
      } else {
        printf("ERROR:  Missing hydraulic water potential of deep storage\n");
        printf("layer value in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%lf", &sitepar->ksatdeep);
        if (flags->debug) {
          printf("ksatdeep: %lf\n", sitepar->ksatdeep);
          printf("Cloud cover:\n");
        }
      } else {
        printf("ERROR:  Missing saturated hydraulic conductivity of deep\n");
        printf("storage layer value in sitepar.in file.\n");
        printf("Ending simulation!\n");
        exit(1);
      }

      for(imo=1; imo<=12; imo++) {
        if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
          sscanf(inptline, "%d %lf", &m, &sitepar->cldcov[imo]);
          if (flags->debug) {
            printf("cloud cover: %d  %6.2lf\n", imo, sitepar->cldcov[imo]);
          }
        } else {
          printf("ERROR:  Missing average cloud cover value\n");
          printf("in sitepar.in file.  Ending simulation!\n");
          exit(1);
        }
      }

      /* The texture parameter is being set in the initlyrs subroutine */
      /* CAK - 05/31/01 */
      /* The texture parameter has been replaced by the minimum and */
      /* maximum temperature values for the bottom soil layer */
      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
/*      sscanf(inptline, "%d", &sitepar->texture); */
        sscanf(inptline, "%lf %lf", &sitepar->tbotmn, &sitepar->tbotmx);
        if (flags->debug) {
          printf("tbotmn: %lf\n", sitepar->tbotmn);
          printf("tbotmx: %lf\n", sitepar->tbotmx);
        }
        if (sitepar->tbotmn > sitepar->tbotmx) {
          fprintf(stderr, "Error in input file %s.\n", sitename);
          fprintf(stderr, "the minimum soil temperature at the bottom of\n");
          fprintf(stderr, "the soil is greater than the maximum soil\n");
          fprintf(stderr, "temperature at the bottom of the soil.\n");
          exit(1);
        }
      } else {
        printf("ERROR:  Missing minimum and maximum temperature for\n");
        printf("bottom soil layer values in sitepar.in file.\n");
        printf("Ending simulation!\n");
        exit(1);
      }

      /* Added dmp parameter to the sitepar.in file, read in the value */
      /* for this parameter, cak - 12/16/02 */
      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%lf", &sitepar->dmp);
        if (flags->debug) {
          printf("damping factor: %lf\n", sitepar->dmp);
        }
      } else {
        printf("ERROR:  Missing damping factor for calculating soil\n");
        printf("temperature value in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      /* Added timlag parameter to the sitepar.in file, read in the value */
      /* for this parameter, cak - 04/24/03 */
      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%lf", &sitepar->timlag);
        if (flags->debug) {
          printf("lag time: %lf\n", sitepar->timlag);
        }
      } else {
        printf("ERROR:  Missing days from Jan 1 to coolest temp at bottom\n");
        printf("of soil value in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      /* Added Ncoeff parameter to the sitepar.in file, read in the value */
      /* for this parameter, cak - 04/08/03 */
      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%lf", &sitepar->Ncoeff);
        if (flags->debug) {
          printf("water/temp coeff: %lf\n", sitepar->Ncoeff);
        }
      } else {
        printf("ERROR:  Missing coefficient for nitrification value\n");
        printf("in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%d %d", &sitepar->jdayStart, &sitepar->jdayEnd);
        if (flags->debug) {
          printf("jdayStart: %d\n", sitepar->jdayStart);
          printf("jdayEnd: %d\n", sitepar->jdayEnd);
        }
      } else {
        printf("ERROR:  Missing days of year to turn off respiration\n");
        printf("restraint on denitrification value in sitepar.in file.\n");
        printf("Ending simulation!\n");
        exit(1);
      }

      /* Added nitrification N2O adjustment factor to sitepar.in file, */
      /* cak - 05/15/2007 */
      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%lf", &sitepar->N2Oadjust);
        if (flags->debug) {
          printf("N2O adjustment factor: %lf\n", sitepar->N2Oadjust);
        }
      } else {
        printf("ERROR:  Missing N2O adjustment factor value\n");
        printf("in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      /* Added elevation (meters) to sitepar.in file, */
      /* cak - 04/15/2009 */
      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%lf", &sitepar->elevation);
        if (flags->debug) {
          printf("Elevation (meters): %lf\n", sitepar->elevation);
        }
      } else {
        printf("ERROR:  Missing elevation (meters) value\n");
        printf("in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      /* Added site slope, aspect, east horizon, and west horizon */
      /* to sitepar.in file, */
      /* cak - 05/08/2009 */
      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%lf", &sitepar->slope);
        if (flags->debug) {
          printf("Site slope (degrees): %lf\n", sitepar->slope);
        }
      } else {
        printf("ERROR:  Missing site slope (degrees) value\n");
        printf("in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%lf", &sitepar->aspect);
        if (flags->debug) {
          printf("Site aspect (degrees): %lf\n", sitepar->aspect);
        }
      } else {
        printf("ERROR:  Missing site aspect (degrees) value\n");
        printf("in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%lf", &sitepar->ehoriz);
        if (flags->debug) {
          printf("Site east horizon (degrees): %lf\n", sitepar->ehoriz);
        }
      } else {
        printf("ERROR:  Missing site east horizon (degrees) value\n");
        printf("in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%lf", &sitepar->whoriz);
        if (flags->debug) {
          printf("Site west horizon (degrees): %lf\n", sitepar->whoriz);
        }
      } else {
        printf("ERROR:  Missing site west horizon (degrees) value\n");
        printf("in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      for(imo=1; imo<=12; imo++) {
        if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
          sscanf(inptline, "%d  %lf", &m, &sradadj[imo-1]);
          if (flags->debug) {
            printf("Solar radiation adjustment factor: %d %lf\n", imo,
                   sradadj[imo-1]);
          }
        } else {
          printf("ERROR:  Missing solar radiation adjustment factor\n");
          printf("value in sitepar.in file.  Ending simulation!\n");
          exit(1);
        }
      }

      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%lf", tminslope);
        if (flags->debug) {
          printf("Slope for VPD min temp adjustment: %lf\n", *tminslope);
        }
      } else {
        printf("ERROR:  Slope for VPD minimum temperature adjustment\n");
        printf("value in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%lf", tminintercept);
        if (flags->debug) {
          printf("Slope for VPD min temp adjustment: %lf\n", *tminintercept);
        }
      } else {
        printf("ERROR:  Intercept for VPD minimum temperature adjustment\n");
        printf("value in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%lf", maxphoto);
        /* Convert ug C to grams C */
        *maxphoto *= 0.000001;
        if (flags->debug) {
          printf("Maximum carbon loss due to photodecomp: %lf\n", *maxphoto);
        }
      } else {
        printf("ERROR:  Missing maximum carbon loss due to photodecomp\n");
        printf("value in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%lf", bioabsorp);
        if (flags->debug) {
          printf("Litter biomass at full absorption of radiation: %lf\n",
                 *bioabsorp);
        }
      } else {
        printf("ERROR:  Missing litter biomass at full absorption of srad\n");
        printf("value in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        printf("ERROR:  Incorrect number of lines in sitepar.in file.\n");
        printf("Ending simulation!\n");
        exit(1);
      }

      fclose(fp_in);

      if (flags->debug) {
        printf("Exiting function initsite\n");
      }

      return;
    }
