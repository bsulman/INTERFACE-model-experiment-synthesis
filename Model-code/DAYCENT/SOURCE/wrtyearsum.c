/*****************************************************************************
**
**  FILE:      wrtyearsum.c
**
**  FUNCTION:  void wrtyearsum()
**
**  PURPOSE:   This function writes out the yearly summation of N2O flux, NO
**             flux, and CH4.
**
**  INPUTS:
**    annppt       - Annual precipitation and irrigation (cm)
**    CH4_year     - Methane oxidation for year (gCH4/m^2)
**    N2_year      - Nitrogen gas for year (gN/m^2)
**    N2O_year     - Nitrous oxide flux for year (gN/m^2)
**    nit_amt_year - Gross nitrification for year ((gN/m^2)
**    NO_year      - Nitric oxide flux for year (gN/m^2)
**    time         - current simulation time (years)
**
**  GLOBAL VARIABLES:
**    None
**
**  EXTERNAL VARIABLES:
**    files->fp_yearsum    - file pointer to co2.out output file
**    files->write_yearsum - flag to indicate if co2.out output file should
**                           be created, 0 = do not create, 1 = create
**
**  LOCAL VARIABLES:
**    None
**
**  OUTPUTS:
**    None
**
**  CALLED BY:
**    main()
**
*****************************************************************************/

#include <math.h>
#include <stdio.h>
#include "soilwater.h"

    void wrtyearsum(double *time, double *N2O_year, double *NO_year,
                    double *N2_year, double *CH4_year, double *nit_amt_year,
                    double *annppt)
    {

      extern FILES_SPT files;

      if (!files->write_yearsum) {
        return;
      }

      fprintf(files->fp_yearsum, "%8.2lf,%12.6lf,%12.6lf,%12.6lf,%12.6lf,",
              *time, *N2O_year, *NO_year, *N2_year, *CH4_year);
      fprintf(files->fp_yearsum, "%12.6lf,%12.6lf\n", *nit_amt_year, *annppt);

      return;
    }
