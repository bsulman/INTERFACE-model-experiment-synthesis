
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      closefiles.c
**
**  FUNCTION:  void closefiles
**
**  PURPOSE:   Close any daily output files that were opened by the model
**              prior to ending the simulation.
**
**  AUTHOR:    Cindy Keough 05/2010
**
**  INPUTS:
**
**  EXTERNAL VARIABLES:
**    files->fp_bio         - file pointer to biowk.out output file
**    files->fp_cflows      - file pointer to cflows.out output file
**    files->fp_co2         - file pointer to co2.out output file
**    files->fp_dcsip       - file pointer to dc_sip.csv output file
**    files->fp_deadc       - file pointer to deadcwk.out output file
**    files->fp_dels        - file pointer to dels.out output file
**    files->fp_dN2lyr      - file pointer to dN2lyr.out output file
**    files->fp_dN2Olyr     - file pointer to the dN2Olyr.out output file
**    files->fp_gresp       - file pointer to gresp.out output file
**    files->fp_harv        - file pointer to the harvest.csv output file
**    files->fp_livec       - file pointer to livecwk.out output file
**    files->fp_mresp       - file pointer to mresp.out output file
**    files->fp_psyn        - file pointer to psyn.out output file
**    files->fp_soilc       - file pointer to soilcwk.out output file
**    files->fp_soiln       - file pointer to soiln.out output file
**    files->fp_soiltavg    - file pointer to soiltavg.out output file
**    files->fp_soiltmax    - file pointer to soiltmax.out output file
**    files->fp_soiltmin    - file pointer to soiltmin.out output file
**    files->fp_stempdx     - file pointer to stemp_dx.out output file
**    files->fp_swc         - file pointer to vswc.out output file
**    files->fp_sysc        - file pointer to syscwk.out output file
**    files->fp_tgmonth     - file pointer to the tgmonth.out output file
**    files->fp_wb          - file pointer to watrbal.out output file
**    files->fp_wflux       - file pointer to wflux.out output file
**    files->fp_wfps        - file pointer to wfps.out output file
**    files->fp_yearsum     - file pointer to year_summary.out output file
**    files->fp_yrcflows    - file pointer to year_cflows.out output file
**    files->write_bio      - flag to indicate if bio.out output file should
**                            be created, 0 = do not create, 1 = create
**    files->write_cflows   - flag to indicate if cflows.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_co2      - flag to indicate if co2.out output file should
**                            be created, 0 = do not create, 1 = create
**    files->write_dcsip    - flag to indicate if dc_sip.csv output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_deadc    - flag to indicate if deadc.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_dels     - flag to indicate if dels.out output file should
**                            be created, 0 = do not create, 1 = create
**    files->write_dN2lyr   - flag to indicate if dN2lyr.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_dN2Olyr  - flag to indicate if dN2Olyr.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_gresp    - flag to indicate if gresp.out output file should
**                            be created, 0 = do not create, 1 = create
**    files->write_harvest  - flag to indicate if harvest.csv output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_livec    - flag to indicate if livec.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_mresp    - flag to indicate if mresp.out output file should
**                            be created, 0 = do not create, 1 = create
**    files->write_psyn     - flag to indicate if psyn.out output file should
**                            be created, 0 = do not create, 1 = create
**    files->write_soilc    - flag to indicate if soilc.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_soiln    - flag to indicate if  output file should
**                            be created, 0 = do not create, 1 = create
**    files->write_soiltavg - flag to indicate if soiltavg.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_soiltmax - flag to indicate if soiltmax.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_soiltmin - flag to indicate if soiltmin.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_stempdx  - flag to indicate if stemp_dx.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_swc      - flag to indicate if vswc.out output file should
**                            be created, 0 = do not create, 1 = create
**    files->write_sysc     - flag to indicate if sysc.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_wb       - flag to indicate if watrbal.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_wflux    - flag to indicate if wflux.out output file should
**                            be created, 0 = do not create, 1 = create
**    files->write_wfps     - flag to indicate if wfps.out output file should
**                            be created, 0 = do not create, 1 = create
**    files->write_yearsum  - flag to indicate if year_summary.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_yrcflows - flag to indicate if year_cflows.out output file
**                            should be created, 0 = do not create, 1 = create
**
**  LOCAL VARIABLES:
**    None
**
**  OUTPUTS:
**    None
**
**  CALLED BY:
**    csa_main()
**
**  CALLS:
**    None
**
*****************************************************************************/

#include "soilwater.h"

    void closefiles()
    {
      extern FILES_SPT files;

      if (files->write_soiln) {
        fclose(files->fp_soiln);
      }

      if (files->write_co2) {
        fclose(files->fp_co2);
      }

      if (files->write_yearsum) {
        fclose(files->fp_yearsum);
      }

      if (files->write_soilc) {
        fclose(files->fp_soilc);
      }

      if (files->write_tgmonth) {
        fclose(files->fp_tgmonth);
      }

      if (files->write_dN2lyr) {
        fclose(files->fp_dN2lyr);
      }

      if (files->write_dN2Olyr) {
        fclose(files->fp_dN2Olyr);
      }

      if (files->write_cflows) {
        fclose(files->fp_cflows);
      }

      if (files->write_yrcflows) {
        fclose(files->fp_yrcflows);
      }

    }
