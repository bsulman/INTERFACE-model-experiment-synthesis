
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


c                           DISCLAIMER
c
c        Neither the Great Plains System Research Unit - USDA (GPSR) nor
c     Colorado State University (CSU) nor any of their employees make
c     any warranty or assumes any legal liability or responsibility for
c     the accuracy, completeness, or usefulness of any information,
c     apparatus, product, or process disclosed, or represents that its
c     use would not infringe privately owned rights.  Reference to any
c     special commercial products, process, or service by tradename,
c     trademark, manufacturer, or otherwise, does not necessarily
c     constitute or imply endorsement, recommendation, or favoring by  
c     the GPSR or CSU.  The views and opinions of the authors do not
c     necessarily state or reflect those of GPSR or CSU and shall not 
c     be used for advertising or product endorsement. 

      program main

c ... Daily Century Stand Alone Soil Organic Matter / Trace Gas Model
c ... Simulation of carbon, nitrogen, phosphorous, and sulfur cycling 
c.....excluding plant growth submodel, soil water submodel, and soil temperature submodel
c ... Project - Community Land Model / DayCent linkage
c ... Scientists - Gordon Bonan and Bill Parton
c ... Programmer - Melannie Hartman

c ... State variables and flows are grams/m2.

      implicit none
      include 'cflows.inc'
      include 'const.inc'
      include 'dovars.inc'
      include 'jday.inc'
      include 'monprd.inc'
      include 'param.inc'
      include 'plot1.inc'
      include 'plot3.inc'
      include 't0par.inc'
      include 'timvar.inc'
      include 'wth.inc'
      include 'zztim.inc'

c ...               (unit 1) = plot/print file used by modaid (unformatted)
c ... <site>.100    (unit 7) = parameter values and initial values for 
c ...                          state variables; see subroutine sitein.
c ... fix.100       (unit 8) = fixed parameter values values for 
c ...                          state variables; see subroutine fixin.
c ...               (unit 9) = a file of weather data read in subroutines 
c ...                          wthini, weathr
c ... cdi file     (unit 18) = file with climate decomposition indices (CDI) for 
c ...                          each decomp timestep
c ... nflux.csv    (unit 70) = N2/N2O fluxes computed by Trace Gas Model
c ... daily.csv    (unit 80) = pet, defac, stemp, and snowpack water content
c ...                          computed by Trace Gas Model
c ... summary.csv  (unit 90) = tmax, tmin, prec, N2O flux, NO flux, CH4, and
c ...                          gross nitrification computed by Trace Gas Model

c ... If you're getting floating point errors mentioned after you exit
c ... Century, uncomment the following lines, recompile, run Century
c ... in dbx with the 'catch FPE' option to find the offending code.
c ... You can also run Century outside of dbx, in which case you will
c ... get messages on your screen giving you an fpe error code (see
c ... the Floating Point Programmer's Guide, p20) and a not-very-
c ... useful-to-3rd-or-4th-generation-language-programmers location. 
c ... The error handler 'mysigal' is a Fortran callable C routine 
c ... written by Martin Fowler; it can be replaced by any user written
c ... handler or any of several library handlers, the most useful 
c ... probably being SIGFPE_ABORT.  The calls to ieee_handler won't 
c ... compile using poa's binaries.

c      external mysigal
c      ieeer=ieee_handler('set','invalid',SIGFPE_ABORT)
c      ieeer=ieee_handler('set','division',mysigal)
c      ieeer=ieee_handler('set','overflow',mysigal)
c      ieeer=ieee_handler('set','underflow',SIGFPE_ABORT)

c ... You probably won't want to uncomment the following line; inexact
c ... floating point signals occur all over the place.

c      ieeer=ieee_handler('set','inexact',mysigal)

c ... Local variables
      double precision month_vals(12), neg_month_vals(12)

      data month_vals /0.08, 0.17, 0.25, 0.33, 0.42, 0.50, 0.58, 0.67,
     &                 0.75, 0.83, 0.92, 1.0/
      data neg_month_vals /-0.92, -0.83, -0.75, -0.67, -0.58, -0.50,
     &                     -0.42, -0.34, -0.25, -0.17, -0.08, 0.0/

c ... Read command line arguments, initialize the model
      call detiv

c ... Write out starting values
      call wrtbin(time)

c ... Update month
20    continue
      month = mod(month,12) + 1

c ... If time is greater than the ending time for the current block,
c ... read the next block
      if ((abs(time - blktnd) .lt. (0.5 * dt)) .and.
     &    (abs(time - tend)   .gt. (0.5 * dt))) then
        call readblk()
      endif

c ... Perform annual tasks
      if (month .eq. 1) then
        call eachyr
      endif

c ... The main driver for the model; call decomposition, growth, etc.
      call simsom()

c ... Write yearly output
c ... Add output for the N2 flux for the year and convert fluxes to
c ... g/m^2, cak - 01/16/03
      if ((time .ge. strplt) .and. (month .eq. 12)) then
        call wrtyearsum(time, N2O_year/10000, NO_year/10000,
     &                  N2_year/10000, CH4_year/10000,
     &                  nit_amt_year/10000, annppt)
        call wrtyrcflows(time, asom11tosom21, asom12tosom22,
     &                   asom12tosom3, asom21tosom11, asom21tosom22,
     &                   asom22tosom12, asom22tosom3, asom3tosom12,
     &                   ametc1tosom11, ametc2tosom12, astruc1tosom11,
     &                   astruc1tosom21, astruc2tosom12,
     &                   astruc2tosom22, awood1tosom11, awood1tosom21,
     &                   awood2tosom11, awood2tosom21, awood3tosom12,
     &                   awood3tosom22)
      endif

c ... Write monthly trace gas output, cak - 05/14/42
      if (time .ge. strplt) then
        call wrttgmonth(time, N2O_month/10000, NO_month/10000,
     &                  N2_month/10000, CH4_month/10000,
     &                  nit_amt_month/10000)
      endif

c ... Update time
c ... Add calculation to handle time drift caused by inexact floating
c ... point addition, cak - 08/23/01
c      time = time + dt
c ... Add calculation to handle incrementing the month for negative years,
c ... cak - 03/29/02
      if (time .ge. 0.0) then
        time = int(time) + month_vals(month)
      else
        time = int(time) + neg_month_vals(month)
        if (month .eq. 1) then
          time = time + 1.0
        endif
      endif
      if (time .ge. -1.0e-07 .and. time .le. 1.0e-07) then
        time = 0.0
      endif
c     time2 = time

c ... Write out values
      if ((tplt - time) .lt. (dt * 0.5)) then
        call wrtbin(time)
        tplt = time + dtpl
      endif

c ... Run for tend years
      if ((tend-time) .gt. (dt*.5)) then
        goto 20
      endif

c ... Write out final values
      call wrtbin(time)

c ... Close data files

c ... Close the weather file
      close(unit=9)
c ... Close the schedule file
      close(unit=15)
c ... Close the cdi file (unit=18)
      close(unit=18)
c ... Close the summary.csv file (unit=90)
      close(unit=90)
c ... Close the cbalance file (unit=123)
      close(unit=123)
c ... Close the nbalance file (unit=124)
      close(unit=124)
c ... Close the decomp.csv (unit=125)
      close(unit=125)
c ... Close the diagnostic file (unit=128)
      close(unit=128)
c ... Close *.out and *.csv files opened by the initsw subroutine
      call closefiles()

c ... Mark end of file
      endfile(unit=1)

c ... Close binary file
      close(unit=1)

      write(*,*) 'Execution success.'
      STOP 'Execution success.'

      end
