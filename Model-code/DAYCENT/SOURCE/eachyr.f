
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine eachyr

      implicit none
      include 'comput.inc'
      include 'const.inc'
      include 'dovars.inc'
      include 'isovar.inc'
      include 'jday.inc'
      include 'ligvar.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfs.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'seq.inc'
      include 'site.inc'
      include 'wth.inc'
      include 'zztim.inc'

c ... Perform tasks that only need to be done once a year.

c ... Local variables
      integer  imo

c ... Correct for the rounding problem with time. The value for time
c ... drifts downward during long runs since dt=1/12 cannot be represented
c ... precisely.  At this point in the run, time should be a whole number.

c ... Changed increment value from .001 to .5 to correct error in time calculation
c ... occuring after year 8192. (mse 3/95).  New code if from Kendrick Killian.

      time = sign(int(abs(time)+.5),int(time))

c ... Reset annual accumulators to zero
      call annacc

c ... Wet-dry fixation of N

c ... Determine annual precipitation and annual PET
c ... For RAMS/Daily Century, prcann = average annual precip. 
c ... petann (output var) is computed in dailymoist. -mdh 12/96

      do 10 imo = 1, MONTHS
        prcann = prcann + precip(imo)
        agdefacm(imo) = -1.
        bgdefacm(imo) = -1.
10    continue

c ... N fixation in atmosphere
      wdfxa = epnfa(INTCPT)+epnfa(SLOPE)*MIN(prcann,80.0)
      if (wdfxa .lt. 0.) then
        wdfxa = 0.0
      endif

c ... Non-symbiotic N fixation in the soil
c ... Use annual ET unless it is the first timestep
c ... No longer using the intercept in the calculation.
c ... wdfxs = epnfs(INTCPT)+epnfs(SLOPE)*MIN(prcann,100.0)

      if (annet .eq. 0.0) then
        wdfxs = epnfs(SLOPE)*MIN(prcann,100.0)
      else 
        wdfxs = epnfs(2) * (annet - epnfs(1))
      endif
c ... Reset annual accumulator for evapotranspiration
      annet = 0
      if (wdfxs .lt. 0.)  then
        wdfxs = 0.0
      endif

      wdfx = wdfxa+wdfxs

c ... Atmospheric S deposition
      satmt = max(0.0, satmos(1) + satmos(2)*prcann)

      cisofr = 0.0
      cisotf = 0.0

c ... Initialize co2 effects
      call co2eff(time)

      if (cursys .ne. FORSYS) then
c ..... Determine what fraction of plant residue added this year will be lignin.
        call cmplig(cursys,fligni,wdlig,pltlig)
      endif

c ... Determine the number of days in each month.  The leap year exception
c ... will be handled in getwth.  

      do 110 imo = 1, 12
        dysimo(imo) = idysimo(imo)
        lstdy(imo) = ilstdy(imo)
        frstdy(imo) = ifrstdy(imo)
110   continue

      return
      end
