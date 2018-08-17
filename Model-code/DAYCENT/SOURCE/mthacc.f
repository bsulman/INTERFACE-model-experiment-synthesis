
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine mthacc(agdefacsum, bgdefacsum)

      implicit none
      include 'const.inc'
      include 'monprd.inc'
      include 'param.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'

c ... Argument declarations
      double precision agdefacsum, bgdefacsum

c ... Reset monthly accumulators.

c ... Local variables
      integer iel, ii, ilyr, iso

c ... Monthly mineralization accumulators for each element.
      do 35 iel = 1, MAXIEL
        do 15 ilyr = 1, 2
          strmnr(ilyr,iel) = 0.0
          metmnr(ilyr,iel) = 0.0
          s1mnr(ilyr,iel) = 0.0
          s2mnr(ilyr,iel) = 0.0
15      continue
        s3mnr(iel) = 0.0
        gromin(iel) = 0.0
        w1mnr(iel) = 0.0
        w2mnr(iel) = 0.0
        w3mnr(iel) = 0.0
35    continue

      do 20 ii = 1, 8
        stream(ii) = 0
20    continue
      pet = 0
      evap = 0
      tran = 0
      pttr = 0
      rain = 0
      agdefacsum = 0.0
      bgdefacsum = 0.0
      irract = 0.0
      runoff = 0.0
      do 102 ilyr = 1,nlayer
        amov(ilyr) = 0
102   continue

c ... Monthly co2 accumlators (10/92)
      do 25 iso = 1, 2
        mt1c2(iso) = 0.0
        mt2c2(iso) = 0.0
        st1c2(iso) = 0.0
        st2c2(iso) = 0.0
        st1uvc2(iso) = 0.0
        s11c2(iso) = 0.0
        s12c2(iso) = 0.0
        s21c2(iso) = 0.0
        s22c2(iso) = 0.0
        s3c2(iso)  = 0.0
        stduvc2(iso) = 0.0
        wd1c2(iso) = 0.0
        wd2c2(iso) = 0.0
        wd3c2(iso) = 0.0
25    continue

c ... Monthly accumulator for volatilization of N during
c ... harvest, senescence, and return from grazing animal waste,

      volpl = 0.0

c ... Monthly accumulator for symbiotic N fixation to track
c ... fixation for both grasses and trees 
      nfix = 0.0

c ... Monthly C production
      cprodc = 0.0
      cprodf = 0.0
      do 30 iel = 1, MAXIEL
        eprodc(iel) = 0.0
        eprodf(iel) = 0.0
30    continue

c ... Monthly accumulator for soil surface temperature,
      stempmth = 0.0

c ... Monthly accumulators for maintenance respiration
      mrspflux(CRPSYS) = 0.0
      mrspflux(FORSYS) = 0.0
      mrspmth(CRPSYS) = 0.0
      mrspmth(FORSYS) = 0.0
      do 40 ii = 1, CPARTS
        cmrspflux(ii) = 0.0
40    continue
      do 45 ii = 1, FPARTS
        fmrspflux(ii) = 0.0
45    continue
      sumrsp = 0.0

c ... Monthly accumulators for growth respiration
      grspflux(CRPSYS) = 0.0
      grspflux(FORSYS) = 0.0
      grspmth(CRPSYS) = 0.0
      grspmth(FORSYS) = 0.0
      do 50 ii = 1, CPARTS
        cgrspflux(ii) = 0.0
50    continue
      do 55 ii = 1, FPARTS
        fgrspflux(ii) = 0.0
55    continue

c ... Monthly respiration from decomposition output variables
      respmth(UNLABL) = 0.0
      respmth(LABELD) = 0.0

c ... Monthly soil respiration output variables
      srspmth(CRPSYS) = 0.0
      srspmth(FORSYS) = 0.0

c ... Monthly autotrophic respiration output variables
      arspmth(CRPSYS,UNLABL) = 0.0
      arspmth(CRPSYS,LABELD) = 0.0
      arspmth(FORSYS,UNLABL) = 0.0
      arspmth(FORSYS,LABELD) = 0.0

c ... Monthly trace gas accumulator output variables,
      N2O_month = 0.0
      NO_month = 0.0
      N2_month = 0.0
      ch4_month = 0.0
      nit_amt_month = 0.0
      tave = 0.0

c ... Monthly wet-dry fixation of N
      wdfxma = 0.0
      wdfxms = 0.0

      return
      end
