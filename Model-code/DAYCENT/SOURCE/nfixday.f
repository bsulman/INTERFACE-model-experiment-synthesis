
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine nfixday(ppt, atmosNdep, nonSymSoilNfix)

c ... Compute the amount of atmospheric and non-symbiotic soil N-fixation 
c ... for the current day and add it to minerl(1,N)

      implicit none
      include 'const.inc'
      include 'npool.inc'
      include 'param.inc'
      include 'parfx.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'timvar.inc'
      include 'zztim.inc'
     
c ... Formal parameters
      double precision ppt, atmosNdep, nonSymSoilNfix

c ... Function declarations
      double precision fsfunc
      external  fsfunc

c ... FORMAL PARAMETERS
c ... Input
c ...   ppt		- precipitation for the current day (cm)
c ... Output
c ...   atmosNdep	- atmospheric N deposition (gN/m^2)
c ...   nonSymSoilNfix	- non-symbiotic soil N fixation (gN/m^2)

c ... N FIXATION VARIABLES (plot1.inc)
c ... "wdfx" stands for wet/dry fixation of N
c ...   wdfxa 	- average annual rate of atmospheric N deposition (gN/m^2)
c ...   wdfxs 	- average annual rate of non-symbiotic soil N fixation (gN/m^2)
c ...   wdfx = wdfxa + wdfxs
c ...   wdfxma 	- actual amount of atmospheric N deposition this month (gN/m^2)
c ...   wdfxms 	- actual amount of non-symbiotic soil N fixation this month (gN/m^2)
c ...   wdfxm = wdfxma + wdfxms
c ...   wdfxaa 	- accumulator for the actual amount of atmospheric N deposition this year (gN/m^2)
c ...   wdfxas 	- accumulator for the actual amount of non-symbiotic soil N fixation this year (gN/m^2)


c ... Local Variables
      integer clyr
      double precision frac_nh4, frac_no3
      double precision biof, fwdfx, fxbiom
      double precision wdfnp, wdfxm
      double precision tbiom
      double precision totNfix
      character subname*10


c ... N Fixation
c ... Does not take into account the effect of irrigation

      if (nsnfix .eq. 1 .and. nelem .ge. P) then

c.......Non symbiotic N fixation should be based on N:P ratio in mineral pool
c ..... Compute mineral N:P ratio for N-Fixation (use suface layer only)
c ..... rnpml1 is used to control soil N-fixation using a regression
c ..... equation based on Kansas data. This ratio is used only if nelem = 2.
c ..... rnpml1 is flagged if either minerl(1,1) or minerl(1,2) is zero.

        rnpml1 = minerl(1,N)/minerl(1,P)*
     &           fsfunc(minerl(1,P),pslsrb,sorpmx)

c ..... Wet-dry fixation of nitrogen -- monthly flow
c ..... Atmospheric fixation is split between monthly dry fall and
c ..... wdfnp is the N:P ratio control function for non-symbiotic
c ..... soil fixation.
c ..... Both wdfnp and fxbiom are relative functions
c ..... which vary from 0 to 1.
c ..... wdfnp computed as a negative natural log of rnpml1
c ..... symnfx is the symbiotic N fixation by legumes derived from Cole and
c ..... Heil (1981) using data from Walker et al. 1959.
        if (rnpml1 .eq. 0) then
          wdfnp = 1.
        else
          wdfnp = min(1., ((-alog(real(rnpml1)))/fxnpb)-.2)
        endif

c ..... The definitions of fxmca and fxmcb originally refered to water,
c ..... not biomass. (Alister 11/91)
        tbiom = aglivc + stdedc + strucc(SRFC)
        biof  = fxmca + tbiom * fxmcb
        fxbiom = 1 - biof
        fxbiom = min(1.,fxbiom)
        if (wdfnp .lt. 0 .or. fxbiom .lt. 0 .or. stemp .lt. 7.5) then
          fwdfx = 0.0
        else
          fwdfx = wdfnp * fxbiom
        endif

c ..... Compute N-fixation for the month

c ..... Wet fall depending on the monthly precipitation (wdfxma)
c ..... 30.42 is the average number of days in a month

        nonSymSoilNfix = fxmxs * fwdfx / 30.42
        atmosNdep = wdfxa * ppt / prcann
        totNfix = nonSymSoilNfix + atmosNdep
        wdfxms = wdfxms + nonSymSoilNfix
        wdfxma = wdfxma + atmosNdep
        wdfxm = wdfxma + wdfxms
        wdfxas = wdfxas + nonSymSoilNfix
        wdfxaa = wdfxaa + atmosNdep

c       wdfxms = fxmxs * fwdfx
c       wdfxma = wdfxa * precip(month) / prcann
c       wdfxas = wdfxas + wdfxms
c       wdfxaa = wdfxaa + wdfxma

c ..... Compute annual N-fixation accumulators for atmosphere and soils
        frac_nh4 = 0.5
        frac_no3 = 0.5

        clyr = 1
        subname = 'nfixation '
        call update_npool(clyr, totNfix, frac_nh4, frac_no3, ammonium,
     &                    nitrate, subname)
        call flow(esrsnk(N),minerl(1,N),time,totNfix)

c.....Non symbiotic N fixation is based on annual precipitation.
      else
        frac_nh4 = 0.5
        frac_no3 = 0.5
        nonSymSoilNfix = wdfxs * ppt / prcann
        atmosNdep = wdfxa * ppt / prcann
        totNfix = nonSymSoilNfix + atmosNdep
        wdfxms = wdfxms + nonSymSoilNfix
        wdfxma = wdfxma + atmosNdep
        wdfxm = wdfxma + wdfxms
        wdfxas = wdfxas + nonSymSoilNfix
        wdfxaa = wdfxaa + atmosNdep

c       wdfxms = wdfxs * precip(month) / prcann
c       wdfxma = wdfxa * precip(month) / prcann
c       wdfxas = wdfxas + wdfxms
c       wdfxaa = wdfxaa + wdfxma

        clyr = 1
        call update_npool(clyr, totNfix, frac_nh4, frac_no3, ammonium,
     &                    nitrate, subname)
        call flow(esrsnk(N),minerl(1,N),time,totNfix)
      endif

      return
      end
