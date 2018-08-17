
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine simsom()

      implicit none
      include 'cflows.inc'
      include 'comput.inc'
      include 'const.inc'
      include 'dovars.inc'
      include 'fertil.inc'
      include 'isovar.inc'
      include 'jday.inc'
      include 'ligvar.inc'
      include 'monprd.inc'
      include 'npool.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfs.inc'
      include 'parfx.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'seq.inc'
      include 'site.inc'
      include 'timvar.inc'
      include 'wth.inc'
      include 'wthdaily.inc'
      include 'zztim.inc'

c ... Simulate flow of carbon, nitrogen, phosphorous, and sulfur.
c ... It calculates N fixation, and runs the decompostion and tracegas models
c ... in the daily timestep loop.

c ... Function declarations
      double precision fsfunc, ramp
      external  fsfunc, ramp

c ... Local variables

      integer iel, lyr, ii
      integer curday
      integer isdecid, isagri
      double precision cmn
      double precision fsol
      double precision agdefacsum, bgdefacsum, tfrac
      double precision sradKJ, soilsrad
      double precision tfunc, agwfunc, bgwfunc
      double precision prevstream(8)
      double precision respsum(2)
      double precision petdly
      double precision oiresp1(2), oiresp2(2), newoiresp(2)
      double precision oeresp1(2), oeresp2(2), newoeresp(2)
      double precision slitrsp1(2), slitrsp2(2), newslitrsp(2)
      double precision mnrlrsp1(2), mnrlrsp2(2), newmnrlrsp(2)
      double precision hetresp1(2), hetresp2(2), newhetresp(2)
      double precision doiresp, doeresp, dslitrsp, dmnrlrsp, dhresp
      double precision grass_lai, tree_lai
      double precision CONVLAI
      double precision rpeff

c ... CONVLAI = biomass needed to produce an LAI of 1 (g/m**2) 

c ... Saved variables
      data isagri /0/
      save isagri
      save respsum

c ... Set aminrl for use in routines called from decomp
c ... move this to daily loop! -MDH 1/9/2012
c     do 10 iel = 1, nelem
c       if (iel .eq. P) then
c         fsol = fsfunc(minerl(1,P), pslsrb, sorpmx)
c       else
c         fsol = 1.0
c       endif
c       aminrl(iel) = minerl(1,iel) * fsol
c10    continue

c ... N Fixation
c ... This call was replaced by nfixday() in dailymoist
c     call nfixmonth()

c ... BEGIN MONTHLY INITIALIZATION

c ... Reset the monthly accumulators
      call mthacc(agdefacsum, bgdefacsum)

      do 20 ii = 1, 8
        prevstream(ii) = 0
20    continue

c ... BEGIN DAILY LOOP...

      do 200 curday = frstdy(month), lstdy(month)

c ..... Calculate the root priming effect on the som2c(2) decomposition
c ..... rate, CAK - 01/28/2014, MDH - 05/02/2016
c ..... Only options 0 and 2 are relevant to the SOM model since root
c ..... respiration is not simulated. dshresp = daily heterotrophic soil 
c ..... respiration from the previous day (calculated in dailymoist).
        rpeff = 1.0
c ..... Crop/grass root priming effect
        if (crpindx .eq. 0) then
c ....... No root priming effect on som2c(2) decomposition
          rpeff = 1.0

c       elseif (crpindx .eq. 1) then
c ....... Root priming effect on som2c(2) decomposition based on total
c ....... soil respiration
c         rpeff = ramp(dsresp, crpcmn, crpmnmul, crpcmx, crpmxmul)

        elseif (crpindx .eq. 2) then
c ....... Root priming effect on som2c(2) decompostion based on soil
c ....... heterotrophic respiration only
          rpeff = ramp(dshresp, crpcmn, crpmnmul, crpcmx, crpmxmul)

c       elseif (crpindx .eq. 3) then
c ....... Root priming effect on som2c(2) decomposition based on fine
c ....... root production 
c         rpeff = ramp(mcprd(BELOWJ)+mcprd(BELOWM), crpcmn, crpmnmul,
c    &                  crpcmx, crpmxmul)
        endif

c ..... Call schedl to determine scheduling options for this day of the year
        call schedl(curday)

        tfrac = 1.0/dble(dysimo(month))

        do 10 iel = 1, nelem
          if (iel .eq. P) then
            fsol = fsfunc(minerl(1,P), pslsrb, sorpmx)
          else
            fsol = 1.0
          endif
          aminrl(iel) = minerl(1,iel) * fsol
10      continue

c ..... Get a day's worth of weather from the weather file
        call getwth(curday, month,  tmn2m, tmx2m, 
     &              tempmax, tempmin, avgtemp, ppt,
     &              solrad, rhumid, windsp, srad, vpd)

        call calcpet(curday, month, tempmin(curday), tempmax(curday), 
     &               avgtemp(curday), solrad(curday), rhumid(curday), 
     &               windsp(curday), snow, usexdrvrs, fwloss, 
     &               sitlat, tmn2m, tmx2m, petdly)

        rain = rain + ppt(curday)
        annppt = annppt + ppt(curday)

c ..... Effect of cultivation on decomposition (used in decomp routine)
c ..... cltfac is this month's effect of cultivation on decomposition
c ..... of som1, som2, som3, and structural.  It is set to clteff
c ..... in months when cultivation occurs; otherwise it equals 1.
c ..... clteff is the effect of cultivation on decomposition read from
c ..... the cult.100 file
        if (docult .and. (cultday .eq. curday) .or.
     &      (cultcnt .gt. 0. .and. cultcnt .lt. 31)) then
          do 34 ii = 1, 4
            cltfac(ii) = clteff(ii)
34        continue
          if (cultday .eq. curday) then
            cultcnt = 0
          endif
          cultcnt = cultcnt + 1
        else
          do 35 ii = 1, 4
            cltfac(ii) = 1.0
35        continue
          cultcnt = 0
        endif

c ..... Reset the reduction factor on nitrification rates due to
c ..... nitrification inhibitors as necessary, cak - 12/01/03
        if ((fertcnt .gt. 0) .and. (fertcnt .lt. (ninhtm * 7.0))) then
          fertcnt = fertcnt + 1
        else
          fertcnt = 0
          nreduce = 1.0
        endif

c ..... Save previous stream values from inorganic and organic leaching.  -mdh 1/97
        do 85 iel = 1, nelem
          prevstream(iel+1) = stream(iel+1)
          prevstream(iel+5) = stream(iel+5)
85      continue
        strm5l = 0.0
        strm5u = 0.0

c ..... Track respiration over the daily timestep
        hetresp1(UNLABL) = mt1c2(UNLABL) + mt2c2(UNLABL) + 
     &                     st1c2(UNLABL) + st2c2(UNLABL) +
     &                     s11c2(UNLABL) + s12c2(UNLABL) +
     &                     s21c2(UNLABL) + s22c2(UNLABL) +
     &                     s3c2(UNLABL)  + wd1c2(UNLABL) +
     &                     wd2c2(UNLABL) + wd3c2(UNLABL)
        hetresp1(LABELD) = mt1c2(LABELD) + mt2c2(LABELD) +
     &                     st1c2(LABELD) + st2c2(LABELD) +
     &                     s11c2(LABELD) + s12c2(LABELD) +
     &                     s21c2(LABELD) + s22c2(LABELD) +
     &                     s3c2(LABELD)  + wd1c2(LABELD) +
     &                     wd2c2(LABELD) + wd3c2(LABELD)
        mnrlrsp1(UNLABL) = mt2c2(UNLABL) + st2c2(UNLABL) +
     &                     s12c2(UNLABL) + s22c2(UNLABL) +
     &                     s3c2(UNLABL)  + wd3c2(UNLABL)
        mnrlrsp1(LABELD) = mt2c2(LABELD) + st2c2(LABELD) +
     &                     s12c2(LABELD) + s22c2(LABELD) +
     &                     s3c2(LABELD)  + wd3c2(LABELD)
        oiresp1(UNLABL) = mt1c2(UNLABL) + st1c2(UNLABL) +
     &                    s11c2(UNLABL)
        oiresp1(LABELD) = mt1c2(LABELD) + st1c2(LABELD) +
     &                    s11c2(LABELD)
        oeresp1(UNLABL) = s21c2(UNLABL)
        oeresp1(LABELD) = s21c2(LABELD)
        slitrsp1(UNLABL) = mt1c2(UNLABL) + st1c2(UNLABL) +
     &                     s11c2(UNLABL) + s21c2(UNLABL) +
     &                     st1uvc2(UNLABL) + stduvc2(UNLABL)
        slitrsp1(LABELD) = mt1c2(LABELD) + st1c2(LABELD) +
     &                     s11c2(LABELD) + s21c2(LABELD) +
     &                     st1uvc2(LABELD) + stduvc2(LABELD)


c ... Solar radiation at the bottom of the plant canopy
c ... Calculate biomass values
      if (cursys .eq. FORSYS) then
c ..... Live biomass
        aglivb = rleavc * 2.5
c ..... Standing dead biomass
        stdead = 0.0
      elseif (cursys .eq. SAVSYS) then
c ..... Live biomass
        aglivb = (rleavc + aglivc) * 2.5
c ..... Standing dead biomass
        stdead = stdedc * 2.5
      else
c ..... Live biomass
        aglivb = aglivc * 2.5
c ..... Standing dead biomass
        stdead = stdedc * 2.5
      endif
      stcrlai = (aglivb + stdead) / 80.0

c ... Solar radiation in KJ
      sradKJ = srad(curday)*W2KJ
      soilsrad = sradKJ * exp(-0.5 * stcrlai)

      CONVLAI = 80.0
      grass_lai = (aglivc * 2.5) / CONVLAI;
      tree_lai = (rleavc * 2.5) * btolai;
      if (grass_lai .lt. 0.0) then
          grass_lai = 0.0
      endif
      if (tree_lai .lt. 0.0) then
          tree_lai = 0.0
      endif

c ... Are we running a decidious forest?
      if ((cursys .eq. FORSYS) .and. (decid .eq. 1)) then
        isdecid = 1
      else
        isdecid = 0
      endif
c ... Once cultivation, fertilization, or harvesting occurs in a system the 
c ... agricultural effect needs to be "turned on".  Once this effect on
c ... methane oxidation has been invoked it stays in place.
      if (isagri .eq. 0) then
        if (docult) isagri = 1
      endif

      call dailymoist(curday, sradKJ, soilsrad, 
     &                grass_lai, tree_lai, isdecid, isagri,
     &                agwfunc, tfunc, bgwfunc, petdly, rpeff)

c ..... Available nutrients
c ..... tminrl is the total amount of each element available in mineral form.

        do 80 iel = 1, nelem
          tminrl(iel) = 0.
          if (iel .eq. P) then
            fsol = fsfunc(minerl(SRFC,P), pslsrb, sorpmx)
          else 
            fsol = 1.0
          endif

          do 70 lyr = 1, nlayer
c ......... Plants can only uptake from a layer with a positive
c ......... value, so only the positive layers are summed here.
            if (minerl(lyr,iel) .gt. 0.)  then
              tminrl(iel) = tminrl(iel) + minerl(lyr,iel) * fsol
            endif
70        continue
80      continue

c ..... Compute the fraction of labile (non-sorbed) P in the surface
c ..... layer available to plants

        favail(2) = max(favail(4),
     &              min(favail(4) + minerl(SRFC,N)*
     &              (favail(5) - favail(4)) / favail(6),
     &              favail(5)))

c ... Accumulate daily stream flow, drainage from each layer, pet, 
c ... evaporation, and transpiration by month
      pet = pet + petdly
      petann = petann + petdly
c     do 104 ilyr = 1,nlayer
c       amov(ilyr) = amov(ilyr) + amovdly(ilyr)
c104   continue

c ... calculate defacm(month) in subroutine simsom. -mdh 10/94
      agdefacsum = agdefacsum + agdefac
      bgdefacsum = bgdefacsum + bgdefac

c ... SCALE: convert gN/m^2 to gN/ha
c     SCALE = 10000.0

c ... Annual Accumulators

c ... Accumulate yearly trace gas output
c     N2O_year = N2O_year + (Nn2oflux+Dn2oflux)*SCALE
c     NO_year = NO_year + NOflux*SCALE
c     N2_year = N2_year + Dn2flux*SCALE
c     CH4_year = CH4_year + CH4
c     nit_amt_year = nit_amt_year + nit_amt*SCALE

c ... Accumulate output variables that are tracking the carbon flows
c ... due to decomposition
      ametc1tosom11  = ametc1tosom11  + metc1tosom11
      ametc2tosom12  = ametc2tosom12  + metc2tosom12
      astruc1tosom11 = astruc1tosom11 + struc1tosom11
      astruc1tosom21 = astruc1tosom21 + struc1tosom21
      astruc2tosom12 = astruc2tosom12 + struc2tosom12
      astruc2tosom22 = astruc2tosom22 + struc2tosom22
      asom11tosom21  = asom11tosom21  + som11tosom21
      asom12tosom22  = asom12tosom22  + som12tosom22
      asom12tosom3   = asom12tosom3   + som12tosom3
      asom21tosom11  = asom21tosom11  + som21tosom11
      asom21tosom22  = asom21tosom22  + som21tosom22
      asom22tosom12  = asom22tosom12  + som22tosom12
      asom22tosom3   = asom22tosom3   + som22tosom3
      asom3tosom12   = asom3tosom12   + som3tosom12
      awood1tosom11  = awood1tosom11  + wood1tosom11
      awood1tosom21  = awood1tosom21  + wood1tosom21
      awood2tosom11  = awood2tosom11  + wood2tosom11
      awood2tosom21  = awood2tosom21  + wood2tosom21
      awood3tosom12  = awood3tosom12  + wood3tosom12
      awood3tosom22  = awood3tosom22  + wood3tosom22

c ... Accumulate monthly trace gas output, cak - 05/14/04
c     N2O_month = N2O_month + (Nn2oflux+Dn2oflux)*SCALE
c     NO_month = NO_month + NOflux*SCALE
c     N2_month = N2_month + Dn2flux*SCALE
c     CH4_month = CH4_month + CH4
c     nit_amt_month = nit_amt_month + nit_amt*SCALE
c ... Growing season accumulator for N2O flux, cak - 06/06/2008
c     n2oacc = n2oacc + (Nn2oflux+Dn2oflux)
c     n2omth(month) = n2omth(month) + (Nn2oflux+Dn2oflux)

c ..... Track respiration over the daily timestep
        hetresp2(UNLABL) = mt1c2(UNLABL) + mt2c2(UNLABL) + 
     &                     st1c2(UNLABL) + st2c2(UNLABL) +
     &                     s11c2(UNLABL) + s12c2(UNLABL) +
     &                     s21c2(UNLABL) + s22c2(UNLABL) +
     &                     s3c2(UNLABL)  + wd1c2(UNLABL) +
     &                     wd2c2(UNLABL) + wd3c2(UNLABL)
        hetresp2(LABELD) = mt1c2(LABELD) + mt2c2(LABELD) +
     &                     st1c2(LABELD) + st2c2(LABELD) +
     &                     s11c2(LABELD) + s12c2(LABELD) +
     &                     s21c2(LABELD) + s22c2(LABELD) +
     &                     s3c2(LABELD)  + wd1c2(LABELD) +
     &                     wd2c2(LABELD) + wd3c2(LABELD)
        newhetresp(UNLABL) = hetresp2(UNLABL) - hetresp1(UNLABL)
        newhetresp(LABELD) = hetresp2(LABELD) - hetresp1(LABELD)
        dhresp = (hetresp2(UNLABL) + hetresp2(LABELD)) -
     &           (hetresp1(UNLABL) + hetresp1(LABELD))

        mnrlrsp2(UNLABL) = mt2c2(UNLABL) + st2c2(UNLABL) +
     &                     s12c2(UNLABL) + s22c2(UNLABL) +
     &                     s3c2(UNLABL)  + wd3c2(UNLABL)
        mnrlrsp2(LABELD) = mt2c2(LABELD) + st2c2(LABELD) +
     &                     s12c2(LABELD) + s22c2(LABELD) +
     &                     s3c2(LABELD)  + wd3c2(LABELD)
        newmnrlrsp(UNLABL) = mnrlrsp2(UNLABL) - mnrlrsp1(UNLABL)
        newmnrlrsp(LABELD) = mnrlrsp2(LABELD) - mnrlrsp1(LABELD)
        dmnrlrsp = (mnrlrsp2(UNLABL) + mnrlrsp2(LABELD)) -
     &             (mnrlrsp1(UNLABL) + mnrlrsp1(LABELD))

        oiresp2(UNLABL) = mt1c2(UNLABL) + st1c2(UNLABL) +
     &                    s11c2(UNLABL)
        oiresp2(LABELD) = mt1c2(LABELD) + st1c2(LABELD) +
     &                    s11c2(LABELD)
        newoiresp(UNLABL) = oiresp2(UNLABL) - oiresp1(UNLABL)
        newoiresp(LABELD) = oiresp2(LABELD) - oiresp1(LABELD)
        doiresp = (oiresp2(UNLABL) + oiresp2(LABELD)) -
     &            (oiresp1(UNLABL) + oiresp1(LABELD))

        oeresp2(UNLABL) = s21c2(UNLABL)
        oeresp2(LABELD) = s21c2(LABELD)
        newoeresp(UNLABL) = oeresp2(UNLABL) - oeresp1(UNLABL)
        newoeresp(LABELD) = oeresp2(LABELD) - oeresp1(LABELD)
        doeresp = (oeresp2(UNLABL) + oeresp2(LABELD)) -
     &            (oeresp1(UNLABL) + oeresp1(LABELD))

        slitrsp2(UNLABL) = mt1c2(UNLABL) + st1c2(UNLABL) +
     &                     s11c2(UNLABL) + s21c2(UNLABL) +
     &                     st1uvc2(UNLABL) + stduvc2(UNLABL)
        slitrsp2(LABELD) = mt1c2(LABELD) + st1c2(LABELD) +
     &                     s11c2(LABELD) + s21c2(LABELD) +
     &                     st1uvc2(LABELD) + stduvc2(LABELD)
        newslitrsp(UNLABL) = slitrsp2(UNLABL) - slitrsp1(UNLABL)
        newslitrsp(LABELD) = slitrsp2(LABELD) - slitrsp1(LABELD)
        dslitrsp = (slitrsp2(UNLABL) + slitrsp2(LABELD)) -
     &             (slitrsp1(UNLABL) + slitrsp1(LABELD))

c ..... Accumulate leached C,N,P,S
        csrsnk(UNLABL) = csrsnk(UNLABL) + strm5u
        csrsnk(LABELD) = csrsnk(LABELD) + strm5l
        stream(5) = stream(5) + strm5u + strm5l
        do 90 iel = 1, nelem
          esrsnk(iel) = esrsnk(iel)+(stream(iel+1)-prevstream(iel+1))+
     &                  (stream(iel+5)-prevstream(iel+5))
90      continue

      if (time .ge. strplt) then

c ...   Write to daily output files

c ...   Write to file daily.csv
        write(80,86)time,curday,petdly,tfunc,agwfunc,bgwfunc,
     &    agdefac,bgdefac,stemp,snow,snlq,srad(curday)
c
c ...   Write to file summary.csv
c       write(90,95) time,curday,tempmax(curday),tempmin(curday),
c    &    ppt(curday),(Nn2oflux+Dn2oflux)*SCALE,NOflux*SCALE,
c    &    CH4, nit_amt*SCALE, CO2resp*SCALE
c
c       call wrtsoiln(time, curday, ammonium, nitrate)
c       call wrtco2(time, curday, co2_conc)
c       call wrtdn2lyr(time, curday, dN2lyr)
c       call wrtdn2olyr(time, curday, dN2Olyr)
        call wrtcflows(time, curday, som11tosom21, som12tosom22,
     &                 som12tosom3, som21tosom11, som21tosom22,
     &                 som22tosom12, som22tosom3, som3tosom12,
     &                 metc1tosom11, metc2tosom12, struc1tosom11,
     &                 struc1tosom21, struc2tosom12, struc2tosom22,
     &                 wood1tosom11, wood1tosom21, wood2tosom11,
     &                 wood2tosom21, wood3tosom12, wood3tosom22)
        call wrtsoilc(time, curday, metabc(2), strucc(2), som1c(1), 
     &                  som1c(2), som2c(1), som2c(2), som3c)


86    format(f10.4,',',i4,',',7(f12.4,','),2(f7.4,','),f12.4)
c95    format(f10.4,',',1x,i4,',',1x,3(f8.2,',',1x),5(f12.4,',',1x))

        endif

200   continue
c ... END DAILY LOOP

      agdefacm(month) = agdefacsum/dble(dysimo(month))
      bgdefacm(month) = bgdefacsum/dble(dysimo(month))
      agdefac = agdefacm(month)
      bgdefac = bgdefacm(month)

c ... Calculate monthly respiration from decomposition for output
      if (month .eq. 1) then
        respmth(1) = resp(1)
        respmth(2) = resp(2)
        respsum(1) = respmth(1)
        respsum(2) = respmth(2)
      else
        respmth(1) = resp(1) - respsum(1)
        respmth(2) = resp(2) - respsum(2)
        respsum(1) = respsum(1) + respmth(1)
        respsum(2) = respsum(2) + respmth(2)
      endif  

c ... Annual co2 accumulators (10/92)
      ast1c2 = ast1c2 + st1c2(UNLABL) + st1c2(LABELD)
      ast2c2 = ast2c2 + st2c2(UNLABL) + st2c2(LABELD)
      amt1c2 = amt1c2 + mt1c2(UNLABL) + mt1c2(LABELD)
      amt2c2 = amt2c2 + mt2c2(UNLABL) + mt2c2(LABELD)
      as11c2 = as11c2 + s11c2(UNLABL) + s11c2(LABELD)
      as12c2 = as12c2 + s12c2(UNLABL) + s12c2(LABELD)
      as21c2 = as21c2 + s21c2(UNLABL) + s21c2(LABELD)
      as22c2 = as22c2 + s22c2(UNLABL) + s22c2(LABELD)
      as3c2  = as3c2  + s3c2(UNLABL)  + s3c2(LABELD)
      ast1uvc2 = ast1uvc2 + st1uvc2(UNLABL) + st1uvc2(LABELD)
      astduvc2 = astduvc2 + stduvc2(UNLABL) + stduvc2(LABELD)

c ... Annual Net Mineralization
      do 100 iel = 1, nelem

c ..... Net mineralization for the mineralizing compartments
c ..... The structural component of litter and the wood compartments
c ..... are not mineralizers.  They should not be added into cmn or sumnrs.
        cmn = metmnr(SRFC,iel) + metmnr(SOIL,iel) +
     &        s1mnr(SRFC,iel) + s1mnr(SOIL,iel) +
     &        s2mnr(SRFC,iel) + s2mnr(SOIL,iel) + s3mnr(iel)
        sumnrs(iel) = sumnrs(iel) + cmn

c ..... Annual soilnm is net mineralization in the soil.
        soilnm(iel) = soilnm(iel) + s1mnr(SOIL,iel) +
     &                s2mnr(SOIL,iel) + s3mnr(iel) +
     &                metmnr(SOIL,iel) + strmnr(SOIL,iel) + w3mnr(iel)

c ..... Annual total net mineralization
        tnetmn(iel) = tnetmn(iel) + cmn + 
     &                strmnr(SRFC,iel) + strmnr(SOIL,iel) +
     &                w1mnr(iel) + w2mnr(iel) + w3mnr(iel)
100   continue

c ... Annual Stream flow accumulators
      do 105 ii = 1, 8
        strmac(ii) = strmac(ii) + stream(ii)
105   continue

c ... Compute output variables for printing or plotting.
      call savarp

      return
      end
