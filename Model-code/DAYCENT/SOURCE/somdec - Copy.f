
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


C ... SOMDEC.F

      subroutine somdec(amovdly, dtm, t1, t2, newminrl, soilsrad)

      implicit none
      include 'cflows.inc'
      include 'comput.inc'
      include 'const.inc'
      include 'param.inc'
      include 'parfx.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'seq.inc'
      include 'zztim.inc'

c ... Argument declarations
      double precision amovdly(MAXLYR)
      double precision dtm, soilsrad
      double precision t1, t2
      double precision newminrl

c ... Soil Organic Matter Decomposition  
c ... Decompose SOM1 (surface and soil), SOM2, and SOM3.
c ... defac = decomposition factor based on water and
c ...         temperature computed in prelim and in cycle
c ... dtm   = time step (fraction of a year)

c ... Function declarations
      double precision agdrat, bgdrat
      logical   candec
      double precision catanf
      external  agdrat, bgdrat, candec, catanf

c ... Local variables
      integer   iel
      double precision accum(ISOS), cf3co2, cfs1s2, cfs1s3, 
     &          cfs2s1, cfs2s3, cfs3s1, cfsfs2, 
     &          cleach, co2los, linten, mnrflo, orgflow,
     &          pheff, rceof1, rceto1(MAXIEL), rceto2(MAXIEL),
     &          rceto3(MAXIEL), rcetob, tcflow
      double precision somc(2), mix
      double precision a, b, mti, x1, x2
      double precision rateDay, rate2

c ... *******************************************************************
c ... Initialize ACCUM even though it is not really used
      accum(UNLABL) = 0.0
      accum(LABELD) = 0.0
      somc(1) = 0.0
      somc(2) = 0.0

c ... Calculate microbial turnover increase (mti) based on incoming solar
c ... radiation, cak - 12/08/2009
      if (soilsrad .le. 0.0) then
        mti = 1.0
      elseif (soilsrad .gt. 30000) then
        mti = 5.0
      else
        a = 1.0
        b = 5.0
        x1 = 0.0
        x2 = 30000.0
        mti = (b - a) / (x2 - x1) * (soilsrad - x1) + a
      endif

c ... Surface SOM1 decomposes to SOM2 with CO2 loss
c ... Now goes to SURFACE SOM2 not SOIL SOM2, cak - 06/14/05
      if (som1c(SRFC) .gt. 1.e-07) then

c ..... Determine C/E ratios for flows to surface som2
        do 10 iel=1,nelem
          rceto2(iel) = agdrat(aminrl,varat21,iel)
10      continue

c ..... Compute pH effect on decomposition
        pheff = catanf(ph, 4.0, 0.5, 1.1, 0.7)
        pheff = min(pheff, 1.0)
        pheff = max(pheff, 0.0)

        rateDay = dec3(SRFC) * dtm
        call rateconv(rateDay, rate2, t1, t2)
c ..... Compute total C flow out of surface microbes.
        tcflow = som1c(SRFC) * agdefac * rate2 * pheff * mti
c ..... where
c .....   som1c(SRFC) =  unlabeled and labeled C in surface microbes
c .....                  (som1ci(SRFC,1)+som1ci(SRFC,2))
c .....   dec3(SRFC)  =  intrinsic decomposition rate of
c .....                  surface microbes

c ..... If decomposition can occur, schedule flows associated with respiration 
c ..... and decomposition 
        if (candec(nelem,aminrl,som1c(SRFC),som1e,2,SRFC,rceto2)) then

c ....... CO2 loss - Compute and schedule respiration flows.
          co2los = tcflow * p1co2(SRFC)
c ....... where 
c .......   p1co2(SRFC)  = set to p1co2a(SRFC) in prelim.
c .......   p1co2a(SRFC) = a fixed parameter;
c .......                  intercept parameter which controls flow from soil
c .......                  organic matter with fast turnover to CO2 (fraction
c .......                  of carbon lost to CO2 when there is no sand in the
c .......                  soil)
c ....... Changed csrsnk to s11c2 (10/92)

          call respir(co2los,2,SRFC,som1c,som1ci,s11c2,resp,
     &                som1e,minerl,gromin,s1mnr,newminrl)

c ....... Decompose Surface SOM1 to SOM2

c ....... cfsfs2 is C Flow from SurFace som1 to Som2
          cfsfs2 = tcflow - co2los
          som11tosom21 = som11tosom21 + cfsfs2

c ....... Partition and schedule C flows by isotope
          call csched(cfsfs2,som1ci(SRFC,LABELD),som1c(SRFC),
     &                som1ci(SRFC,UNLABL),som2ci(SRFC,UNLABL),
     &                som1ci(SRFC,LABELD),som2ci(SRFC,LABELD),
     &                1.0,accum)

c ....... Compute and schedule N, P, and S flows.

c ....... Update mineralization accumulators.
          do 20 iel=1,nelem
            call esched(cfsfs2,som1c(SRFC),rceto2(iel),
     &                  som1e(SRFC,iel),som2e(SRFC,iel),
     &                  minerl(SRFC,iel),mnrflo)
            call mnracc(mnrflo,gromin(iel),s1mnr(SRFC,iel))
            if (iel .eq. N) then
              newminrl = newminrl + mnrflo
            endif
20        continue
        endif
      endif

c ... End of SOM1 (surface layer) Decomposition

c ... *******************************************************************

c ... Soil SOM1 decomposes to soil SOM2 and SOM3 with CO2 loss and
c ... possible leaching of organics.
      if (som1c(SOIL) .gt. 1.e-07) then

c ..... Determine C/E ratios for flows to soil som2
        do 30 iel=1,nelem
          rceto2(iel) = bgdrat(aminrl,varat22,iel)
30      continue

c ..... Compute total C flow out of soil microbes.

c ..... Compute pH effect on decomposition
        pheff = catanf(ph, 4.8, 0.5, 1.14, 0.7)
        pheff = min(pheff, 1.0)
        pheff = max(pheff, 0.0)

        rateDay = dec3(SOIL) * dtm
        call rateconv(rateDay, rate2, t1, t2)
        tcflow = som1c(SOIL) * bgdefac * rate2 * 
     &           cltfac(1) * eftext * anerb * pheff
c ..... where
c .....   som1c(SOIL) = unlabeled and labeled C in soil microbes
c .....                 (som1ci(SOIL,1)+som1ci(SOIL,2))
c .....   dec3(SOIL)  = intrinsic decomposition rate of
c .....                 soil microbes
c .....   cltfac(1)   = cultivation factor for som1
c .....                 (set in cycle)
c .....   eftext      = effect of soil texture on the soil microbe
c .....                 decomposition rate (computed in prelim)

c ..... If soil som1 can decompose to som2, it will also go to som3.
c ..... If it can't go to som2, it can't decompose at all.

c ..... If decomposition can occur,
        if (candec(nelem,aminrl,som1c(SOIL),som1e,2,SOIL,rceto2)) then

c ....... CO2 Loss - Compute and schedule respiration flows
          co2los = tcflow * p1co2(SOIL)
c ....... where 
c .......   p1co2(SOIL) is computed in prelim as a function of fixed parameters 
c .......   p1co2a(SOIL) and p1co2b(SOIL) and soil texture.
          call respir(co2los,2,SOIL,som1c,som1ci,s12c2,resp,
     &                som1e,minerl,gromin,s1mnr,newminrl)

c ....... Decompose Soil SOM1 to SOM3
c ....... The fraction of tcflow that goes to SOM3 is a function of 
c ....... clay content (fps1s3 is computed in prelim).
          cfs1s3 = tcflow * fps1s3 * (1.0 + animpt * 
     &             (1.0 -anerb))
          som12tosom3 = som12tosom3 + cfs1s3

c ....... Partition and schedule C flows by isotope
          call csched(cfs1s3,som1ci(SOIL,LABELD),som1c(SOIL),
     &                som1ci(SOIL,UNLABL),som3ci(UNLABL),
     &                som1ci(SOIL,LABELD),som3ci(LABELD),
     &                1.0,accum)

c ....... Compute and schedule N, P, and S flows and update mineralization 
c ....... accumulators.
          do 40 iel=1,nelem
            rceto3(iel) = bgdrat(aminrl,varat3,iel)
            call esched(cfs1s3,som1c(SOIL),rceto3(iel),
     &                  som1e(SOIL,iel),som3e(iel),
     &                  minerl(SRFC,iel),mnrflo)
            call mnracc(mnrflo,gromin(iel),s1mnr(SOIL,iel))
            if (iel .eq. N) then
              newminrl = newminrl + mnrflo
            endif
40        continue

c ....... Leaching of Organics
c ....... This only occurs when the water flow out of water layer 2
c ....... exceeds a critical value.  Use the same C/N, C/P, and C/S
c ....... ratios as for the flow to SOM3.

          if(amovdly(2) .gt. 0.0) then
            linten = min(1.0-(omlech(3)-amovdly(2))/
     &                   omlech(3), 1.0)
            cleach = tcflow * orglch * linten
c ......... Partition and schedule C flows by isotope
            call csched(cleach,som1ci(SOIL,LABELD),som1c(SOIL),
     &                  som1ci(SOIL,UNLABL),strm5u,
     &                  som1ci(SOIL,LABELD),strm5l,
     &                  1.0,accum)

c ......... Compute and schedule N, P, and S flows and update mineralization 
c ......... accumulators.
            do 50 iel=1,nelem

c ........... Need to use the ratio for som1 for organic leaching     -rm 3/92
c              rceof1 = som1c(SOIL) / som1e(SOIL,iel)
c ........... Dissolved organic matter leaching rates for N and P contents
c ........... were too high, add adjustment to make these rates be
c ........... consistent with observed data, cak - 04/07/03
              if (iel .eq. P) then
                rceof1 = (som1c(SOIL) / som1e(SOIL,iel)) * 35.0
              else
                rceof1 = (som1c(SOIL) / som1e(SOIL,iel)) * 2.0
              endif
              orgflow = cleach / rceof1
              call flow(som1e(SOIL,iel),stream(iel+5),time,
     &                  orgflow)
50          continue

c ....... No leaching this time step
          else
            cleach = 0.
          endif

c ....... Decompose Soil SOM1 to SOM2.
c ....... SOM2 gets what's left of tcflow.
          cfs1s2 = tcflow - co2los - cfs1s3 - cleach
          som12tosom22 = som12tosom22 + cfs1s2

c ....... Partition and schedule C flows by isotope
          call csched(cfs1s2,som1ci(SOIL,LABELD),som1c(SOIL),
     &                som1ci(SOIL,UNLABL),som2ci(SOIL,UNLABL),
     &                som1ci(SOIL,LABELD),som2ci(SOIL,LABELD),
     &                1.0,accum)

c ....... Compute and schedule N, P, and S flows and update mineralization 
c ....... accumulators.
          do 60 iel=1,nelem
            call esched(cfs1s2,som1c(SOIL),rceto2(iel),
     &                  som1e(SOIL,iel),som2e(SOIL,iel),
     &                  minerl(SRFC,iel),mnrflo)
            call mnracc(mnrflo,gromin(iel),s1mnr(SOIL,iel))
            if (iel .eq. N) then
              newminrl = newminrl + mnrflo
            endif
60        continue
        endif
      endif

c ... End of Soil SOM1 decomposition

c ... *****************************************************************

c ... Soil SOM2 decomposes to soil SOM1 and SOM3 with CO2 loss
      if (som2c(SOIL) .gt. 1.e-07) then

c ..... Determine C/E ratios for flows to soil SOM1
        do 70 iel=1,nelem
          rceto1(iel) = bgdrat(aminrl,varat12,iel)
70      continue

c ..... Compute total C flow out of SOM2C

c ..... Compute pH effect on decomposition
        pheff = catanf(ph, 4.0, 0.5, 1.1, 0.7)
        pheff = min(pheff, 1.0)
        pheff = max(pheff, 0.0)

        rateDay = dec5(SOIL) * dtm
        call rateconv(rateDay, rate2, t1, t2)
        tcflow = som2c(SOIL) * bgdefac * rate2 * 
     &           cltfac(2) * anerb * pheff
c ..... where
c .....   dec5(SOIL) = intrinsic decomposition rate of soil som2
c .....   cltfac(2)  = cultivation factor for som2 (set in cycle)

c ..... If som2 can decompose to som1, it will also go to som3.
c ..... If it can't go to som1, it can't decompose at all.

c ..... If decomposition can occur,
        if (candec(nelem,aminrl,som2c(SOIL),som2e,2,SOIL,rceto1)) then

c ....... CO2 loss - Compute and schedule respiration flows
          co2los = tcflow * p2co2(SOIL)

          call respir(co2los,2,SOIL,som2c,som2ci,s22c2,resp,
     &                som2e,minerl,gromin,s2mnr,newminrl)

c ....... Decompose SOM2 to SOM3, SOM3 gets what's left of tcflow.
          cfs2s3 = tcflow * fps2s3 * (1.0 + animpt * (1.0 - anerb))
          som22tosom3 = som22tosom3 + cfs2s3

c ....... Partition and schedule C flows by isotope
          call csched(cfs2s3,som2ci(SOIL,LABELD),som2c(SOIL),
     &                som2ci(SOIL,UNLABL),som3ci(UNLABL),
     &                som2ci(SOIL,LABELD),som3ci(LABELD),
     &                1.0,accum)

c ....... Compute and schedule N, P, and S flows and update mineralization 
c ....... accumulators.
          do 80 iel=1,nelem
            rcetob = bgdrat(aminrl,varat3,iel)
            call esched(cfs2s3,som2c(SOIL),rcetob,
     &                  som2e(SOIL,iel),som3e(iel),
     &                  minerl(SRFC,iel),mnrflo)
            call mnracc(mnrflo,gromin(iel),s2mnr(SOIL,iel))
            if (iel .eq. N) then
              newminrl = newminrl + mnrflo
            endif
80        continue

c ....... Decompose SOM2 to SOM1

          cfs2s1 = tcflow  - co2los - cfs2s3
          som22tosom12 = som22tosom12 + cfs2s1

c ....... Partition and schedule C flows by isotope
          call csched(cfs2s1,som2ci(SOIL,LABELD),som2c(SOIL),
     &                som2ci(SOIL,UNLABL),som1ci(SOIL,UNLABL),
     &                som2ci(SOIL,LABELD),som1ci(SOIL,LABELD),
     &                1.0,accum)

c ....... Compute and schedule N, P, and S flows and update mineralization 
c ....... accumulators.
          do 90 iel=1,nelem
            call esched(cfs2s1,som2c(SOIL),rceto1(iel),
     &                  som2e(SOIL,iel),som1e(SOIL,iel),
     &                  minerl(SRFC,iel),mnrflo)
            call mnracc(mnrflo,gromin(iel),s2mnr(SOIL,iel))
            if (iel .eq. N) then
              newminrl = newminrl + mnrflo
            endif
90        continue
        endif
      endif

c ... End of soil SOM2 Decompositon

c ... *****************************************************************

c ... Surface SOM2 decomposes to surface SOM1 with CO2 loss
c ... This section added by pulliam

      if (som2c(SRFC) .gt. 1.e-07) then

c ..... Determine C/E ratios for flows to surface SOM1
        do 75 iel=1,nelem
          rceto1(iel) = bgdrat(aminrl,varat11,iel)
75      continue

c ..... Compute total C flow out of surface SOM2C

c ..... Compute pH effect on decomposition
        pheff = catanf(ph, 4.0, 0.5, 1.1, 0.7)
        pheff = min(pheff, 1.0)
        pheff = max(pheff, 0.0)

        rateDay = dec5(SRFC) * dtm 
        call rateconv(rateDay, rate2, t1, t2)
        tcflow = som2c(SRFC) * agdefac * rate2 * pheff * mti
c ..... where
c .....   dec5(SRFC) = intrinsic decomposition rate of surface som2

c ..... No flow to som3 from surface.

c ..... If decomposition can occur,
        if (candec(nelem,aminrl,som2c(SRFC),som2e,2,SRFC, rceto1)) then

c ....... CO2 loss - Compute and schedule respiration flows
          co2los = tcflow * p2co2(SRFC)

          call respir(co2los,2,SRFC,som2c,som2ci,s21c2,resp,
     &                som2e,minerl,gromin,s2mnr,newminrl)

c ....... Decompose surface SOM2 to surface SOM1

          cfs2s1 = tcflow  - co2los
          som21tosom11 = som21tosom11 + cfs2s1

c ....... Partition and schedule C flows by isotope
          call csched(cfs2s1,som2ci(SRFC,LABELD),som2c(SRFC),
     &                som2ci(SRFC,UNLABL),som1ci(SRFC,UNLABL),
     &                som2ci(SRFC,LABELD),som1ci(SRFC,LABELD),
     &                1.0,accum)

c ....... Compute and schedule N, P, and S flows and update mineralization 
c ....... accumulators.
          do 95 iel=1,nelem
            call esched(cfs2s1,som2c(SRFC),rceto1(iel),
     &                  som2e(SRFC,iel),som1e(SRFC,iel),
     &                  minerl(SRFC,iel),mnrflo)
            call mnracc(mnrflo,gromin(iel),s2mnr(SRFC,iel))
            if (iel .eq. N) then
              newminrl = newminrl + mnrflo
            endif
95        continue
        endif
      endif

c ... End of surface SOM2 Decompositon

c ... **********************************************************************

c ... SOM3 decomposes to soil SOM1 with CO2 loss.
      if (som3c .gt. 1.e-07) then

c ..... Determine C/E ratios for flows to SOM1.
        do 100 iel=1,nelem
          rceto1(iel) = bgdrat(aminrl,varat12,iel)
100     continue

c ..... Compute pH effect on decomposition
        pheff = catanf(ph, 3.0, 0.5, 1.1, 0.7)
        pheff = min(pheff, 1.0)
        pheff = max(pheff, 0.0)

c ..... Compute total C flow out of SOM3C
        rateDay = dec4 * dtm 
        call rateconv(rateDay, rate2, t1, t2)
        tcflow = som3c * bgdefac * rate2 * cltfac(3) *  
     &           anerb * pheff
c ..... where
c .....   dec4      = intrinsic decomposition rate of som3
c .....   cltfac(3) = cultivation factor for som3 (set in cycle)

c ..... If decomposition can occur,
        if (candec(nelem,aminrl,som3c,som3e,1,1,rceto1)) then

c ....... CO2 loss - Compute and schedule respiration flows.
          cf3co2 = tcflow * p3co2 * anerb

c ....... Changed csrsnk to s3c2 (10/92)
          somc(1) = som3c
          call respir(cf3co2,1,SRFC,somc,som3ci,s3c2,resp,
     &                som3e,minerl,gromin,s3mnr,newminrl)

c ....... Decompose SOM3 to soil SOM1
          cfs3s1 = tcflow - cf3co2
          som3tosom12 = som3tosom12 + cfs3s1

c ....... Partition and schedule C flows by isotope
          call csched(cfs3s1,som3ci(LABELD),som3c,
     &                som3ci(UNLABL),som1ci(SOIL,UNLABL),
     &                som3ci(LABELD),som1ci(SOIL,LABELD),
     &                1.0,accum)

c ....... Compute and schedule N, P, and S flows and update mineralization 
c ....... accumulators.
          do 110 iel=1,nelem
            call esched(cfs3s1,som3c,rceto1(iel),
     &                  som3e(iel),som1e(SOIL,iel),
     &                  minerl(SRFC,iel),mnrflo)
            call mnracc(mnrflo,gromin(iel),s3mnr(iel))
            if (iel .eq. N) then
              newminrl = newminrl + mnrflo
            endif
110       continue
        endif
      endif

c ... End of SOM3 Decomposition

c ... **********************************************************************

c ... Transfer of surface SOM2 to soil SOM2, mixing (added by pulliam)

      if (som2c(SRFC) .gt. 1.0 e-7) then

        if (CURSYS .eq. CRPSYS .or. CURSYS .eq. SAVSYS) then
          mix = cmix
        endif
        if (CURSYS .eq. FORSYS) then
          mix = tmix
        endif

c ..... No pH effect on this mixing flow
        rateDay = mix * dtm
        call rateconv(rateDay, rate2, t1, t2)
        tcflow = som2c(SRFC) * rate2 * agdefac 
        som21tosom22 = som21tosom22 + tcflow

c ..... Partition and schedule C flows by isotope
        call csched(tcflow,som2ci(SRFC,LABELD),som2c(SRFC),
     &              som2ci(SRFC,UNLABL),som2ci(SOIL,UNLABL),
     &              som2ci(SRFC,LABELD),som2ci(SOIL,LABELD),
     &              1.0,accum)

c ..... Compute and schedule N, P, and S flows 
        do 120 iel=1,nelem
          call esched(tcflow,som2c(SRFC),(som2c(SRFC)/som2e(SRFC,iel)),
     &                som2e(SRFC,iel),som2e(SOIL,iel),
     &                minerl(SRFC,iel),mnrflo)
120     continue
      endif
 
c ... **********************************************************************

      return
      end
