
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


C ... LITDEC.F

      subroutine litdec(dtm, t1, t2, newminrl, soilsrad)

      implicit none
      include 'cflows.inc'
      include 'comput.inc'
      include 'const.inc'
      include 'param.inc'
      include 'parfx.inc'
      include 'plot1.inc'
      include 'plot2.inc'

c ... Argument declarations
      double precision dtm, soilsrad
      double precision t1, t2
      double precision newminrl

c ... Litter Decomposition
c ... Decompose structural and metabolic material for surface and soil.

c ... Function declarations
      logical  candec
      double precision agdrat, bgdrat, catanf
      external agdrat, bgdrat, candec, catanf

c ... Local variables
      integer  iel, lyr
      double precision accum(ISOS), biocnv, cfmes1, co2los, mnrflo, 
     &         pheff, rceto1(3), tcflow
      double precision a, b, mdr, x1, x2
      double precision rateDay, rate2

c ... Factor to convert C to biomass is 2.5 for everything but wood.
      parameter (biocnv = 2.5)

      accum(LABELD) = 0.0
      accum(UNLABL) = 0.0

c ... Calculate metabolic decomposition reducer (mdr) based on incoming solar
c ... radiation, cak - 12/08/2009
      if (soilsrad .le. 0.0) then
        mdr = 1.0
      elseif (soilsrad .gt. 30000) then
        mdr = 0.2
      else
        a = 1.0
        b = 0.2
        x1 = 0.0
        x2 = 30000.0
        mdr = abs((b - a) / (x2 - x1)) * (x2 - soilsrad) + b
      endif

c ... Surface STRUCTURAL Material
      if(strucc(SRFC) .gt. 1.e-07) then

c ..... Compute pH effect on decomposition
        pheff = catanf(ph, 4.0, 0.5, 1.1, 0.7)
        pheff = min(pheff, 1.0)
        pheff = max(pheff, 0.0)

c ..... Compute total C flow out of structural in layer SRFC
        rateDay = dec1(SRFC) * dtm
        call rateconv(rateDay, rate2, t1, t2)
        tcflow = min(strucc(SRFC),strmax(SRFC)) * agdefac *
     &               rate2 * exp(-pligst(SRFC)*strlig(SRFC))
     &               * pheff
c ..... where
c .....   tcflow       = grams C
c .....   strucc(SRFC) = the current value for total structural C
c .....                  (grams C) in layer SRFC
c .....   strmax(SRFC) = the maximum amount of structural C in 
c .....                  layer SRFC that will decompose (grams C).
c .....   agdefac      = the aboveground decomposition factor based on
c .....                  water and temperature computed in prelim and
c .....                  in cycle
c .....   dec1(SRFC)   = the intrinsic decomposition rate of 
c .....                  structural C (a fixed parameter)
c .....   pligst       = a fixed parameter that represents the effect 
c .....                  of lignin-to-structural-ratio on structural 
c .....                  decomposition
c .....   strlig(SRFC) = the lignin content of surface structural 
c .....                  residue (grams lignin/grams biomass)
c .....   dtm          = the time step in years
c .....   pheff        = pH effect on decomposition

c ..... Decompose structural into som1 and som2 with CO2 loss.
c ..... Changed csrsnk to st1c2 (10/92)

        call declig(aminrl,strlig(SRFC),SRFC,nelem,2,ps1co2,ratnew,
     &              rsplig,tcflow,strucc,st1c2,strcis,struce,
     &              gromin,minerl,strmnr,resp,som1ci,som1e,som2ci,
     &              som2e,struc1tosom11,struc1tosom21,newminrl)

      endif

c ... Soil STRUCTURAL Material
      if (strucc(SOIL) .gt. 1.e-07) then

c ..... Compute pH effect on decomposition
        pheff = catanf(ph, 4.0, 0.5, 1.1, 0.7)
        pheff = min(pheff, 1.0)
        pheff = max(pheff, 0.0)

c ..... Compute total C flow out of structural in layer SOIL
        rateDay = dec1(SOIL) * dtm
        call rateconv(rateDay, rate2, t1, t2)
        tcflow = min(strucc(SOIL),strmax(SOIL)) * bgdefac *
     &               rate2 * exp(-pligst(SOIL)*strlig(SOIL))
     &               * cltfac(4) * anerb * pheff
c ..... where
c .....   tcflow       = grams C
c .....   strucc(SOIL) = the current value for total structural C
c .....                  (grams C) in layer SOIL
c .....   strmax(SOIL) = the maximum amount of structural C in 
c .....                  layer SOIL that will decompose (grams C).
c .....   bgdefac      = the belowground decomposition factor based on
c .....                  water and temperature computed in prelim and
c .....                  in cycle
c .....   dec1(SOIL)   = the intrinsic decomposition rate of 
c .....                  structural C (a fixed parameter)
c .....   pligst       = a fixed parameter that represents the effect 
c .....                  of lignin-to-structural-ratio on structural 
c .....                  decomposition
c .....   strlig(SOIL) = the lignin content of soil structural 
c .....                  residue (grams lignin/grams biomass)
c .....   cltfac(4)    = the cultivation factor for soil
c .....                  structural material (set in cycle)
c .....   anerb        = impact of soil anaerobic conditions on
c .....                  decomposition
c .....   dtm          = the time step in years 
c .....   pheff        = pH effect on decomposition

c ..... Decompose structural into som1 and som2 with CO2 loss.
c ..... Changed csrsnk to st2c2 (10/92)

        call declig(aminrl,strlig(SOIL),SOIL,nelem,2,ps1co2,ratnew,
     &              rsplig,tcflow,strucc,st2c2,strcis,struce,
     &              gromin,minerl,strmnr,resp,som1ci,som1e,som2ci,
     &              som2e,struc2tosom12,struc2tosom22,newminrl)

      endif
c ... End of Structural Decomposition

c ... METABOLIC Material

c ... Process each layer
      do 30 lyr = SRFC, SOIL

        if (metabc(lyr) .gt. 1.e-07) then

c ....... Determine C/E ratios for flows to SOM1
          do 10 iel=1,nelem

c ......... Compute ratios for surface metabolic residue
            if (lyr .eq. SRFC) then
              rceto1(iel) = agdrat(aminrl,varat11,iel)

c ......... Compute ratios for soil metabolic residue
            else
              rceto1(iel) = bgdrat(aminrl,varat12,iel)
            endif
10        continue

c ....... Compute pH effect on decomposition
          pheff = catanf(ph, 4.8, 0.5, 1.14, 0.7)
          pheff = min(pheff, 1.0)
          pheff = max(pheff, 0.0)

c ....... Compute total C flow out of metabolic in layer lyr
          rateDay = dec2(lyr) * dtm
          call rateconv(rateDay, rate2, t1, t2)
          if (lyr .eq. SRFC) then
            tcflow = metabc(lyr) * agdefac * rate2 * pheff *
     &               mdr
          else
            tcflow = metabc(lyr) * bgdefac * rate2 * pheff
          endif

c ....... Added impact of soil anerobic conditions -rm 12/91
          if (lyr .eq. SOIL) then
            tcflow = tcflow * anerb
          endif
c ....... where:
c .......   tcflow      = grams C
c .......   metabc(lyr) = the current value for total metabolic C
c .......                 (grams C) in layer lyr
c .......   agdefac     = the aboveground decomposition factor
c .......   bgdefac     = the belowground decomposition factor
c .......   dec2(lyr)   = the intrinsic decomposition rate of 
c .......                 metabolic C
c .......   dtm         = the time step in years
c .......   pheff       = pH effect on decomposition
c .......   mdr         = metabolic decomposition reducer
c .......   anerb       = impact of soil anaerobic conditions on
c .......                 decomposition

c ....... Make sure metab does not go negative.
          if (tcflow .gt. metabc(lyr)) then
            tcflow = metabc(lyr)
          endif

c ....... If decomposition can occur,
          if (candec(nelem,aminrl,metabc(lyr),metabe,2,lyr,rceto1))
     &        then

c ......... CO2 loss
            co2los = tcflow * pmco2(lyr)

c ......... Changed csrsnk to mt1c2, mt2c2 (10/92)
            if (lyr .eq. SRFC) then
              call respir(co2los,2,lyr,metabc,metcis,mt1c2,
     &                    resp,metabe,minerl,gromin,metmnr,newminrl)
            else
              call respir(co2los,2,lyr,metabc,metcis,mt2c2,
     &                    resp,metabe,minerl,gromin,metmnr,newminrl)
            endif

c ......... Decompose metabolic into som1
            cfmes1 = tcflow - co2los
            if (lyr .eq. 1) then
              metc1tosom11 = cfmes1
            else
              metc2tosom12 = cfmes1
            endif

c ......... Partition and schedule C flows by isotope
            call csched (cfmes1,metcis(lyr,LABELD),metabc(lyr),
     &                   metcis(lyr,UNLABL),som1ci(lyr,UNLABL),
     &                   metcis(lyr,LABELD),som1ci(lyr,LABELD),
     &                   1.0,accum)

c ......... Compute and schedule N, P, and S flows and update
c ......... mineralization accumulators.
            do 20 iel = 1, nelem

              call esched(cfmes1,metabc(lyr),rceto1(iel),
     &                    metabe(lyr,iel),som1e(lyr,iel),
     &                    minerl(SRFC,iel),mnrflo)
              call mnracc(mnrflo,gromin(iel),metmnr(lyr,iel))
              if (iel .eq. N) then
                newminrl = newminrl + mnrflo
              endif
20          continue
          endif
        endif

c ..... End of Metabolic Decomposition

c ... Next layer
30    continue

      return
      end
