
c               Copyright 1993 Colorado State University
c                       All Rights Reserved

      subroutine addomad(omadC, omadN)

      implicit none
      include 'const.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'timvar.inc'
      include 'zztim.inc'

c ... Function Arguments
      double precision omadC, omadN

c ... Organic matter addition
c ... Allow for two types of organic matter addition.  Litter, for
c ... example wheat straw, is added to the structural and metabolic
c ... pools by the partit subroutine.  Partially decomposed organic
c ... matter, for example compost, is added directly into the surface
c ... slow pool (som2c(1))

c ... Local variables
      double precision accum
      double precision k2, fraclabl
      double precision caddm, cadds
      double precision fcadds, fligst
      double precision orglte(MAXIEL)
      integer iel, ii

      accum = 0.0

c ..... Calculate fraction of labeled C as 14C if 14C labeling
        if (labtyp .eq. 1) then
c          fraclabl = ((astlbl / 1000.0) + 1.0) / 1000.0
          k2 = FRAC_C14 * (1.0 + (astlbl / 1000.0))
          fraclabl = k2 / (1.0 + k2)
        endif
c ..... Calculate fraction of labeled C as 13C if 13C labeling
        if (labtyp .eq. 2) then
          fraclabl = astlbl * PEEDEE * 1.0e-03 + PEEDEE
          fraclabl = 1 / (1/fraclabl + 1)
        endif
c ....... Add as fresh plant litter to above-ground litter pools
        if (omadtyp .eq. 1) then
          call partit(astgc,astrec,SRFC,csrsnk,esrsnk,
     &                astlig,astlbl)
        elseif (omadtyp .eq. 2 .or. omadtyp .eq. 4) then
c ....... Add compost to som2c(SRFC)
          call csched(astgc, astlbl, 1.0,
     &                csrsnk(UNLABL), som2ci(SRFC,UNLABL),
     &                csrsnk(LABELD), som2ci(SRFC,LABELD),
     &                1.0, accum)
          do 11 iel = 1, nelem
            call flow(esrsnk(iel), som2e(SRFC,iel), time,
     &                astgc*astrec(iel))
11        continue
        elseif (omadtyp .eq. 3) then
c ....... Add as fresh plant litter to below-ground litter pools
          call partit(astgc,astrec,SOIL,csrsnk,esrsnk,
     &                astlig,astlbl)
        elseif (omadtyp .eq. 5) then
c ....... omadtyp == 5.0 -mdh 5/5/2016
c ....... fcadds is the fraction of plant residue that goes to structural
c ....... (1.0 - fcadds) is the fraction that goes to metabolic
c ....... cadds = amount of structural material added to above-ground litter
c ....... caddm = amount of metabolic material added to above-ground litter
c ....... CURRENTLY 100% metabolic
          fcadds = 0.0
          cadds = astgc*fcadds
          caddm = astgc*(1.0-fcadds)
          call csched(cadds, fraclabl, 1.0,
     &                csrsnk(UNLABL), strcis(SRFC,UNLABL),
     &                csrsnk(LABELD), strcis(SRFC,LABELD),
     &                1.0, accum)
          call csched(caddm, fraclabl, 1.0,
     &                csrsnk(UNLABL), metcis(SRFC,UNLABL),
     &                csrsnk(LABELD), metcis(SRFC,LABELD),
     &                1.0, accum)

          do 40 iel = 1, nelem
c ......... orglte(iel) = amount of omad E in litter addition (gE/m2)
            orglte(iel) = astgc*astrec(iel)
            call flow(esrsnk(iel), struce(SRFC,iel), time,
     &                fcadds*orglte(iel))
            call flow(esrsnk(iel), metabe(SRFC,iel), time,
     &                (1.0-fcadds)*orglte(iel))
40        continue

c ....... Code from partit. 
c ....... Adjust lignin content of structural.
c ....... fligst = fraction of incoming structural residue
c ....... fligst = frlign/(cadds/cpart)
c ....... ATTENTION: What about direct absorption?
          if ((astgc .gt. 0.0) .and. (cadds .gt. 0.0)) then
            fligst = astlig/(cadds/astgc)
          else
            fligst = 0.0
          endif
          fligst = min(1.0, fligst)

c ....... Adjust lignin
          call adjlig(strucc(SRFC),fligst,cadds,strlig(SRFC))

        elseif (omadtyp .eq. 6) then
c ....... omadtyp == 6.0 -mdh 5/14/2016
c ....... fcadds is the fraction of plant residue that goes to structural
c ....... (1.0 - fcadds) is the fraction that goes to metabolic
c ....... cadds = amount of structural material added to below-ground litter
c ....... caddm = amount of metabolic material added to below-ground litter
c ....... CURRENTLY 100% metabolic
          fcadds = 0.0
          cadds = astgc*fcadds
          caddm = astgc*(1.0-fcadds)
          call csched(cadds, fraclabl, 1.0,
     &                csrsnk(UNLABL), strcis(SOIL,UNLABL),
     &                csrsnk(LABELD), strcis(SOIL,LABELD),
     &                1.0, accum)
          call csched(caddm, fraclabl, 1.0,
     &                csrsnk(UNLABL), metcis(SOIL,UNLABL),
     &                csrsnk(LABELD), metcis(SOIL,LABELD),
     &                1.0, accum)

          do 45 iel = 1, nelem
c ......... orglte(iel) = amount of omad E in litter addition (gE/m2)
            orglte(iel) = astgc*astrec(iel)
            call flow(esrsnk(iel), struce(SOIL,iel), time,
     &                fcadds*orglte(iel))
            call flow(esrsnk(iel), metabe(SOIL,iel), time,
     &                (1.0-fcadds)*orglte(iel))
45        continue

c ....... Code from partit. 
c ....... Adjust lignin content of structural.
c ....... fligst = fraction of incoming structural residue
c ....... fligst = frlign/(cadds/cpart)
c ....... ATTENTION: What about direct absorption?
          if ((astgc .gt. 0.0) .and. (cadds .gt. 0.0)) then
            fligst = astlig/(cadds/astgc)
          else
            fligst = 0.0
          endif
          fligst = min(1.0, fligst)

c ....... Adjust lignin
          call adjlig(strucc(SOIL),fligst,cadds,strlig(SOIL))

        else
          write(*,*) 'omadtyp out of range omad.100.'
          write(*,*) 'Must be equal to 1, 2, 3, 4, or 5'
          write(*,*) 'omadtyp currently set to: ', omadtyp
          STOP
        endif

c ..... Update OMAD accumulator output variables, cak - 07/13/2006
        omadC = astgc
        omadac = omadac + astgc
        omadmth(month) = omadmth(month)  + astgc
        omadtot = omadtot + astgc
        do 21 ii = 1, nelem
          omadae(ii) = omadae(ii) +
     &                 (astgc * astrec(ii))
          omadmte(month, ii) = omadmte(month, ii) +
     &                         (astgc * astrec(ii))
          omaetot(ii) = omaetot(ii) +
     &                  (astgc * astrec(ii))
21      continue

      omadN = astgc * astrec(N)

      return
      end
