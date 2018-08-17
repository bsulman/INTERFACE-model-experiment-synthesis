
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


C ... PSCHEM.F

      subroutine pschem(dtm, t1, t2)

      implicit none
      include 'const.inc'
      include 'npool.inc'
      include 'param.inc'
      include 'parfx.inc'
      include 'plot1.inc'
      include 'site.inc'
      include 'zztim.inc'

c ... Argument declarations
      double precision dtm, t1, t2

c ... This subroutine calculates the Phosophorus and Sulfur chemistry
c ... - decomposition flows.  These calculations were removed from
c ... the DECOMP subroutine, and slightly modified to include the
c ... calculation for the fraction of mineral P in solution.  Also
c ... removed the flow from secondary to parent material.
c ...
c ... Date:            6/91
c ... Coded by:        McKeown
c ...
c ... Called From:       SIMSOM
c ... Calls:             FLOW
c ...                    FLOW_DOUBLE
c ...                    FLOW_DOUBLE_IN
c ...                    FLOW_DOUBLE_OUT
c ...                    FSFUNC
c ...
c ... Local Variables:   FMNSEC        - flow from mineral to secondary
c ...                    FOCSEC        - flow from occluded P to secondary
c ...                    FPARNT        - flow from parent to mineral
c ...                    FSECND        - flow from secondary to mineral
c ...                    FSECOC        - flow from secondary to occluded P
c ...                    FSOL          - fraction of mineral P in solution
c ...
c ... ********************************************************************

c ... Function declarations
      double precision fsfunc
      external  fsfunc

c ... Local variables
      integer          iel, lyr
      double precision fparnt, fsol
      double precision dely, delx, xslope, yint
      double precision fmnsec, focsec, fsecnd, fsecoc
      double precision rateDay, rate2
      character        subname*10

      subname = 'pschem    '

c ... Code added to modify the intercept for the texture equation of
c ... secondary P dependant upon pH input and include the effect of
c ... texture on psecmn(2), copied from prelim subroutine, cak 10/17/05
      if (ph .le. phesp(1)) then
        texesp(2) = phesp(2)
      else if (ph .ge. phesp(3)) then
        texesp(2) = phesp(4)
      else
        dely = phesp(4) - phesp(2)
        delx = phesp(3) - phesp(1)
        xslope = dely / delx
        yint = phesp(2) - (xslope*phesp(1))
        texesp(2) = (xslope*ph) + yint
      endif
      if (texesp(1) .eq. 1.0) then
c ..... Calculate psecmn(2)
c ..... Include effect of texture
c ..... Note that this code changes the value of a 'fixed' parameter
c ..... (psecmn(2))
        psecmn(2) = 12.0 * (texesp(2) + texesp(3) * sand)
      endif

c ... For Phosphorus and Sulfur...
      do 20 iel = 2, nelem

c ..... Determine the fraction of mineral P in solution
        if (iel .eq. P) then
          fsol = fsfunc(minerl(SRFC,iel), pslsrb, sorpmx)
        else
          fsol = 1.0
        endif

c ..... Flow from parent material to mineral compartment.
c ..... Soil Texture may affect weathering of Phosophorus.

c ..... This calculation is actually done in prelim, it is
c ..... shown here for clarity.
c .....   if ((iel .eq. P) .and. (texepp(1) .eq. 1.0)) then
c .....     Include effect of texture
c .....     textur = clay + silt
c .....     Weathering factor should be per year
c .....     wfact = 12.0 * catanf(textur, texepp(2), texepp(3),
c ..... +                         texepp(4), texepp(5))
c .....   else
c .....     wfact = pparmn(iel)
c .....   endif

        rateDay = pparmn(iel) * dtm
        call rateconv(rateDay, rate2, t1, t2)
        fparnt = parent(iel) * rate2 * bgdefac 
        call flow(parent(iel), minerl(1,iel), time, fparnt)

c ..... Flow from secondary to mineral compartment.
c ..... Soil texture may affect mineralization of secondary P.

c ..... This calculation is actually done in prelim, it is
c ..... shown here for clarity.
c .....   if ((iel .eq. P) .and. (texesp(1) .eq. 1.0)) then
c .....     Include effect of texture
c .....     wfact = 12.0 * (texesp(2) + texesp(3) * sand)
c .....   else
c .....     wfact = psecmn(iel)
c .....   endif

c ..... Use double precision variables for computing the flows to/from
c ..... secondary P to occluded P, cak - 03/20/02
        rateDay = psecmn(iel) * dtm
        call rateconv(rateDay, rate2, t1, t2)
        fsecnd = secndy(iel) * rate2 * bgdefac 
        call flow(secndy(iel), minerl(1,iel), time, fsecnd)

c ..... Flow from mineral to secondary
        do 10 lyr = 1, nlayer
          rateDay = pmnsec(iel) * dtm
          call rateconv(rateDay, rate2, t1, t2)
          fmnsec = minerl(lyr,iel) * (1 - fsol) *
     &             rate2 * bgdefac
          call flow(minerl(lyr,iel), secndy(iel),
     &                     time, fmnsec)
10      continue
20    continue

c ... Flow from secondary Phosophorus to occluded Phosophorus.
      rateDay = psecoc1 * dtm
      call rateconv(rateDay, rate2, t1, t2)
      fsecoc = secndy(P) * rate2 * bgdefac 
      call flow(secndy(P), occlud, time, fsecoc)

c ... Flow from occluded Phosophorus to secondary Phosophorus, cak - 03/20/02
      rateDay = psecoc2 * dtm
      call rateconv(rateDay, rate2, t1, t2)
      focsec = occlud * rate2 * bgdefac 
      call flow(occlud, secndy(P), time, focsec)

      return
      end
