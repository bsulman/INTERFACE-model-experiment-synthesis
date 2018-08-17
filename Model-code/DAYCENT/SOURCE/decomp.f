
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


C ... DECOMP.F

      subroutine decomp(dtm,t1,t2,decsys,amovdly,newminrl,
     &                  soilsrad,rpeff)

      implicit none
      include 'comput.inc'
      include 'const.inc'
      include 'param.inc'
      include 'parfx.inc'
      include 'plot1.inc'

c ... Argument declarations
      double precision dtm
      double precision t1, t2
      integer   decsys
      double precision amovdly(MAXLYR), soilsrad, rpeff
      double precision newminrl

c ... Decomposition Submodel (rewritten by vek 04/91)


c ... Function declarations
      double precision agdrat
      external agdrat

c ... Local variables
      integer   iel

c ... ratnew(iel,1) - the C/E ratio for new material created when a
c ...                 lignin component decomposes to som1.
c ... ratnew(iel,2) - the C/E ratio for new material created when a
c ...                 lignin component decomposes to som2.

c ... Determine C/E ratios for flows from structural material to
c ... surface som1 and surface som2
c      do 30 iel=1,nelem
c        ratnew(iel,SRFC) = agdrat(aminrl,varat11,iel)
c        ratnew(iel,SOIL) = agdrat(aminrl,varat21,iel)
c30    continue
c ... The second index of ratnew(iel,*) should be 1 or 2 (not SRFC or SOIL) 
c ... for som1c or som2c. -MDH 1/23/2012
      do 30 iel=1,nelem
        ratnew(iel,1) = agdrat(aminrl,varat11,iel)
        ratnew(iel,2) = agdrat(aminrl,varat21,iel)
30    continue

c                   ********** LITTER **********
c ... Decompose structural and metabolic components for surface and soil.
      call litdec(dtm, t1, t2, newminrl, soilsrad)

c                   *********** WOOD ***********
c ... If the system is a forest or savanna...
c ... Decompose dead fine branches, large wood, and coarse roots.
c ... Dead fine roots are in the soil structural compartment.
      if (decsys .eq. FORSYS) then
        call woodec(dtm, t1, t2, newminrl)
      endif

c                 ***** SOIL ORGANIC MATTER *****
c ... Decompose som1 and som2 (surface and soil) and som3.
c ... Added amovdly parameter for daily version. -mdh 10/10/94
      call somdec(amovdly, dtm, t1, t2, newminrl, rpeff, soilsrad)

      return
      end
