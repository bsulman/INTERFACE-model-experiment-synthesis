
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


C ... ADDLIT

      subroutine addlit()

      implicit none
      include 'const.inc'
      include 'forrem.inc'
      include 'isovar.inc'
      include 'param.inc'
      include 'parfs.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'site.inc'
      include 'timvar.inc'
      include 'zztim.inc'

c ... Death of leaves, fine branches, large wood, fine roots, and coarse roots.

c ... Local variables
      integer iel
      double precision accum(ISOS), ctodie, etodie, fr14, recres(MAXIEL)

      accum(LABELD) = 0.0
      accum(UNLABL) = 0.0

      recres(N) = 0.04
      recres(P) = 0.004
      recres(S) = 0.004
      fr14 = 0.0

c ... Leave litter

      ctodie = 0.1
      call partit(ctodie, recres, SRFC, rlvcis, rleave, 
     &            wdlig(LEAF), fr14)


c ... Fine Root litter

      ctodie = 0.1
      call partit(ctodie, recres, SOIL, frtcisj, frootej, 
     &            wdlig(FROOTJ), fr14)


c ... Fine Branches, Large Wood, and Coarse Roots go to the dead wood
c ... compartments: WOOD1, WOOD2, WOOD3

      recres(N) = 0.01
      recres(P) = 0.001
      recres(S) = 0.001

c ... Death of fine branches

      ctodie = 0.1
      call csched(ctodie, fbrcis(LABELD), fbrchc,
     &            fbrcis(UNLABL), wd1cis(UNLABL),
     &            fbrcis(LABELD), wd1cis(LABELD),
     &            1.0, accum)

      do 30 iel = 1, nelem
        etodie = ctodie * recres(iel)
        call flow(fbrche(iel), wood1e(iel), time, etodie)
30    continue


c ... Death of large wood

      ctodie = 0.1
      call csched(ctodie, rlwcis(LABELD), rlwodc,
     &            rlwcis(UNLABL), wd2cis(UNLABL),
     &            rlwcis(LABELD), wd2cis(LABELD),
     &            1.0, accum)

      do 40 iel = 1, nelem
        etodie = ctodie * recres(iel)
        call flow(rlwode(iel), wood2e(iel), time, etodie)
40    continue


c ... Death of coarse roots

      ctodie = 0.1
      call csched(ctodie, crtcis(LABELD), crootc,
     &            crtcis(UNLABL), wd3cis(UNLABL),
     &            crtcis(LABELD), wd3cis(LABELD),
     &            1.0, accum)

      do 50 iel = 1, nelem
        etodie = ctodie * recres(iel)
        call flow(croote(iel), wood3e(iel), time, etodie)
50    continue

      return
      end
