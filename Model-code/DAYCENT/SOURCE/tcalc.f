
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine tcalc(avgstemp, teff, tfunc)

      implicit none

c ... Argument declarations
      double precision avgstemp, teff(4), tfunc

c ... This function computes the effect of temperature on
c ... decomposition.  It is an exponential function.  Older
c ... versions of Century used a density function.
c ... Created 10/95 - rm
c ...
c ... The temperature effect is now being computed using an
c ... arctangent curve.  CAK - 03/16/01
c ...
c ... Called From:  calcdefac
c ...
c ... Variables
c ...   AVGSTEMP:  weighted average of the average soil temperature in
c ...              the second and third soil layers
c ...   TEFF(1):   "x" location of inflection point
c ...   TEFF(2):   "y" location of inflection point
c ...   TEFF(3):   step size (distance from the maximum point to the
c ...              minimum point)
c ...   TEFF(4):   slope of line at inflection point

c ... Allow C code to call Fortran function
c      !MS$ATTRIBUTES C :: tcalc

c ... Function declarations
      double precision catanf
      external  catanf

c ... Local variables
      double precision normalizer

c ... The normalizer is the value of the numerator at 30 deg C
      normalizer = catanf(30.0, teff(1), teff(2), teff(3), teff(4))

      tfunc = max(0.01,
     &            catanf(avgstemp, teff(1), teff(2), teff(3), teff(4)) /
     &            normalizer)

c     write(*,*) 'normalizer = ', normalizer
c     write(*,*) 'tfunc = ', tfunc
      return
      end
