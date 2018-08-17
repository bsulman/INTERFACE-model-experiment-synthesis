

c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine rateconv(rate1, rate2, t1, t2)

      implicit none

c ...

c ... FORMAL PARAMETERS
      double precision rate1, rate2, t1, t2

c ... Local Variables
      double precision kcont

      if (rate1 .gt. 1.0) then
          rate2 = rate1 * t2 / t1
      else
          kcont = -alog(1.0 - rate1)
          rate2 = 1.0 - exp(-kcont * t2 / t1)
      endif

      return
      end

