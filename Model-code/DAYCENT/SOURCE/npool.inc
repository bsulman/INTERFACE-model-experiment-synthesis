c************************************************************************
c  npool.inc
c
c  Melannie D. Hartman
c
c************************************************************************

      common /npool/ nitrate(21), ammonium, frac_nh4_fert,
     &               frac_no3_fert, texture, iminerl(10,3)

      double precision nitrate, ammonium, frac_nh4_fert,
     &                 frac_no3_fert, iminerl
      integer          texture

      save /npool/
