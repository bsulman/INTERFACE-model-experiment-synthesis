

c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine addfert(nelem, month, fertN)

      implicit none
      include 'const.inc'
      include 'dovars.inc'
      include 'fertil.inc'
      include 'npool.inc'
      include 'plot1.inc'
      include 'plot2.inc'

c ... FORMAL PARAMETERS
      integer nelem
      integer month
      double precision fertN

c ... Local Variables
      integer iel, clyr
      character subname*10

      subname = 'addfert   '

c ... Fertilization option

      fertN = 0

        do 60 iel = 1, nelem
          if (iel .eq. N) then
            clyr = SRFC
            esrsnk(iel) = esrsnk(iel) - feramt(iel)
            call update_npool(clyr, feramt(iel),
     &                        frac_nh4_fert, frac_no3_fert,
     &                        ammonium, nitrate, subname)
            minerl(SRFC,iel) = minerl(SRFC,iel) + feramt(iel)
            fertot(iel) = fertot(iel) + feramt(iel)
            fertac(iel) = fertac(iel) + feramt(iel)
            fertmth(month,iel) = fertmth(month,iel) + feramt(iel)
            fertN = feramt(N)
          else
            esrsnk(iel) = esrsnk(iel) - feramt(iel)
            minerl(SRFC,iel) = minerl(SRFC,iel) + feramt(iel)
            fertot(iel) = fertot(iel) + feramt(iel)
            fertac(iel) = fertac(iel) + feramt(iel)
            fertmth(month,iel) = fertmth(month,iel) + feramt(iel)
          endif
60      continue
        fertcnt = 1
        nreduce = ninhib

      return
      end

