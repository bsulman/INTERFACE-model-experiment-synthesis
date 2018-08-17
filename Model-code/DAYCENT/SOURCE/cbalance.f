
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine cbalance(time, curday, omadC, hrespCO2, co2losUV,
     &                     CH4, totCstate)

      implicit none
      include 'const.inc'
      include 'param.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'

c ... Argument declarations
      integer curday
      double precision time
      double precision omadC
      double precision hrespCO2, co2losUV
      double precision CH4
      double precision totCstate

c ... Compute the total C in the system

c ... Local variables
      double precision Closses, Cgains, totCstate1, balance

c ... Save the previous sum of C in state variables

      totCstate1 = totCstate

c ... omadC - C from organic matter addition (gC/m^2)
c ... hrespCO2 	- sum the co2 heterotrophic respiration losses 
c ... co2losUV 	- CO2 loss from photodegradation of litter
c ... DOC = strm5u + strm5l
c ... Note: the source of C for CH4 is not transferred to/from any
c ...   C pool, therefore it is not included as a gain or loss

      Cgains = omadC
      Closses = hrespCO2 + co2losUV + strm5u + strm5l 

      somsc = som1c(SOIL) + som2c(SOIL) + som3c
      somtc = somsc + strucc(SOIL) + metabc(SOIL)
      woodc = wood1c + wood2c + wood3c
      frstc = rleavc + frootcj + frootcm + fbrchc + rlwodc + crootc
      fsysc = somtc + woodc + frstc + strucc(SRFC) + metabc(SRFC) +
     &        som1c(SRFC) + som2c(SRFC)
      totsysc = fsysc + aglivc + bglivcj + bglivcm + stdedc +
     &          carbostg(CRPSYS,UNLABL) + carbostg(CRPSYS,LABELD) +
     &          carbostg(FORSYS,UNLABL) + carbostg(FORSYS,LABELD)

      totCstate = totsysc

      balance = totCstate1 - totCstate + Cgains - Closses

      if (totCstate1 .gt. 0.0) then
          write(123,100)time, curday, balance, totCstate1, totCstate,
     &                  Cgains, Closses, omadC, hrespCO2, co2losUV, 
     &                  strm5u, strm5l, CH4
      end if

100   format(f10.4,',',i4,',',11(f12.4,','))

      return
      end
