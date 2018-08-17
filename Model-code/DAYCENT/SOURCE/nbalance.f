
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine nbalance(time, curday, omadN, fertN, 
     &              orgNlch, inorgNlch,
     &              NOflux, Nn2oflux, Dn2oflux, Dn2flux, 
     &              NOabsorp_grass, NOabsorp_tree, 
     &              atmosNdep, nonSymSoilNfix, totNstate)

      implicit none
      include 'const.inc'
      include 'param.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'

c ... Argument declarations
      integer curday
      double precision time
      double precision omadN, fertN
      double precision orgNlch, inorgNlch
      double precision atmosNdep, nonSymSoilNfix
      double precision NOflux, Nn2oflux, Dn2oflux, Dn2flux 
      double precision NOabsorp_grass, NOabsorp_tree
      double precision totNstate

c ... Local variables
      integer iel, lyr 
      double precision Nlosses, Ngains, totNstate1, balance

c ... Save the previous sum of N in state variables

      totNstate1 = totNstate

c ... atmosNdep is atmospheric N deposition (gN/m^2)
c ... nonSymSoilNfix is non-symbiotic soil N fixation (gN/m^2)
c ... orgNlch is organic N leaching (gN/m^2)
c ... inorgNlch is inorganic N leaching (gN/m^2)
c ... Trace gas fluxes: NOflux, Nn2oflux, Dn2oflux, Dn2flux
c ...                   NOabsorp_grass, NOabsorp_tree (gN/m^2)

      Ngains = atmosNdep + nonSymSoilNfix + omadN + fertN
      Nlosses = orgNlch + inorgNlch +
     &          NOflux + Nn2oflux + Dn2oflux + Dn2flux + 
     &          NOabsorp_grass + NOabsorp_tree 

c ... Compute the total C in the system

      do 10 iel = 1, nelem
        somse(iel) = som1e(SOIL,iel) + som2e(SOIL,iel) + som3e(iel)
        somte(iel) = somse(iel) + struce(SOIL,iel) + metabe(SOIL,iel)
        woode(iel) = wood1e(iel) + wood2e(iel) + wood3e(iel)
        frste(iel) = rleave(iel) + frootej(iel) + frootem(iel) +
     &               fbrche(iel) + rlwode(iel) + croote(iel)
        fsyse(iel) = somte(iel) + woode(iel) + frste(iel) +
     &               struce(SRFC,iel) + metabe(SRFC,iel) +
     &               som1e(SRFC,iel) + som2e(SRFC,iel)
        totsyse(iel) = fsyse(iel) + aglive(iel) + bglivej(iel) +
     &                 bglivem(iel) + stdede(iel)
10    continue

c ... Calculate tminrl
      do 50 iel = 1, nelem
        tminrl(iel) = 0.0
        do 40 lyr = 1, nlayer
          if (minerl(lyr,iel).gt.0.0) then
            tminrl(iel) = tminrl(iel) + minerl(lyr,iel)
          endif
40      continue
50    continue

      do 60 iel = 1, nelem
        totale(iel) = totsyse(iel) + tminrl(iel) + 
     &                minerl(nlayer+1,iel) + parent(iel) + secndy(iel) +
     &                crpstg(iel) + forstg(iel)
        if (iel .eq. P) then
          totale(iel) = totale(iel) + occlud
        endif
60    continue

      totNstate = totale(N)

      balance = totNstate1 - totNstate + Ngains - Nlosses

      if (totNstate1 .gt. 0.0) then
          write(124,100) time, curday, balance, totNstate1, totNstate,
     &                   Ngains, Nlosses, atmosNdep, nonSymSoilNfix,
     &                   omadN, fertN,
     &                   orgNlch, inorgNlch, NOflux, Nn2oflux, Dn2oflux, 
     &                   Dn2flux, NOabsorp_grass, NOabsorp_tree 
      end if

100   format(f10.4,',',i4,',',17(f12.4,','))

      return
      end
