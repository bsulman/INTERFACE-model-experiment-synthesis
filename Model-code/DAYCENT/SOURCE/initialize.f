
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine initialize(ext)

      implicit none
      include 'chrvar.inc'
      include 'const.inc'
      include 'dovars.inc'
      include 'fertil.inc'
      include 'npool.inc'
      include 'param.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'seq.inc'
      include 'timvar.inc'
      include 'wth.inc'
      include 'zztim.inc'

c ... Argument declarations
      logical ext

c ... This subroutine initializes several common block variables that were
c ... being used in the code without initialization, - cak - 09/26/00

c ... Local variables
      integer ii, jj, kk

c ... Initialize weather labels
      wlabel(1) = 'prec'
      wlabel(2) = 'tmin'
      wlabel(3) = 'tmax'

c ... Set to cause initialization of annual production accumulators to be done
c ... in inprac subroutine
      month = 1
      dofone = .true.
      dofrst = .true.
      doplnt = .true.
c ... Used in schedl subroutine, set to .true. so that crpgrw, msplt, and
c ... forgrw are initialized
      dolast = .true.
      doflst = .true.
      call inprac(CRPSYS)
      call inprac(FORSYS)
c ... If this is an extend simulation call inprac a second time for
c ... each system at the start of the simulation so that the production
c ... output does not carry over from the previous simulation
      if (ext) then
        call inprac(CRPSYS)
        call inprac(FORSYS)
      endif

c ... Used in cropin when reading initial crop which passes it to cmplig
c ... without having been initialized
      if (initcp .ne. ' ' .and. initre .ne. ' ') then
        cursys = SAVSYS 
      else if (initcp .ne. ' ') then
        cursys = CRPSYS
      else
        cursys = FORSYS
      endif

c ... Used in partit, flow, csched, calciv, and flowup without having been
c ... initialized
      time = strtyr

c ... Initialize new variables added to fertil common block
      nreduce = 1.0

c ... The co2*(*) and cltfac(*) parameters from the plot1 and plot2
c ... common blocks should not be inititalized on an extend
      if (.not. ext) then
c ..... These variables are set only if the proper conditions are met,
c ..... initialize with default values
        cltfac(1) = 1.0
        cltfac(2) = 1.0
        cltfac(3) = 1.0
        cltfac(4) = 1.0

c ..... Initialize these variables for first print of plot commons
        do 80 ii = 1, 2
          co2crs(ii) = 1.0
          co2cpr(ii) = 1.0
          co2ctr(ii) = 1.0
          do 85 jj = 1, 2
            do 90 kk = 1, MAXIEL
              co2cce(ii,jj,kk) = 1.0
90          continue
85        continue
80      continue
      endif

c ... Initialize the frac_nh4_fert and frac_no3_fert npool commons for
c ... use by the detiv subroutine, cak - 04/16/2007
      frac_nh4_fert = 0.5
      frac_no3_fert = 0.5

      if (ext) then
        do 230 ii = 1, MAXIEL
          do 220 jj = 1, MAXLYR
            iminerl(jj,ii) = minerl(jj,ii)
220       continue
230     continue
      endif

      return
      end
