
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine detiv

      implicit none
      include 'chrvar.inc'
      include 'const.inc'
      include 'dovars.inc'
      include 'fertil.inc'
      include 'jday.inc'
      include 'npool.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfs.inc'
      include 'parfx.inc'
      include 'plot1.inc'
      include 'plot3.inc'
      include 'seq.inc'
      include 'site.inc'
      include 't0par.inc'
      include 'timvar.inc'
      include 'wth.inc'
      include 'zztim.inc'

c ... Determine name of schedule file, which contains the
c ... name of the site file, values of timing variables,
c ... and order of events

c ... Function declarations
      integer getlen, iargc

      common /libpath/filpath
      character*100 filpath
      character*50 sitnam
      character*10 subname

c ... Local variables
      integer          clen, ii, nargs
      character*100    extflag, newbin, oldbin, schnam
      logical          ext, goahead
      character*100    iname
      integer*2        status
      character*80     string
      integer          iel, numlyrs
      double precision swcinit(SWMAXLYR)
      character*100    cdifile

      integer          wthinput, wthstart

c ... Add new routine to do a "brute force" initialization of all common block
c ... variables, cak - 06/04/02
      call default

      subname = 'detiv     '

c ... Initialize potential command line arguments
      ext = .false.
      schnam = ' '
      newbin = ' '
      oldbin = ' '
      filpath = ' '
      cdifile = ' '

c ... VAX NOTE: Need to get information interactively from user

c ... Get command line arguments
      nargs = iargc()
      if (nargs .eq. 0) then
        call message(' ')
        call message('              DAYCENT SOIL ORGANIC MATTER')
        call message('                 AND TRACEGAS MODEL')
        call message('                 STAND-ALONE VERSION')
        call message('                     DayCent_UV')
        call message('                     10/31/2011')
        call message(' ')
        call message('  Invalid command line arguments were supplied.')
        call message('  To run DayCent, please supply these arguments')
        call message('  as needed:')
        call message('     -l   directory to search for library files')
        call message('     -n   name of binary output file (no .bin)')
        call message('     -s   name of schedule file (no .sch)')
        call message('     -e   name of old binary file to extend')
        call message('          from, if extending (no .bin)')
        call message('     -c   name of cdi file')
        call message('  Example:')
        call message('  daycent -l /lib -e oldata -s schnam '//
     &               '-n newoutput')
        call message(' ')
        STOP 'Execution error.'
      endif

c ... Process command line arguments
      ii = 1
10    if (ii .lt. nargs) then
        call getarg(ii, extflag)
        ii = ii + 1
        call getarg(ii, iname)
        ii = ii + 1
        clen = getlen(iname)

        if (extflag .eq. '-s') then
          schnam = iname(1:clen)
          if (index(schnam,'.sch').eq.0) schnam(clen+1:clen+4) = '.sch'
          inquire(file=schnam,exist=goahead)
          if (.not. goahead) then
            call message(' ')
            call message('   The schedule file could not be read.')
            call message(' ')
            STOP 'Execution error.'
          endif

        else if (extflag .eq. '-n') then
          newbin = iname(1:clen)
          if (index(newbin,'.bin').eq.0) newbin(clen+1:clen+4) = '.bin'
          inquire(file=newbin,exist=goahead)
          if (goahead) then
            call message(' ')
            call message('   The new binary file already exists.')
            call message(' ')
            STOP 'Execution error.'
          endif

        else if (extflag .eq. '-e') then
          ext = .true.
          oldbin = iname
          if (index(oldbin,'.bin').eq.0) oldbin(clen+1:clen+4) = '.bin'
          inquire(file=oldbin,exist=goahead)
          if (.not. goahead) then
            call message(' ')
            call message('   The old binary file could not be read.')
            call message(' ')
            STOP 'Execution error.'
          endif

        else if (extflag .eq. '-l') then
          filpath = iname(1:clen)

        else if (extflag .eq. '-c') then
          cdifile = iname(1:clen)
          write(*,*) 'cdifile = ', cdifile

        else
          call message('   Unknown argument skipped.')
        endif
        goto 10
      endif

c ... Check that minimal information was entered
      if (schnam .eq. ' ') then
        call message(' ')
        call message('   No schedule file name was given.')
        call message(' ')
        STOP 'Execution error.'
      endif

      if (newbin .eq. ' ') then
        if (ext) then
          newbin = oldbin
        else
          call message(' ')
          call message('   No binary output file name was given.')
          call message(' ')
          STOP 'Execution error.'
        endif
      endif

c ... Open binary file to write to
      if (ext .and. newbin .eq. oldbin) then
        open(unit=1,file=newbin,form='UNFORMATTED',status='OLD')
      else
        open(unit=1,file=newbin,form='UNFORMATTED',status='NEW')
      endif

c ... Open cdi file if it has been specified
      if (cdifile .ne. ' ') then
        open(unit=18,file=cdifile,status='OLD')
      endif

c ... Open the schedule file and read the header lines
      open(unit=15,file=schnam,status='OLD')

c ... Allow comments at the top of a schedule file for documentation
c ... purposes.  All of the comment lines have a # character at the
c ... start and are stored at the top of the schedule file and are 
c ... ignored here.  There is no blank line permitted between the last
c ... comment line and the line of the schedule file containing the
c ... start year information.  cak - 12/30/02
      read(15,100) string
100   format(a80)
      do while (string(1:1) .eq. '#') 
        read(15,100) string
      end do

c      read(15,*) strtyr
      read(string,*) strtyr

      read(15,*) tend
      tend = tend + 1

      read(15,*) sitnam

      read(15,*) labtyp
      read(15,*) labyr

      read(15,*) mctemp
      micosm = 0
      if (mctemp .ge. 0) then
        call message(' ')
        call message('Microcosms are not implemented in this version')
        call message('of Daily Century')
        micosm = 1
        STOP
      endif

      read(15,*) co2sys
      if (co2sys .gt. 0) then
        read(15,*) co2tm(1), co2tm(2)
      endif

c ... New header line in schedule file to handle pH shift, cak - 08/02/02
c ... Change pH shift implementation to use scalar values, cak - 10/17/05
      read(15,*) phsys
      if (phsys .gt. 0) then
        read(15,*) phtm
      endif

c ... New header lines in schedule file to handle soil temperature warming
c ... experiments, cak - 07/02/03
      read(15,*) stsys
      if (stsys .gt. 0) then
        read(15,*) ststart
        read(15,*) stamt
      endif

c ... New header lines in schedule file to handle the N input scalars,
c ... cak - 04/06/04
      read(15,*) Ninput
      if (Ninput .gt. 0) then
        read(15,*) Nstart
      endif

c ... New header lines in schedule file to handle the OMAD input scalars,
c ... cak - 04/06/04
      read(15,*) OMADinput
      if (OMADinput .gt. 0) then
        read(15,*) OMADstart
      endif

c ... New header lines in schedule file to handle the weather input scalars,
c ... cak - 10/18/05
      read(15,*) wthinput
      if (wthinput .gt. 0) then
        read(15,*) wthstart
      endif

      read(15,*) decsys
      if (decsys .eq. SAVSYS) then
        decsys = FORSYS
      endif
      read(15,40) initcp
40    format(a5)
      if (initcp .eq. 'Initi') then
        initcp = ' '
      endif
      read(15,40) initre
      if (initre .eq. 'Initi') then
        initre = ' '
      endif

      read(15,*)
      read(15,*)

c ... Read starting values from fixed parameter file
      call fixin

c ... Read starting values from site-specific file
      open(unit=7,file=sitnam,status='OLD',err=1000)
      call sitein(ext)

c ... Moved the read calls for the initial tree and crop to inside the extend
c ... if statement.  This is done to prevent a rather subtle bug that occurs
c ... when the initial crop/tree do not match the final values in the 
c ... original schedule.  In that case, the derived output values, (crpval ...)
c ... do not match the current crop values.
c ... The crop/tree reads must occur before the calciv call on a normal run.
c ... 7/20/95  K. Killian

c ... Determine initial values
      if (ext) then
        if (oldbin .ne. newbin) then
          open(unit=3,file=oldbin,form='UNFORMATTED',status='OLD')
          call extend(3,.TRUE.)
          close(unit=3)
        else
          call extend(1,.FALSE.)
        endif
      endif

c ... Add new routine to initialize common block variables to other than
c ... default values as necessary
      call initialize(ext)

c ... Initialize the Daily Soil Water Model
      call initsw(sitlat, swcinit, usexdrvrs, numlyrs, texture,
     &            sradadj, tminslope, tminintercept,
     &            maxphoto, bioabsorp)

c ... Initialize soil properties based on structure of Daily Soil Water Model
      call setlyrs(adep,nlayer,numlyrs, sand, silt, clay, bulkd, ph, 
     &             awilt, afiel, swflag)

      ammonium = 0.0
      do 45 ii = 1, SWMAXLYR
        nitrate(ii) = 0.0
45    continue

c ... Zero out ALL layers below nlayer. 
      do 19 iel = 1, 3
        do 18 ii = nlayer+1, MAXLYR
          minerl(ii, iel) = 0.0
18      continue
19    continue
 
      do 110 ii=1,nlayer
        if (minerl(ii,N) .lt. 0.05) then
          minerl(ii,N) = 0.1
        endif
        call update_npool(ii, minerl(ii,N), frac_nh4_fert,
     &                    frac_no3_fert, ammonium, nitrate, subname)
110   continue

c ... Obtain the initial values for the crop or forest system
c ... Initialize the fine root pools based on the information read from
c ... the site file and the crop/tree parameterization, cak - 05/24/2007
      cursys = 0
      if (initcp .ne. ' ') then
        call cropin(initcp)
        cursys = CRPSYS
        bglcism(UNLABL) = bglcisj(UNLABL) * mrtfrac
        bglcism(LABELD) = bglcisj(LABELD) * mrtfrac
        bglcisj(UNLABL) = bglcisj(UNLABL) * (1.0 - mrtfrac)
        bglcisj(LABELD) = bglcisj(LABELD) * (1.0 - mrtfrac)
        do 120 ii = 1, MAXIEL
          bglivem(ii) = bglivej(ii) * mrtfrac
          bglivej(ii) = bglivej(ii) * (1.0 - mrtfrac)
120     continue
      endif
      if (initre .ne. ' ') then
        call treein(initre)
        cursys = cursys + FORSYS
        frtcism(UNLABL) = frtcisj(UNLABL) * wmrtfrac
        frtcism(LABELD) = frtcisj(LABELD) * wmrtfrac
        frtcisj(UNLABL) = frtcisj(UNLABL) * (1.0 - wmrtfrac)
        frtcisj(LABELD) = frtcisj(LABELD) * (1.0 - wmrtfrac)
        do 130 ii = 1, MAXIEL
          frootem(ii) = frootej(ii) * wmrtfrac
          frootej(ii) = frootej(ii) * (1.0 - wmrtfrac)
130     continue
      endif
      if (.not. ext) then
        call calciv
      endif

c ... Sum up isotopes
      call sumcar

c ... Do preliminary initializations and calculations
      call prelim

c ... Initialize initial soil moisture (asmos) based on structure of 
c ... the Daily Soil Water Model 
      call setasmos(asmos, nlayer, swcinit, numlyrs, avh2o, rwcf)

c ... Read the first block of events
      call readblk

      call message(' ')
      call message('   Model is running...')

      return
1000    call message(' Fatal error: unknown site file :'//sitnam)
        stop ' Abnormal Termination'
      end


      integer function getlen(name)

      implicit none
      character*(*) name
      integer jj

C -----------------------------------------------------------------------------
C     this subroutine left justifies the file name and determines the length
C
C Variables
C      Input
C   name    character (*)  the input and processed file name
C
C  Modified by K. Killian 8/11/94
C              included the left justification on a subroutine coded by Laura
C
C -----------------------------------------------------------------------------
 
15    getlen = index(name,' ')-1

      if (getlen .eq. -1) then
        getlen = len(name)
      else if (getlen .eq. 0) then
        do 20 jj= 1,len(name)
          if (name(jj:jj) .ne. ' ') then
            name = name(jj:)
            goto 15
          endif
20      continue
        getlen = 0
      endif

      return
      end
