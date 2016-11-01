      SUBROUTINE newXS(species,nw,wavl,wav,T,DEN,j,sq,columndepth,zy,
     &    IO2)
      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)

      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc' ! photolabel

*************************************************************************************
* newXS takes in a species name (called from Photo.f while iterating over ISPEC)    *
* and standard parameters from old Xsections. Reactions.rx files are iterated over  *
* in this code, the new versions of which have updated PHOTO labels.                *
* PHOTS = ISOS, PHOTP = prefactoring PHOTT = temperature dependent PHOTO = standard *
* Photolysis reactions that can be done with standard code will be handled by       *
* XS_Standard. Anything else will be passed to old Xsections (XS). Standard code is *
* handled by reading the reaction list, interpolating the cross section file once   *
* then passing that data to XS_Standard along with the file names for quantum yield *
* data. XS_Standard interpolates QY and uses standard temperature dependence if     *
* needed in combination with interpolated XS to solve for SQ. photolabel is updated *
* the same as it was in old XSections. Eventually format reading here will need to  *
* be permanently updated, once species.dat/reactions.rx files are updated to match  *
* Titan formatting.                                                                 *
*************************************************************************************

!!! The format reading will need to change to match Titan's reaction.rx
!!! This will involve using prod4 (I think, double check with Ryan)

! TO RUN NEW: Change the call on line 93 of Initphoto to CALL newXS
! Make sure PHOTO labels are updated in Reactions.rx (I have a python script to do this)

!------- EWS - these are the variables that need to be declared here-------!
      dimension columndepth(KJ,NZ)  !testing this as an optional parameter
      INTEGER nw,j,IO2
      REAL*8 wavl(nw+1), wav(nw)
      REAL*8 T(NZ), DEN(NZ)
      REAL*8 sq(kj,nz,kw)
      REAL*8 zy
      integer n1, n2, n3 ! number of lines of data
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      INTEGER kdata, i, stat, labnum
      REAL*8, allocatable :: x1(:), x2(:), x3(:), y1(:), y2(:), y3(:)
      REAL*8 yg1(nw), yg2(nw), yg3(nw)

!     File name stuff
      character*8 species
      character*11 photolabel
      character*10 reac1, reac2, prod1, prod2, prod3
      character*5 label
      character*60 xsname, qyname
      character*25 reacs, prods

      character*60 line
      integer temps
      character*5 dumtemp, tempA, tempB, tempC

! Open reactions.rx and find current species
! This read style relies on a single species' reaction paths being together in reactions.rx
      open(12, file='PHOTOCHEM/INPUTFILES/reactions.rx', status='OLD')

! use this format for now:
667   FORMAT(A10, A10, A10, A10, A8, A5)
! new formatting (commented out for now, will update):
!667   FORMAT(A10,A10,A10,A10,A10,10X,A5)

!  Read until first photo species match (error if species not in .rx)
!  Matching first 4 characters to PHOT for photolysis
      stat = 0
      READ(12,667,IOSTAT=stat) reac1,reac2,prod1,prod2,prod3,label
      ! print *,label
      do while((trim(reac1).ne.trim(species).OR.label(1:4).ne.'PHOT')
     &  .AND.stat == 0)
        READ(12,667,IOSTAT=stat) reac1,reac2,prod1,prod2,prod3,label
      end do
      if(stat.ne.0) then
        print*, " *** Species not in reactions.rx!!!! *** "
        print*, reac1,reac2,prod1,prod2,prod3,label
        STOP
      endif

      ! If this reaction has been marked as non-standard call old xsections
      if(label(5:5).ne.'O') then
        CALL XS(species,nw,wavl,wav,T,DEN,j,sq,columndepth,zy,IO2)
      else
! Now at line of reactions.rx where first occurence of spcies is
        labnum = 0
* Determine species cross section
        xsname = 'PHOTOCHEM/DATA/XSECTIONS/'//trim(species)//'/'//
     &    trim(species)//'.XS.dat'
        open(11, file=trim(xsname), status='OLD')
  66    FORMAT(E12.6, 4X, E12.6)
  67    FORMAT(E12.6, 4X, E12.6, 4X, E12.6)
  68    FORMAT(E12.6, 4X, E12.6, 4X, E12.6, 4X, E12.6)
  ! count lines
        kdata = 0
        stat = 0
        !!!! dumtemp will get filled by "temp:" due to file style
        do while(stat == 0)
          kdata = kdata + 1
          ! on temperature line
          if(kdata == 3) then
            READ(11, '(a)', IOSTAT=stat) line
            CALL countColumns(line, temps)
            if(temps == 1) then
              READ(line, *) dumtemp, tempA
            else if(temps == 2) then
              READ(line, *) dumtemp, tempA, tempB
            else if(temps == 3) then
              READ(line, *) dumtemp, tempA, tempB, tempC
            endif
          else
            READ(11, *, IOSTAT=stat)
          endif
        end do

        close(11)
        kdata = kdata - 1 ! stat will /= 0 when it reaches EOF, so kdata was incremented once too many
        n1 = kdata - 4 ! subtract the four headers
        n2 = kdata - 4
        n3 = kdata - 4
        kdata = kdata !need +4 for inter but - 4 for headers so leave it

        allocate(x1(kdata))
        allocate(x2(kdata))
        allocate(x3(kdata))
        allocate(y1(kdata))
        allocate(y2(kdata))
        allocate(y3(kdata))

  ! read data
        open(11, file=xsname, status='OLD')
        DO i = 1, 4
          READ(11, *) ! skip header info
        ENDDO

! Get wavelength (x) and cross sections (y)
        i = 1 ! reset i to 1, but still on 5th line of reactions.rx
        stat = 0
        do while(i <= n1)
          if(temps == 2) then
            READ(11, 67, IOSTAT=stat) x1(i), y1(i), y2(i)
          else if(temps == 3) then
            READ(11, 68, IOSTAT=stat) x1(i), y1(i), y2(i), y3(i)
          else
            ! temps == 0 or 1 wil still give 1 other column
            READ(11, 66, IOSTAT=stat) x1(i), y1(i)
          endif
          i = i + 1
        enddo
        close(11)

        x2 = x1
        x3 = x1

! XS GRID WORK

        CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)
        CALL addpnt(x1,y1,kdata,n1,               zero,zero)
        CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
        CALL addpnt(x1,y1,kdata,n1,            biggest,zero)

        CALL inter2(nw+1,wavl,yg1,n1,x1,y1,ierr)   !inter2 is used for discrete points -> bins

        IF (ierr .NE. 0) THEN
          WRITE(*,*) ierr,' ***Something wrong in XS_Standard*** '
          STOP
        ENDIF

        ! need to interpolate at each temperature
        if(temps > 1) then

          CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),zero)
          CALL addpnt(x2,y2,kdata,n2,               zero,zero)
          CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),zero)
          CALL addpnt(x2,y2,kdata,n2,            biggest,zero)

          CALL inter2(nw+1,wavl,yg2,n2,x2,y2,ierr)   !inter2 is used for discrete points -> bins
          IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr,' ***Something wrong in XS_Standard*** '
            STOP
          ENDIF
        endif

        if(temps > 2) then
          CALL addpnt(x3,y3,kdata,n3,x3(1)*(1.-deltax),zero)
          CALL addpnt(x3,y3,kdata,n3,               zero,zero)
          CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),zero)
          CALL addpnt(x3,y3,kdata,n3,            biggest,zero)

          CALL inter2(nw+1,wavl,yg3,n3,x3,y3,ierr)   !inter2 is used for discrete points -> bins
          IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr,' ***Something wrong in XS_Standard*** '
            STOP
          ENDIF
        endif

! QY LOOP
        ! This loop will call subroutines on each reaction path
        ! Continues reading reactions.rx
        do while(trim(reac1).eq.trim(species).AND.
     &    label(1:4).eq.'PHOT')

          if(label(5:5).ne.'O') then
           CALL XS(species,nw,wavl,wav,T,DEN,j,sq,columndepth,zy,IO2)
          else

            labnum = labnum + 1
            reacs = trim(reac1) // " " // trim(reac2)
            prods = trim(prod1) // " " // trim(prod2) // " "
     &        // trim(prod3)

            do i=1, len(trim(reacs))
              if(reacs(i:i) == ' ') then
                reacs(i:i) = '_'
              endif
            end do

            do i=1, len(trim(prods))
              if(prods(i:i) == ' ') then
                prods(i:i) = '_'
              endif
            enddo

! qyname is the quantum yield file
            qyname = 'PHOTOCHEM/DATA/XSECTIONS/'//trim(species)//'/'//
     &        trim(reacs)//"_"//trim(prods)//".QY.dat"
! Call new xs subroutine
            CALL XS_Standard(qyname,yg1,yg2,yg3,nw,wavl,T,tempA,
     &       tempB,tempC,temps,species,j,sq,labnum)
           endif

!  move to next line in file, use label to exit loop if EOF
          stat = 0
          READ(12, 667, IOSTAT=stat) reac1,reac2,prod1,prod2,prod3,
     &        label

          if(stat.ne.0) then
            label = 'DONE ' !breaks while
          endif
        enddo !closes do while for current species in reactions.rx
        deallocate(x1)
        deallocate(x2)
        deallocate(x3)
        deallocate(y1)
        deallocate(y2)
        deallocate(y3)
      endif ! closes if(label(5) /= O)
      close(12)

      RETURN
      END SUBROUTINE

!!!!  Dont make changes yet, but keep in mind this loop will need to be changed
!!!!  so that qy is summed as you go over wavelength.

***************************************************************
* qyname = quantum yield file name                            *
* yg1,2,3 = interpolated cross sections                       *
* nw, wl, tlev, jn, sq  = same old stuff                      *
* tempA,B,C are temperature values (if specified in XS file)  *
* temps = count of temperatures (0 - 3)                       *
***************************************************************

      SUBROUTINE XS_Standard(qyname,yg1,yg2,yg3,nw,wl,tlev,tempA,
     &   tempB,tempC,temps,species,jn,sq,labnum)
      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PHOTABLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc' ! photolabel

      INTEGER nw
      REAL*8 yg1(nw), yg2(nw), yg3(nw) ! xs
      REAL*8 yq1(nw) ! interpolated qy
      REAL*8 sq(kj,nz,kw)
      REAL*8 tlev(NZ)
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      REAL*8, allocatable :: qy (:), x2 (:) ! qy is uninterpolated qy, x2 is wavelength from qy file
      INTEGER i, iw, kdata, temps, jn, labnum
      INTEGER ierr
      INTEGER stat
      INTEGER n1 ! number of lines of data
      character*8 species
      CHARACTER*11 plab, photolabel
      CHARACTER*60 qyname
      CHARACTER*5 tempA, tempB, tempC
      integer ta, tb, tc
      REAL*8 m, x, b
      REAL*8 wl(nw+1)

      ierr = 0

66    FORMAT(E12.6, 4X, E12.6)
! READ QY FILE
! count lines
      open(13, file=trim(qyname), status='OLD')
      kdata = 0
      stat = 0
      do while(stat == 0)
        READ(13, *, IOSTAT=stat)
        kdata = kdata + 1
      end do
      kdata = kdata - 1 ! quick fix to stat check EOF
      n1 = kdata - 4 ! subtract the four headers
      kdata = kdata !need +4 -4
      close(13)

      open(13, file=trim(qyname), status='OLD')
! Ignore headers
      do i = 1, 4
        READ(13, *)
      enddo
! Get QY
      allocate(qy(kdata))
      allocate(x2(kdata))

      i = 1
      stat = 0
      do while(i <= n1)
        READ(13, 66, IOSTAT=stat) x2(i), qy(i)
        i = i + 1
      enddo
      close(13)

! QY GRID WORK
! adding quantum yield points and interpolating them to yq1

      CALL addpnt(x2,qy,kdata,n1,x2(1)*(1.-deltax),zero)
      CALL addpnt(x2,qy,kdata,n1,               zero,zero)
      CALL addpnt(x2,qy,kdata,n1,x2(n1)*(1.+deltax),zero)
      CALL addpnt(x2,qy,kdata,n1,            biggest,zero)

      CALL inter2(nw+1,wl,yq1,n1,x2,qy,ierr)   !inter2 is used for discrete points -> bins
      IF (ierr .NE. 0) THEN
        print*, "FAILED IN QY"
        WRITE(*,*) ierr, ' ***Something wrong in XS_Standard*** '
        STOP
      ENDIF

! sq work
      if(temps <= 1) then
        DO iw = 1, nw
          DO i = 1, nz
            sq(jn,i,iw) = yg1(iw)*yq1(iw)
          ENDDO
        ENDDO
      else if(temps == 2) then
        read(tempA(1:len_trim(tempA)-1), *) ta ! Remove the K and put tempA into integer ta
        read(tempB(1:len_trim(tempB)-1), *) tb

        DO iw = 1, nw ! interpolated vector range
          DO i = 1, nz ! tlev range
            if(tlev(i) <= ta) then
              sq(jn,i,iw) = yg1(iw)*yq1(iw)
            else if(tlev(i) <= tb) then
              b = yg1(iw)
              m = (yg2(iw) - yg1(iw)) / (tb - ta)
              x = tlev(i) - ta
              sq(jn,i,iw) = (m*x + b)*(yq1(iw))
            else
              sq(jn,i,iw) = yg2(iw)*yq1(iw)
            endif
          ENDDO
        ENDDO
      else
        read(tempA(1:len_trim(tempA)-1), *) ta ! Remove the K and put tempA into integer ta
        read(tempB(1:len_trim(tempB)-1), *) tb
        read(tempC(1:len_trim(tempC)-1), *) tc

        DO iw = 1, nw
          DO i = 1, nz
            if(tlev(i) <= ta) then
              sq(jn,i,iw) = yg1(iw)*yq1(iw)
            else if(tlev(i) <= tb) then
              b = yg1(iw)
              m = (yg2(iw) - yg1(iw)) / (tb - ta)
              x = tlev(i) - ta
              sq(jn,i,iw) = (m*x + b)*(yq1(iw))
            else if(tlev(i) <= tc) then
              b = yg2(iw)
              m = (yg3(iw) - yg2(iw)) / (tc - tb)
              x = tlev(i) - tb
              sq(jn,i,iw) = (m*x + b)*(yq1(iw))
            else
              sq(jn,i,iw) = yg3(iw)*yq1(iw)
            endif
          enddo
        enddo
      endif

! photolabel
      write(plab, '(a, a, i1)') 'P',trim(species)//'_',labnum
      photolabel(jn) = plab
      jn = jn + 1

      DEALLOCATE(x2)
      DEALLOCATE(qy)

      RETURN
      END SUBROUTINE

! Helper function
      SUBROUTINE countColumns(line, cols)
      character*60 :: line
      integer i, n, cols

      n = len_trim(line)
      cols = 0
      i = 1

      do while(i <= n)
        ! Iterate the spaces
        do while(line(i:i) == ' ')
          i = i + 1
          if(i > n) exit
        enddo
        ! No more spaces and i still <= n
        ! we are in a column
        cols = cols + 1
        ! Iterate i until end of line (n) is reached or another space
        do
          i = i + 1
          if(i > n) exit
          if (line(i:i) == ' ') exit
        enddo
      enddo

      cols = cols - 1 ! Dont count column = "Temp:"
      RETURN
      END SUBROUTINE
