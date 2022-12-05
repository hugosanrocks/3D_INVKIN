      Program Km2Deg

! Computes lat and long for a specified distance in km N and E of a
! reference point

! If want coords corresponding to a distance D and azimuth AZ from a reference
! point, use this program with N = D*cos(AZ) and E = D*sin(AZ)

! Sample input (col. 1 starts after the initial "*")

!Name of output file:
!  km2deg.out
!Coordinates (title, ref. coord. (+east), n, e dist to point (km)):
!  JAB
! 34.31149 -118.49675  Ref for JAB
! 0.00250  0.00963    dist to JAB
!  JGB 
! 34.31295 -118.49890  Ref for JGB
! 0.00395  0.00588    dist to JGB
! stop


! Dates: 05/27/98 - Written by D.M. Boore
!        06/01/98 - Changed to double precision
!        06/03/98 - Obtain information from a control file
!        09/16/00 - Made east long positive, added simple computation
!                   for coordinates
!        04/30 02 - Suggest name of control file
!        02/26/04 - Cleanup for use in LF95
!        02/14/09 - Change input to km2deg_f to single
!                   precision, and convert to double precision
!                   in the subprogram
!        11/18/15 - Some modernization of code

      character buf*70, f_in*30, f_out*30, buf4*4
!      double precision vn, ve, alat_ref, along_ref, vlat, vlong

      logical f_in_exist

      f_in_exist = .false.
      do while (.not. f_in_exist)
        f_in = ' '
        write(*, '(a)') 
     :    ' Enter name of control file (cr = KM2DEG.CTL;'//
     :    ' CTL-BRK to quit): '
        read(*, '(a)') f_in
        if (f_in(1:4) == '    ') f_in = 'km2deg.ctl'
        call trim_c(f_in,nc_f_in)
        inquire(file=f_in(1:nc_f_in), exist=f_in_exist)
        if (.not. f_in_exist) then
          write(*,'(a)') ' ******* FILE DOES NOT EXIST ******* '
        end if
      end do
      call get_lun(nu_in)
      open(nu_in, file=f_in(1:nc_f_in), status = 'unknown')
      write(*, '(1x,a)') f_in(1:nc_f_in)

! Obtain info from control file:

      read(nu_in,*)                            ! skip comment line
      f_out = ' '
      read(nu_in, '(a)') f_out
      call trim_c(f_out,nc_f_out)
      call get_lun(nu_out)
      open(nu_out, file=f_out(1:nc_f_out), status = 'unknown')
      
      read(nu_in,*)                            ! skip comment line

      loop_over_input: DO

        buf = ' '
        read(nu_in, '(a)') buf                   ! read title or stop
        call trim_c(buf,nc_buf)
        buf4 = ' '
        buf4 = buf(1:4)
        call upstr(buf4)
        if(buf4(1:4) == 'STOP') then
          close(nu_in)
          close(nu_out)
          stop
        end if
      
        read(nu_in,*) alat_ref, along_ref        ! reference coordinates
        read(nu_in,*) vn, ve                     ! dist. to point

        call km2deg_f( vn, ve, alat_ref, along_ref, 
     :               vlat, vlong )

        !RICHTER TABLES
        !write(nu_out,'(2(1x,f10.5))') vlat, vlong

        pi = 4.0 * atan(1.0)
        d2r = pi/180.0

        vlat = alat_ref + vn/(6371.0 * d2r)
        vlong = along_ref + ve/(6371.0 * d2r * 
     :          cos(0.5 * (alat_ref + vlat) * d2r))

        !SPHERIC PROJECTION
        write(nu_out,'(2(1x,f10.5))') vlat, vlong

      END DO loop_over_input

      end

!      include '\forprogs\km2deg_f.for'
!      include '\forprogs\locate_d.for'
!      include '\forprogs\upstr.for'
!      include '\forprogs\trim_c.for'
!      include '\forprogs\get_lun.for'




!-------------------- BEGIN KM2DEG_F ----------------------------
      subroutine km2deg_f( vn_in, ve_in, alat_ref_in, along_ref_in, 
     :               vlat_out, vlong_out )
        
! convert km north and east from a reference point into lat, long

! assumes positive latitude between 0 and 70 degrees
! assumes east longitude is positive
! assumes angles in degrees

! WARNING: NEEDS DOUBLE PRECISION VERSION OF LOCATE (ATTACHED HERE)

! Dates:  10/01/95 - written by D. Boore
!         05/27/98 - Name changed to km2deg_f
!         06/01/98 - Changed to double precision
!         02/14/09 - Changed input to single precision
!        04/28/15 - Replaced comment characters * or C with ! (The Fortran 95 standard)
               
      double precision alat_tbl(71), b_tbl(71), adcoslat_tbl(71)
      double precision vn, ve, alat_ref, along_ref, vlat, vlong
      Data alat_tbl /
     : 0.000000, 1.000000, 2.000000, 3.000000, 4.000000, 5.000000,
     : 6.000000, 7.000000, 8.000000, 9.000000,10.000000,11.000000,
     :12.000000,13.000000,14.000000,15.000000,16.000000,17.000000,
     :18.000000,19.000000,20.000000,21.000000,22.000000,23.000000,
     :24.000000,25.000000,26.000000,27.000000,28.000000,29.000000,
     :30.000000,31.000000,32.000000,33.000000,34.000000,35.000000,
     :36.000000,37.000000,38.000000,39.000000,40.000000,41.000000,
     :42.000000,43.000000,44.000000,45.000000,46.000000,47.000000,
     :48.000000,49.000000,50.000000,51.000000,52.000000,53.000000,
     :54.000000,55.000000,56.000000,57.000000,58.000000,59.000000,
     :60.000000,61.000000,62.000000,63.000000,64.000000,65.000000,
     :66.000000,67.000000,68.000000,69.000000,70.000000
     :/
      Data b_tbl /
     : 1.842808, 1.842813, 1.842830, 1.842858, 1.842898, 1.842950,
     : 1.843011, 1.843085, 1.843170, 1.843265, 1.843372, 1.843488,
     : 1.843617, 1.843755, 1.843903, 1.844062, 1.844230, 1.844408,
     : 1.844595, 1.844792, 1.844998, 1.845213, 1.845437, 1.845668,
     : 1.845907, 1.846153, 1.846408, 1.846670, 1.846938, 1.847213,
     : 1.847495, 1.847781, 1.848073, 1.848372, 1.848673, 1.848980,
     : 1.849290, 1.849605, 1.849992, 1.850242, 1.850565, 1.850890,
     : 1.851217, 1.851543, 1.851873, 1.852202, 1.852531, 1.852860,
     : 1.853188, 1.853515, 1.853842, 1.854165, 1.854487, 1.854805,
     : 1.855122, 1.855433, 1.855742, 1.856045, 1.856345, 1.856640,
     : 1.856928, 1.857212, 1.857490, 1.857762, 1.858025, 1.858283,
     : 1.858533, 1.858775, 1.859008, 1.859235, 1.859452
     :/
      Data adcoslat_tbl /
     : 1.855365, 1.855369, 1.855374, 1.855383, 1.855396, 1.855414,
     : 1.855434, 1.855458, 1.855487, 1.855520, 1.855555, 1.855595,
     : 1.855638, 1.855683, 1.855733, 1.855786, 1.855842, 1.855902,
     : 1.855966, 1.856031, 1.856100, 1.856173, 1.856248, 1.856325,
     : 1.856404, 1.856488, 1.856573, 1.856661, 1.856750, 1.856843,
     : 1.856937, 1.857033, 1.857132, 1.857231, 1.857331, 1.857435,
     : 1.857538, 1.857643, 1.857750, 1.857858, 1.857964, 1.858074,
     : 1.858184, 1.858294, 1.858403, 1.858512, 1.858623, 1.858734,
     : 1.858842, 1.858951, 1.859061, 1.859170, 1.859276, 1.859384,
     : 1.859488, 1.859592, 1.859695, 1.859798, 1.859896, 1.859995,
     : 1.860094, 1.860187, 1.860279, 1.860369, 1.860459, 1.860544,
     : 1.860627, 1.860709, 1.860787, 1.860861, 1.860934
     :/

      pi = 4.0*atan(1.0)
      d2r = pi/ 180.0

      vn = dble(vn_in)
      ve = dble(ve_in)
      alat_ref =  dble(alat_ref_in)
      along_ref =  dble(along_ref_in) 

! interpolate to find proper arc distance:

      call locate_d( alat_tbl, 71, alat_ref, j)
      b = b_tbl(j) + (alat_ref-alat_tbl(j))*
     :  (b_tbl(j+1)-b_tbl(j))/
     :  (alat_tbl(j+1)-alat_tbl(j))

      adcoslat = adcoslat_tbl(j) + (alat_ref-alat_tbl(j))*
     :  (adcoslat_tbl(j+1)-adcoslat_tbl(j))/
     :  (alat_tbl(j+1)-alat_tbl(j))

      a = adcoslat * cos(d2r*alat_ref)

      dlambda = +ve/a ! version with minus used if assume west long is +
!      dlambda = -ve/a ! minus; positve ve corresponds to decrease in long
      dphi    =  vn/b

! convert from minutes of arc to degrees:
      dlambda = dlambda / 60.0
      dphi    = dphi    / 60.0

      vlat  = alat_ref  + dphi
      vlong = along_ref + dlambda

! Consider using the simpler sphere approximation:
!      vlat = alat_ref + vn/(6371.0 * d2r)
!      vlong = along_ref + ve/(6371.0 * d2r * 
!     :        cos(0.5 * (alat_ref + vlat) * d2r))

      vlat_out = sngl(vlat)
      vlong_out = sngl(vlong)
      
      return
      end
!-------------------- END KM2DEG_F ----------------------------





! --------------------- BEGIN UPSTR ----------------------------------
      Subroutine UPSTR ( text )
! Converts character string in TEXT to uppercase
! Dates: 03/12/96 - Written by Larry Baker
!        04/28/15 - Replaced comment characters * or C with ! (The Fortran 95 standard)

!
      Implicit   None
!
      Character  text*(*)
!
      Integer    j
      Character  ch
!
      Do 1000 j = 1,LEN(text)
         ch = text(j:j)
         If ( LGE(ch,'a') .and. LLE(ch,'z') ) Then
            text(j:j) = CHAR ( ICHAR(ch) - ICHAR('a') + ICHAR('A') )
         End If
 1000    Continue
!
      Return
      End
! --------------------- END UPSTR ----------------------------------




! --------------------------- BEGIN TRIM_C -----------------------
      subroutine trim_c(cstr, nchar)

! strips leading and trailing blanks from cstr, returning the
! result in cstr, which is now nchar characters in length

! Strip off tabs also.

! Here is a sample use in constructing a column header, filled out with 
! periods:

!* Read idtag:
!        idtag = ' '
!        read(nu_in, '(1x,a)') idtag
!        call trim_c(idtag, nc_id)
!* Set up the column headings:
!        colhead = ' '
!        colhead = idtag(1:nc_id)//'......' ! nc_id + 6 > length of colhead

! Dates: 12/23/97 - written by D. Boore
!        12/08/00 - pad with trailing blanks.  Otherwise some comparisons
!                   of the trimmed character string with other strings
!                   can be in error because the original characters are left
!                   behind after shifting.  For example, here is a string
!                   before and after shifting, using the old version:
!                      col:12345
!                           MTWH  before
!                          MTWHH  after (but nc = 4).
!        03/21/01 - Check for a zero length input string
!        11/09/01 - Change check for zero length string to check for all blanks
!        10/19/09 - Strip off tabs
!        04/28/15 - Replaced comment characters * or C with ! (The Fortran 95 standard)

      character cstr*(*)

      if(cstr .eq. ' ') then
        nchar = 0
        return
      end if

      nend = len(cstr)

! Replace tabs with blanks:

      do i = 1, nend
        if(ichar(cstr(i:i)) .eq. 9) then
           cstr(i:i) = ' '
        end if
      end do



!      if(nend .eq. 0) then
!        nchar = 0
!        return
!      end if

      do i = nend, 1, -1
        if (cstr(i:i) .ne. ' ') then
           nchar2 = i
           goto 10
        end if
      end do

10    continue

      do j = 1, nchar2
        if (cstr(j:j) .ne. ' ') then
          nchar1 = j
          goto 20
        end if
      end do

20    continue
   
      nchar = nchar2 - nchar1 + 1
      cstr(1:nchar) = cstr(nchar1: nchar2)
      if (nchar .lt. nend) then
        do i = nchar+1, nend
          cstr(i:i) = ' '
        end do
      end if

      return
      end
! --------------------------- END TRIM_C -----------------------




! --------------------------- BEGIN GET_LUN ----------------
      subroutine get_lun(lun)

! Finds a logical unit number not in use; returns
! -1 if it cannot find one.

! Dates -- 05/19/98 - Written by D. Boore, following
!                     Larry Baker's suggestion
!        04/28/15 - Replaced comment characters * or C with ! (The Fortran 95 standard)

      logical isopen
      do i = 99,10,-1
        inquire (unit=i, opened=isopen)
        if(.not.isopen) then
          lun = i
          return
        end if
      end do
      lun = -1

      return
      end
! --------------------------- END GET_LUN ----------------
     



!* --------------------- BEGIN LOCATE -----------------
      SUBROUTINE locate(xx,n,x,j)
      
! Comments added by D. Boore on 26feb2010:
!  finds j such that xx(j) < x <= xx(j+1)
!  EXCEPT if x = xx(1), then j = 1 (logically it would be 0 from
!  the above relation, but this is the same returned value of j
!  for a point out of range).
!  Also, if x = xx(n), j = n-1, which is OK
!  Note that j = 0 or j = n indicates that x is out of range.
!
! See the program test_locate.for to test this routine.

! Dates: 04/28/15 - Replaced comment characters * or C with ! (The Fortran 95 standard)

      INTEGER j,n
      REAL x,xx(n)
      INTEGER jl,jm,ju
      jl=0
      ju=n+1
10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      goto 10
      endif
      if(x.eq.xx(1))then
        j=1
      else if(x.eq.xx(n))then
        j=n-1
      else
        j=jl
      endif
      return
      END
!* --------------------- END LOCATE -----------------

*-------------------- BEGIN LOCATE_D ----------------------------
      SUBROUTINE locate_d(xx,n,x,j)   ! Double precision version of locate
      INTEGER j,n
      double precision x,xx(n)
      INTEGER jl,jm,ju
      jl=0
      ju=n+1
10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      goto 10
      endif
      if(x.eq.xx(1))then
        j=1
      else if(x.eq.xx(n))then
        j=n-1
      else
        j=jl
      endif
      return
      END
* --------------------- END LOCATE_D -----------------

