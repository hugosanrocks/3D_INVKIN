      Program Dst4List

* Compute distance and azimuths for pairs of coordinates:

* Sample control file:
*!Control file for program dst4list
*!Name of output file:
*  dst4list.out
*!Sign of west longitude:
* -1.0
*!Coordinates of reference station:
* 35.706  -121.102    San Simeon eq
*!List of coordinates for which distance from reference station
*!will be computed (stacode a string up to 30 characters long):
*  paso 35.633  -120.700 
*  pkdb 35.94523 -120.541557 
* stop

* Reserve first 5 characters for station (should use Chuck Mueller's rcc, rcf)
 
* Dates: 07/25/06 - Written by D. Boore (modified from dst4pair.for;
*                   eliminated double precision computations).
*        08/07/06 - Allow longer file names
*        11/24/07 - Trap for colocated pair in subroutine distaz
*        08/04/09 - Use rcc, rcf to read station coordinate line
*                   and remove tabs
!        06/13/10 - Allow stacode to be 30 characters
!        06/16/12 - In output, allow longer comment lines and comment after 
!                   reference station coordinates                   
!        11/18/15 - Some modernization of code

      character comment_line(50)*200

      character f_ctl*100, f_out*100, buf*200, buf4*5,
     :          stacode*30
      
      integer nc_comment_line
      
      logical f_ctl_exist

      f_ctl_exist = .false.
      do while (.not. f_ctl_exist)
        f_ctl = ' '
        write(*, '(a\)') 
     :    ' Enter name of control file ("Enter"=DST4LIST.CTL;'//
     :    ' CTL-BRK to quit): '
        read(*, '(a)') f_ctl
        if (f_ctl(1:4) == '    ') f_ctl = 'dst4list.ctl'
        call trim_c(f_ctl,nc_f_ctl)
        inquire(file=f_ctl(1:nc_f_ctl), exist=f_ctl_exist)
        if (.not. f_ctl_exist) then
          write(*,'(a)') ' ******* FILE DOES NOT EXIST ******* '
        end if
      end do

      call get_lun(nu_ctl)
      open(nu_ctl, file=f_ctl(1:nc_f_ctl), status = 'unknown')
      write(*, '(1x,a)') ' Control file: '//f_ctl(1:nc_f_ctl)

* Obtain info from control file:

      call skipcmnt(nu_ctl, comment_line, n_comment_lines)
      f_out = ' '
      read(nu_ctl, '(a)') f_out
      call trim_c(f_out,nc_f_out)
      call get_lun(nu_out)
      open(nu_out, file=f_out(1:nc_f_out), status = 'unknown')
      
      call skipcmnt(nu_ctl, comment_line, n_comment_lines)
      read(nu_ctl,*) sign_west_longitude
      write(*,'(1x,a,1x,f5.1)') 
     : ' sign_west_longitude = ', sign_west_longitude
      write(nu_out,'(1x,a,1x,f5.1)') 
     : ' sign_west_longitude = ', sign_west_longitude
          
      call skipcmnt(nu_ctl, comment_line, n_comment_lines)
      do i = 1, n_comment_lines
        call trim_c(comment_line(i), nc_comment_line)
        write(nu_out,'(a)') comment_line(i)(1:nc_comment_line)
      end do
      buf = ' '
      read(nu_ctl,'(a)') buf
      call trim_c(buf,nc_buf)
      read(buf(1:nc_buf),*) alat_ref, along_ref
      write(*,'(1x,a)') buf(1:nc_buf)
      write(nu_out,'(1x,a)') buf(1:nc_buf)
      
      write(nu_out,'(26x,a, 5x,a,      5x,a,    6x,a,    6x,a, 
     :               7x,a,   8x,a,  4x,a, 3x,a)') 
     :               'scode', 'alatr', 'alongr', 'alat', 'along', 
     :               'rdeg', 'rkm', 'az', 'baz'
            
      call skipcmnt(nu_ctl, comment_line, n_comment_lines)

      loop_over_input: DO

        buf = ' '
        read(nu_ctl, '(a)') buf
        call trim_c(buf,nc_buf)

        write(*, '(1x,a)') buf(1:nc_buf)
        
        buf4 = ' '
        buf4 = buf(1:4)
        call trim_c(buf4,nc_buf4)
        call upstr(buf4)
        if(buf4(1:4) == 'STOP') then
          write(*,'(a)') ' Results written to '//f_out(1:nc_f_out)
          close(nu_ctl)
          close(nu_out)
          stop
        end if

        stacode = ' '
        call RCC (1,buf,nc_buf,stacode,nc_stacode,ibgn,iend)
      
        call RCF (2,buf,nc_buf,alat,ibgn,iend, ISTAT)
        call RCF (3,buf,nc_buf,along,ibgn,iend, ISTAT)
        print *,' stacode, alat, along: ', 
     :          stacode(1:nc_stacode), alat, along 
 
        call distaz( sign_west_longitude, alat_ref, along_ref, 
     :                   alat, along, 
     :                   rdeg, rkm, az, baz)
        write(nu_out, '(1x,es10.3)') rkm

      END DO loop_over_input

      end

      include '\forprogs\distazdp.for'
      include '\forprogs\distaz.for'
      include '\forprogs\get_lun.for'
      include '\forprogs\upstr.for'
      include '\forprogs\trim_c.for'
      include '\forprogs\skip.for'
      include '\forprogs\skipcmnt.for'
      include '\forprogs\rc_subs.for'
