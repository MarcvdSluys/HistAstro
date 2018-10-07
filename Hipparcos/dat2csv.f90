!> \file dat2csv.f90  Convert the selection from the Hipparcos catalog from .dat (space separated) to csv
!!                    (comma separated values) and test whether the new file can be read

!***********************************************************************************************************************************
program dat2csv
  implicit none
  integer :: totLines
  
  ! Convert the .dat file to a .csv file:
  call convert(totLines)
  
  ! Test the newly created .csv file:
  call testcsv(totLines)
  
  write(*,*)
end program dat2csv
!***********************************************************************************************************************************


!***********************************************************************************************************************************
subroutine convert(totLines)
  implicit none
  integer, intent(out) :: totLines
  integer, parameter :: nCols=9
  integer :: status,ip,op,ln, iDat
  real(8) :: dat(nCols)
  character :: inFile*(99),outFile*(99), dumStr*(999)
  
  inFile = 'combihip.dat'
  ip = 91
  open(unit=ip,form='formatted',status='old',action='read',position='rewind',file=trim(inFile),iostat=status)
  if(status.ne.0) then
     write(0,'(A)') 'Error opening file: '//trim(inFile)//' - aborting'
     stop
  end if
  
  outFile = 'combihip.csv'
  op = 92
  open(unit=op,form='formatted',status='replace',action='write',position='rewind',file=trim(outFile),iostat=status)
  if(status.ne.0) then
     write(0,'(A)') 'Error opening file: '//trim(outFile)//' - aborting'
     stop
  end if
  
  write(*,'(A)') 'Converting '//trim(inFile)//' to '//trim(outFile)//'...'
  
  ! Add an explanatory header to the output file:
  write(op,'(A6,",", A6,",", 2(A14,","), 2(A9,","), 4(A7,","))') '# Hip', 'mag', 'RA','Dec', 'pmRA','pmDec', &
       'ErRA','ErDec','ErPa','ErPd'
  
  ! Process file header:
  do ln=1,1
    read(ip,'(A)') dumStr
    write(op,'(A)') trim(dumStr)
  end do
  
  ! Process file body:
  ln = 0
  do
     ln = ln+1
     
     read(ip,*, iostat=status) iDat, dat
     if(status.lt.0) exit
     if(status.gt.0) then
        write(0,'(A,I0,A)') 'Error reading file: '//trim(inFile)//' line ',ln,' - aborting'
        stop
     end if
     
     write(op,'(I6,",", F6.2,",", 2(F14.10,","), 2(F9.2,","), 4(F7.2,","))') iDat, dat
  end do  ! ln
  
  close(ip)
  close(op)
  
  totLines = ln-1
  write(*,'(I0,A)') totLines, ' lines converted'
  
end subroutine convert
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Test .csv file

subroutine testcsv(totLines)
  implicit none
  integer, intent(in) :: totLines
  integer, parameter :: nCols=9
  integer :: status,ip,ln, iDat
  real(8) :: dat(nCols), limMag
  character :: inFile*(99),dumStr*(999)
  
  inFile = 'combihip.csv'
  write(*,'(/,A)') 'Testing csv file '//trim(inFile)//'...'
  
  ip = 93
  open(unit=ip,form='formatted',status='old',action='read',position='rewind',file=trim(inFile),iostat=status)
  if(status.ne.0) then
     write(0,'(A)') 'Error opening file: '//trim(inFile)//' - aborting'
     stop
  end if
  
  ! Read file header:
  read(ip,'(A)') dumStr
  
  ! Read and report limiting magnitude:
  read(ip,*) limMag
  write(*,'(/,A,F6.2)') 'Limiting magnitude: ', limMag
  
  ! Read file body and print first and last ten entries:
  ln = 0
  write(*,'(/,A)') 'First ten lines:'
  do
     ln = ln+1
     
     read(ip,*, iostat=status) iDat, dat  ! 'Unformatted'
     !read(ip,'(I6,1x, F6.2,1x, 2(F14.10,1x), 2(F9.2,1x), 4(F7.2,1x))', iostat=status) iDat, dat  ! 'Formatted'
     if(status.lt.0) exit
     if(status.gt.0) then
        write(0,'(A,I0,A)') 'Error reading file: '//trim(inFile)//' line ',ln,' - aborting'
        stop
     end if
     
     if(ln .le. 10 .or. ln .gt. totLines-10) write(*,'(I6, F6.2, 2(F14.10), 2(F9.2), 4(F7.2))') iDat, dat
     if(ln .eq. 11) write(*,'(/,A)') 'Last ten lines:'
  end do  ! ln
  
  close(ip)
end subroutine testcsv
!***********************************************************************************************************************************
