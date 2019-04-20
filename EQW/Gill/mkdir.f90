subroutine mkdir_f90(outdir)
! Description:
!
! Author: am
!
! Host: aofd30
! Directory: /work2/am/teaching/21.Fortran/03.Makedir
!
! Revision history:
!  2011-05-10 11:42
!    Initial Version

!  use
  implicit none

  character(len=*),intent(in) :: outdir
  character(len=1000) :: comm ! unix command
  integer::is,ie,idxodr
  integer lnblnk

!  write(*,'(a)')'Program makedir starts.'
!  write(*,*)''
!
!  outdir="output/temp/"
  idxodr=lnblnk(outdir)

!
  is=1;  ie=is+len('if [ ! -d ')
  write(comm(is:ie),'(A)') 'if [ ! -d '
  is=ie+1; ie=is+idxodr
  write(comm(is:ie),'(A)') outdir(1:idxodr)
  is=ie+1
  ie=is+len(' ]; then (echo "Create directory,')
  write(comm(is:ie),'(A)')&
&' ]; then (echo "Create directory,'
  is=ie+1; ie=is+idxodr-1
  write(comm(is:ie),'(A)') outdir(1:idxodr)
  is=ie+1
  ie=is+len('"; mkdir -p ')-1
  write(comm(is:ie),'(A)')&
& '"; mkdir -p '
  is=ie+1; ie=is+idxodr-1
  write(comm(is:ie),'(A)') outdir(1:idxodr)
  is=ie+1;   ie=is+4
  write(comm(is:ie),'(A)') '); fi'

!  print '(A)',comm(1:ie)

  call system(comm(1:ie))

!  write(*,'(a)')'Done program makedir.'
!  write(*,*)
end subroutine mkdir_f90
