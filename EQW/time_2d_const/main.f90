!============================================!

module setup

implicit real(a-h,o-z)

include 'Pcon.h'
include 'Grid.h'

!include 'FR.h'

!! BG field variables

real(4)::ubg(nlev) !!! zonal  wind
real(4)::tbg(nlev) !!! temperature 

!! diagnosed variables

real(4)::thbg(nlev)!!! potential temperature 
real(4)::bn2(nlev) !!! buoyancy freq.
real(4)::alp(nlev) !!! 1/(2H) -- H: density scale height
real(4)::alp0 !!! constant alpha

!! topography

real(4)::h(nx,nt,ny)
real(4)::h_diag(nx,nt,ny)

!! disturbance field variables

real(4)::dispx(nx,nt,ny,nlev) !!! displacement x
real(4)::dispz(nx,nt,ny,nlev) !!! displacement z
real(4)::u(nx,nt,ny,nlev)     !!! zonal wind
real(4)::w(nx,nt,ny,nlev)     !!! vertical wind
real(4)::t(nx,nt,ny,nlev)     !!! templerature
real(4)::p(nx,nt,ny,nlev)     !!! pressure
real(4)::th(nx,nt,ny,nlev)    !!! potential templerature

character*10::cbgmode

namelist /init_nml/ cbgmode

namelist /param_nml/ vnu

!!namelist /FR_nml/

end module setup

!============================================!
program main
use setup
use spl_psp

vnu=0.0

open (11,file='calc.conf')
 read(11,nml=init_nml)
 read(11,nml=param_nml)
close(11)

call set_grid_hermite(nh,ylen,yl,yr,axy_gauss,wty_gauss)

!do ih=1,nh
! write(*,*) ih,axy_gauss(ih),wty_gauss(ih)
!end do
!stop

axy=(/( yl + (yr-yl) * (real(iy)-0.5)/real(ny) , iy=1,ny)/)

call spl_init

call read_BG

call read_init

call FRmodel_time_EQ_const

call output

call QuickLook

stop
end program main
!============================================!
