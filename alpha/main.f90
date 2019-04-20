!============================================!

module setup

implicit real(a-h,o-z)

include 'Pcon.h'
include 'Grid.h'

!include 'FR.h'

integer::inlay(nlev) 

!! BG field variables

real(4)::ubg(nlev) !!! zonal  wind
real(4)::vbg(nlev) !!! merid. wind
real(4)::tbg(nlev) !!! temperature
real(4)::pbg(nlev) !!! pressure
real(4)::rhobg(nlev) !!! density

!! diagnosed variables

real(4)::thbg(nlev)!!! potential temperature 
real(4)::bn2(nlev) !!! buoyancy freq.
real(4)::alp(nlev) !!! 1/(2H) -- H: density scale height
real(4)::alp0 !!! constant alpha

real(4)::ubg_bnd(nlayer,2) !!! zonal  wind
real(4)::vbg_bnd(nlayer,2) !!! merid. wind
real(4)::bn2_bnd(nlayer,2) !!! buryancy freq.
real(4)::tbg_bnd(nlayer,2) !!! temprature


real(4)::zbnds(nlayer+1)

!! Corioli parameter

real(4),parameter::axlat_0=0.0 !!! fixed latitude
real(4),parameter::f0=2.0*eomg*sin(drad*axlat_0)

!! topography

real(4)::h(nx,ny)

real(4)::uvdir !!! surface wind directon (deg)

!! disturbance field variables

real(4)::dispx(nx,ny,nlev) !!! displacement x
real(4)::dispy(nx,ny,nlev) !!! displacement y 
real(4)::dispz(nx,ny,nlev) !!! displacement z
real(4)::u(nx,ny,nlev)     !!! zonal wind
real(4)::v(nx,ny,nlev)     !!! merid. wind
real(4)::w(nx,ny,nlev)     !!! vertical wind
real(4)::t(nx,ny,nlev)     !!! templerature
real(4)::p(nx,ny,nlev)     !!! pressure
real(4)::th(nx,ny,nlev)    !!! potential templerature


character*10::cbgmode, chmode

namelist /init_nml/ cbgmode, chmode
namelist /param_nml/ vnu

!!namelist /FR_nml/

end module setup

!============================================!
program main
use setup

open (11,file='calc.conf')
 read(11,nml=init_nml)
 read(11,nml=param_nml)
close(11)

call read_BG

call read_init

call FRmodel_3d_general

call output

call QuickLook_yz
call QuickLook_xy

stop
end program main
!============================================!
