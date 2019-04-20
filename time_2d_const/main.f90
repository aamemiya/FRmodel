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

!! heating

real(4)::h(nx,nt)

!! disturbance field variables

real(4)::dispx(nx,nt,nlev) !!! displacement x
real(4)::dispz(nx,nt,nlev) !!! displacement z
real(4)::u(nx,nt,nlev)     !!! zonal wind
real(4)::w(nx,nt,nlev)     !!! vertical wind
real(4)::t(nx,nt,nlev)     !!! templerature
real(4)::p(nx,nt,nlev)     !!! pressure
real(4)::th(nx,nt,nlev)    !!! potential templerature

character*10::cbgmode

namelist /init_nml/ cbgmode

namelist /param_nml/ vnu

!!namelist /FR_nml/

end module setup

!============================================!
program main
use setup

vnu=0.0

open (11,file='calc.conf')
 read(11,nml=init_nml)
 read(11,nml=param_nml)
close(11)

call read_BG

call read_init

call FRmodel_time_2d_const

call output

call QuickLook

stop
end program main
!============================================!
