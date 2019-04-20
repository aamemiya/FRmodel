!============================================!

module setup

implicit real(a-h,o-z)

include 'Pcon.h'
include 'Grid.h'

!include 'FR.h'

integer::inlay(nlev) 

!! BG field variables

real(4)::ubg(nlayer) !!! zonal  wind
real(4)::bn2(nlayer) !!! buoyancy freq.

!real(4)::ubg(nlev) !!! zonal  wind
!real(4)::tbg(nlev) !!! temperature 

!! diagnosed variables

real(4)::thbg(nlev)!!! potential temperature 
!real(4)::bn2(nlev) !!! buoyancy freq.
real(4)::alp(nlev) !!! 1/(2H) -- H: density scale height
real(4)::alp0 !!! constant alpha

!! heating

real(4)::heat(nx,nt,ny,nlev)
real(4)::heat_diag(nx,nt,ny,nlev)


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


if (real(nj+1)*(ztop-zbot)/real(nlev).ne.(zbnds(2)-zbnds(1)))then
 write(*,*) 'error: zbnds(2) and zlev is not consistent'
 write(*,*) 'delta(zlev)=',real(nj+1)*(ztop-zbot)/real(nlev),'delta(zbnds)=',(zbnds(2)-zbnds(1))
end if

call read_BG

call read_init

call decomp_test

call output

call QuickLook_yz
call QuickLook_xy

stop
end program main
!============================================!
