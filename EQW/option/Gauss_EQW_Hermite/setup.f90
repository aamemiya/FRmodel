!==============================================================!
module param

integer,parameter::ny=64
integer,parameter::nh=16

real(4),parameter::ylen=2.0e6
real(4),parameter::yl=-8.0e6,yr=8.0e6

integer::nhe
real(4)::ylene   !!! effecitve
real(4)::yle,yre !!! effecitve
real(4)::axy_gauss(nh)
real(4)::wty_gauss(nh)
real(4)::axy(ny)

end module param
!====================================!
module output
use param

real(4)::f(ny)
real(4)::f_new(ny)
real(4)::f_new_eqw(ny)
real(4)::f_new_comp(ny,nh)
real(4)::f_new_comp_eqw(ny,nh)

real(4)::fh(nh)

character*20::title1,title2
character*20::xlabel,ylabel
real(4)::vmin,vmax
real(4)::astics,amtics
real(4)::bstics,bmtics
integer::iout

end module output
!==============================================================!
