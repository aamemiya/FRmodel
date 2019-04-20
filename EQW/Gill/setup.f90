!======================================!
module setup
implicit real(a-h,o-z)
include 'Pcon.h'
include 'Grid.h'

!!! input
 real(4)::heat(nt,nx,ny,nlev)

!!! output
 real(4)::var_yz(ny,nlev)
 real(4)::var_yz_2(ny,nlev)
 real(4)::var_y(ny)
 real(4)::var_z(nlev)
 real(4)::var_xy(nx,ny)
 real(4)::var_xy_2(nx,ny)

 real,parameter::t00=300.0

 real,parameter::bn2=0.012**2
 real,parameter::beta=1.1e-11

 namelist /gill_nml/ &
 iout,       & 
 ifile,      & 
 igill,      &
 efold,      &
 h0_day,     &
 zloc_h,     &
 zlen_h,     &  
 xloc_h,     &
 xlen_h,     &
 yloc_h,     &
 ylen_h,     &
 xlen_nkcut, &
 period0,    &
 peri_nocut, &
 nhcut,      &
 nncut,      &
 icuth       

 character*30::cname

 namelist /output_nml/ &
 cname,      &
 itimev
 !!! common 

 real(4)::eps0, H
 integer::nkcut,nocut,nhcut,nncut

end module setup
!======================================!
