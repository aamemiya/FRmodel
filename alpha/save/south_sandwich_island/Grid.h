!!! Grid configuration !!!


integer,parameter::nx=128,ny=128,nlev=100
integer,parameter::nx2=nx/2,ny2=ny/2

integer,parameter::nx_data=64,ny_data=64
real(4),parameter::vxlen=6.0e5,vylen=6.0e5

integer,parameter::nlayer=2

real(4),parameter::xr=0.5*vxlen*real(nx)/real(nx_data),xl=-xr
real(4),parameter::yr=0.5*vylen*real(ny)/real(ny_data),yl=-yr

real(4),parameter::zbot=0.0,ztop=5.0e4


real(4),parameter::axx_data(nx_data)   &
  = (/(-0.5*vxlen+(real(ix-1))*vxlen/real(nx_data),ix=1,nx_data)/)
real(4),parameter::axy_data(ny_data)   &
  = (/(-0.5*vylen+(real(iy-1))*vylen/real(ny_data),iy=1,ny_data)/)


real(4),parameter::axx(nx)   &
  = (/(xl+(real(ix-1))*(xr-xl)/real(nx),ix=1,nx)/)
real(4),parameter::axy(ny)   &
  = (/(yl+(real(iy-1))*(yr-yl)/real(ny),iy=1,ny)/)
real(4),parameter::zlev(nlev)    &
  = (/(zbot+(real(ilev-1))*(ztop-zbot)/real(nlev),ilev=1,nlev)/)

real(4),parameter::axk(nx)   &
  = (/ (/(2.0*pi/(xr-xl)*ik,ik=0,nx2)/) , (/(-2.0*pi/(xr-xl)*(nx2-ik),ik=1,nx2-1)/)  /)
real(4),parameter::axl(ny)   &
  = (/ (/(2.0*pi/(yr-yl)*il,il=0,ny2)/) , (/(-2.0*pi/(yr-yl)*(ny2-il),il=1,ny2-1)/)  /)


