!!! Grid configuration !!!


integer,parameter::nx=64,nlev=64
integer,parameter::nx2=nx/2

real(4),parameter::xl=-1.0e5,xr=1.0e5
real(4),parameter::zbot=0.0,ztop=1.0e4

real(4),parameter::axx(nx)   &
  = (/(xl+(real(ix-1))*(xr-xl)/real(nx),ix=1,nx)/)
real(4),parameter::zlev(nlev)    &
  = (/(zbot+(real(ilev-1))*(ztop-zbot)/real(nlev),ilev=1,nlev)/)

  real(4),parameter::axk(nx)   &
  = (/ (/(2.0*pi/(xr-xl)*ik,ik=0,nx2)/) , (/(-2.0*pi/(xr-xl)*(nx2-ik),ik=1,nx2-1)/)  /)
