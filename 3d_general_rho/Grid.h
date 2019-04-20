!!! Grid configuration !!!


integer,parameter::nx=120,ny=120,nlev=60
integer,parameter::nx2=nx/2,ny2=ny/2

integer,parameter::nlayer=2

real(4),parameter::xl=-2.0e5,xr=2.0e5
real(4),parameter::yl=-2.0e5,yr=2.0e5
real(4),parameter::zbot=0.0,ztop=3.0e4

!!!real(4),parameter::zbnds(nlayer+1)=(/ zbot,0.5e4,1.0e4,ztop /)
real(4),parameter::zbnds(nlayer+1)=(/ zbot,1.5e4,ztop /)

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


real(4),parameter::p00=1.0e5 !!! surface pressure (dif. from ps0)
