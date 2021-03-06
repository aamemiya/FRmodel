!!! Grid configuration !!!

integer,parameter::nx=64,nt=64,ny=40,nlev=40
integer,parameter::nx2=nx/2,nt2=nt/2

real(4),parameter::xl=-2.0e6,xr=2.0e6
real(4),parameter::tl=-1.0e5,tr=1.0e5

real(4),parameter::yl=-2.0e6,yr=2.0e6

real(4),parameter::zbot=0.0,ztop=2.0e4

integer,parameter::nlayer=2
real(4),parameter::zbnds(nlayer+1)=(/ zbot,1.0e4,ztop /)

real(4),parameter::axx(nx)   &
  = (/(xl+(real(ix-1))*(xr-xl)/real(nx),ix=1,nx)/)
real(4),parameter::axy(ny)   &
  = (/(yl+(real(iy-1))*(yr-yl)/real(ny),iy=1,ny)/)
real(4),parameter::axt(nt)   &
  = (/(tl+(real(it-1))*(tr-tl)/real(nt),it=1,nt)/)
real(4),parameter::zlev(nlev)    &
  = (/(zbot+(real(ilev-1))*(ztop-zbot)/real(nlev),ilev=1,nlev)/)

real(4),parameter::axk(nx)   &
  = (/ (/(2.0*pi/(xr-xl)*ik,ik=0,nx2)/) , (/(-2.0*pi/(xr-xl)*(nx2-ik),ik=1,nx2-1)/)  /)
real(4),parameter::axo(nt)   &
  = - (/ (/(2.0*pi/(tr-tl)*io,io=0,nt2)/) , (/(-2.0*pi/(tr-tl)*(nt2-io),io=1,nt2-1)/)  /)

  integer,parameter::nj = (zbnds(2)-zbnds(1)) / ((ztop-zbot)/real(nlev)) -1

  integer,parameter::nh=30 !!! tori aezu
