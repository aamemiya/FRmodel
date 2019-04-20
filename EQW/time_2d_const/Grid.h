!!! Grid configuration !!!


integer,parameter::nx=64,nt=64,nlev=30
integer,parameter::nx2=nx/2,nt2=nt/2

real(4),parameter::xl=-1.0e7,xr=1.0e7
real(4),parameter::zbot=0.0,ztop=2.0e4
real(4),parameter::tl=-1.0e6,tr=1.0e6

real(4),parameter::axx(nx)   &
  = (/(xl+(real(ix-1))*(xr-xl)/real(nx),ix=1,nx)/)
real(4),parameter::axt(nt)   &
  = (/(tl+(real(it-1))*(tr-tl)/real(nt),it=1,nt)/)
real(4),parameter::zlev(nlev)    &
  = (/(zbot+(real(ilev-1))*(ztop-zbot)/real(nlev),ilev=1,nlev)/)

real(4),parameter::axk(nx)   &
  = (/ (/(2.0*pi/(xr-xl)*ik,ik=0,nx2)/) , (/(-2.0*pi/(xr-xl)*(nx2-ik),ik=1,nx2-1)/)  /)
real(4),parameter::axo(nt)   &
  = - (/ (/(2.0*pi/(tr-tl)*io,io=0,nt2)/) , (/(-2.0*pi/(tr-tl)*(nt2-io),io=1,nt2-1)/)  /)


integer,parameter::ny=64
integer,parameter::nh=16 !!! WN

real(4),parameter::ylen=1.0e6

  real(4)::axy_gauss(nh)
  real(4)::wty_gauss(nh)
  real(4)::axy(ny)
  real(4):: yl,yr
