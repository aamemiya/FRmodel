
integer,parameter::nlon=240,nlat=121,nlev=37
integer,parameter::ntmax=124

character*5::chour(4)=(/'00:00','06:00','12:00','18:00'/)

real(4),parameter::axlon(nlon)=(/(1.5*real(ilon-1),ilon=1,nlon)/)
  real(4),parameter::axlatSN(nlat)=(/(-90.0+real(ilat-1)*1.5,ilat=1,nlat)/)
  real(4),parameter::axlatNS(nlat)=(/(90.0-real(ilat-1)*1.5,ilat=1,nlat)/)

real(4),parameter::plev(nlev)=(/&
       1000., 975., 950., 925., 900., 875., 850., 825., 800., 775., &
        750., 700., 650., 600., 550., 500., 450., 400., 350., 300., &
        250., 225., 200., 175., 150., 125., 100.,  70.,  50.,  30., &
         20.,  10.,   7.,   5.,   3.,   2.,   1./)

real(4),parameter::plev_inv(nlev)=plev(nlev:1:-1)
