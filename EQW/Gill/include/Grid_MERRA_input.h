
integer,parameter::mlon=288,mlat=144,mlev=42
integer,parameter::mtmax=248

character*5::chour(8)=(/'00:00','03:00','06:00','09:00','12:00','15:00','18:00','21:00'/)

  real(4),parameter::axmlon(mlon)=(/(1.25*real(ilon-1),ilon=1,mlon)/)
  real(4),parameter::axmlatSN(mlat)=(/(-90.0+(real(ilat)-0.5)*1.25,ilat=1,mlat)/)
  real(4),parameter::axmlatNS(mlat)=(/(90.0-(real(ilat)-0.5)*1.25,ilat=1,mlat)/)

real(4),parameter::pmlev(mlev)=(/&
       1000., 975., 950., 925., 900., 875., 850., 825., 800., 775., &
        750., 725., 700., 650., 600., 550., 500., 450., 400., 350., &
        300., 250., 200., 150., 100.,  70.,  50.,  40.,  30.,  20., &
	 10.,   7.,   5.,   4.,   3.,   2.,   1.,  0.7,  0.5,  0.4, &
	 0.3,  0.1   &
       /)

real(4),parameter::pmlev_inv(mlev)=pmlev(mlev:1:-1)
