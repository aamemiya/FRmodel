program vertical 

implicit real(a-h,o-z)

integer,parameter::nlev=100
real(4),parameter::zbot=0.0,ztop=1.0e5

real(4),parameter::zlev(nlev+1)=(/(zbot+(ztop-zbot)*real(ilev-1)/real(nlev)  ,ilev=1,nlev+1)/)

real(4),parameter::p00=1.0e5,t00=300.0
real(4),parameter::rd=287.1,cp=1004.0,grav=9.81


real(4)::bn2(nlev+1),p(nlev+1),rho(nlev+1),theta(nlev+1),t(nlev+1),tmp(nlev+1)


character*20::title1,title2


!!!bn2=(1.6e-2)**2 !!! const

bn_0=1.0e-2
bn_1=2.1e-2
bn_2=1.7e-2
bn_3=2.5e-2

ilev1=12 !!! tropopause
ilev2=50 !!! stratopause
ilev3=80 !!! mesopause

bn2(1:nlev+1) = (/ (/(bn_0**2, ilev=1,ilev1)/), &
                   (/(bn_1**2, ilev=ilev1+1,ilev2)/),&
                   (/(bn_2**2, ilev=ilev2+1,ilev3)/),  &
                   (/(bn_3**2, ilev=ilev3+1,nlev+1)/)   /)

!!! smoothing
 tmp=bn2
do ilev=3,nlev-1
 bn2(ilev)=sum(tmp(ilev-2:ilev+2))/5.0
end do

p(1)=p00
t(1)=t00
theta(1)=t00
rho(1)=p(1)/rd/t(1)


tmp(1)=log(theta(1))
do ilev=2,nlev+1
 tmp(ilev)=tmp(ilev-1)+ bn2(ilev) * (zlev(ilev)-zlev(ilev-1)) / grav
end do

theta=exp(tmp)

tmp(1)=1.0
do ilev=2,nlev+1
 tmp(ilev)=tmp(ilev-1)- grav/cp/theta(ilev) * (zlev(ilev)-zlev(ilev-1))
end do

p=p00*tmp**(cp/rd)

t=theta*(p/p00)**(rd/cp)
rho=p/rd/t

write(*,*) ' == result == '
write(*,*) 'ilev    z    p(hPa)   T(K)   th(K)  rho*1e3'

do ilev=nlev+1,1,-1
 write(*,'(I4, F6.1, F9.3, F7.1, F8.1, F10.4)') ilev, 0.001*zlev(ilev), 0.01*p(ilev), t(ilev), theta(ilev), 1000.*rho(ilev)
end do
write(*,*) 'ilev    z    p(hPa)   T(K)   th(K)  rho*1e3'



vmin=(int(0.1*minval(t))-1)*10.0
vmax=min((int(0.1*maxval(t))+2)*10.0,500.0)

iout=1

! output
iouts=1
ioute=2
if(iout.eq.1) ioute=1
if(iout.eq.2) iouts=2
do ioutl=iouts,ioute
! *** general settings ***
      call sgiset ('IFONT',1)
      call swistx ('ICLRMAP',14) ! colormap blue-white-red 
      call swcmll
      call swcset ('FNAME','figure')
      call swiset ('IWIDTH',1000) ! output window size
      call swiset ('IHEIGHT',800)     
      call gropn(ioutl) 
      call sglset ('LFULL',.TRUE.) ! using fullsize
      call swlset ('LSEP',.FALSE.) ! psfile numbering
      call slmgn (0.0,0.0,0.0,0.0) ! margin 
      call grfrm 
      call grsvpt (0.17,0.77,0.13,0.73) ! set new window
      call grswnd (vmin,vmax,zbot,ztop) ! set viewport
      call grstrn (1) ! linear or log
      call grstrf

      call uuslnt (1)
      call uuslni (2)

      call uulin  (nlev,t,zlev)
      call uzinit
      call uziset ('INDEXT2',3)
      call uziset ('INDEXL1',5)
      call uzrset ('RSIZEL1',0.014)
      call uzrset ('RSIZEC1',0.014)
      call uzrset ('RSIZET2',0.010)
      call uzrset ('RSIZET1',0.004)
      call uxaxdv ('B',10.0,50.0)
      call uxsttl ('B','T(K)',0.0)

      call uzlset ('LOFFSET',.true.)
      call uzrset ('YOFFSET',0.0)
      call uzrset ('YFACT',1.e-3)
      call uyaxdv ('L',10.0,50.0)
      call uysttl ('L','z (km)',0.0)
      call uysfmt ('(I2)')
      call uzlset ('LABELYR',.FALSE.)
      call uzlset ('LABELYL',.TRUE.)


      call sgtxzv (0.45,0.65,trim(title1),0.020,0,0,5) ! title
      call sgtxzv (0.17,0.65,trim(title2),0.016,0,-1,5) ! title
      call grcls

end do


stop
end program vertical
