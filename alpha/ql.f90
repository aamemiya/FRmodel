
!============================================!

subroutine quicklook_yz
use setup

integer::iwork(nx_data,nlev)

character*3::cvarc,cvarl,cvars, cvarv
character*20::title1,title2
namelist /quicklook_yz_nml/ iout, cvarc, cvarl, cvars, cvarv, title1, title2, cvmin,cvmax,cdv


ixs=nx/2-nx_data/2+1
ixe=nx/2+nx_data/2
iys=ny/2-ny_data/2+1
iye=ny/2+ny_data/2


iysmp=ny/2 !!! TORI AEZU

open(11,file='calc.conf')
 read(11,nml=quicklook_yz_nml)
close(11)


 write(*,*) 'max,min z',maxval(dispz),minval(dispz)
!do ilev=1,nlev
! write(*,*) 'max,min z',maxval(dispz(:,:,ilev)),minval(dispz(:,:,ilev)),ilev
!end do
!stop

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
      call swlset ('LSEP',.TRUE.) ! psfile numbering
      call slmgn (0.0,0.0,0.0,0.0) ! margin 
      call grfrm 
      call grswnd (-1.0e-3*0.5*vxlen,1.0e-3*0.5*vxlen,1.0e-3*zbot,1.0e-3*ztop) ! set window
      call grsvpt (0.15,0.75,0.13,0.63) ! set viewport !!! 1:1
      call grstrn (1) ! linear or log
      call grstrf


  zasp= (ztop-zbot) / (xr-xl) 

  call uwsgxa ((1.0e-3)*axx(ixs:ixe),nx_data) ! yaxis value
  call uwsgya ((1.0e-3)*zlev,nlev) ! zaxis value (km)

! *** Shade ***
 call ueitlv
 if ((cvmax-cvmin)/20.0.eq.cdv)then
 do ip=1,9
  call uestlv(real(ip)*cdv,real(ip+1)*cdv,55999+ip*4000)
  call uestlv(-real(ip+1)*cdv,-real(ip)*cdv,55999-ip*4000)
 end do
 huge=1.0e10
 call uestlv (cvmax,huge,95999)
 call uestlv (-huge,cvmin,15999)

 select case(trim(cvars))
 case('u');  call uetone (u(ixs:ixe,iysmp,:),nx_data,nx_data,nlev)
 case('w');  call uetone (w(ixs:ixe,iysmp,:),nx_data,nx_data,nlev)
 end select
 end if


! **** contour ****

      call udlset ('LMSG',.FALSE.) ! 'contour interval'
      call udrset ('RSIZEL',0.015) ! label
      call udgcla (cvmin,cvmax,cdv)
      call uddclv (0.0)
 select case(trim(cvarc))
 case('u');  call udcntr (u(ixs:ixe,iysmp,:),nx_data,nx_data,nlev)
 case('w');  call udcntr (w(ixs:ixe,iysmp,:),nx_data,nx_data,nlev)
 end select

! *** vector ***

select case (trim(cvarv))
 case ('uw'); call ugvect (u(ixs:ixe,iysmp,:),nx_data,w(ixs:ixe,iysmp,:)/zasp,nx_data,nx_data,nlev)
end select

! *** lines on xy plane ***

  call uuslnt(1) ! line type
  call uuslni(1) ! line index (color*10 + index*1)
  do ilev=1,nlev,2
   select case (trim(cvarl))
   case ('xz');call uulin (nx_data,0.001*(axx(ixs:ixe)+dispx(ixs:ixe,iysmp,ilev)),0.001*(zlev(ilev)+dispz(ixs:ixe,iysmp,ilev)))
   case ('z'); call uulin (nx_data,0.001*(axx(ixs:ixe)),0.001*(zlev(ilev)+dispz(ixs:ixe,iysmp,ilev)))
   end select
  end do

! *** topo ***

  call uuslnt(1) ! line type
  call uuslni(3) ! line index (color*10 + index*1)
  call uulin (nx_data,0.001*axx(ixs:ixe),0.001*h(ixs:ixe,iysmp))

! **** x ,y axis ****
      call uzinit
      call uziset ('INDEXT2',3)
      call uziset ('INDEXL1',5)
      call uzrset ('RSIZEL1',0.018)
      call uzrset ('RSIZEC1',0.018)
      call uzrset ('RSIZET2',0.010)
      call uzrset ('RSIZET1',0.004)
      call uxsfmt ('(I4)')
      call uzlset ('LOFFSET',.false.)
      call uzlset ('LABELXT',.FALSE.)
      call uzlset ('LABELXB',.TRUE.)
      call uxaxdv ('B',50.0,100.0)
      call uxaxdv ('T',50.0,100.0)
      call uxsttl ('B','x(km)',0.0)
      call uysfmt ('(I2)')
      call uzlset ('LABELYR',.FALSE.)
      call uzlset ('LABELYL',.TRUE.)
      call uyaxdv ('L',5.0,10.0)
      call uyaxdv ('R',5.0,10.0)
      call uziset ('IROTCYL',1)
      call uysttl ('L','Height(km)',0.0)


! **** background U and N **** 

      umin=10.*real(int(0.1*minval(ubg)))-15.0
      umax=10.*real(int(0.1*maxval(ubg)))+15.0
      bn2min=minval(bn2) - 0.1 * abs(minval(bn2))
      bn2max=maxval(bn2) + 0.1 * abs(maxval(bn2))
      
      fac = ( umax - umin ) / ( bn2max - bn2min )
!      bn2min=
!      bn2max=


      write(*,*) maxval(ubg),minval(ubg)
      write(*,*) umax,umin
      write(*,*) bn2max,bn2min
      write(*,*) fac
! stop

      call grsvpt (0.80,0.92,0.13,0.63) ! set new window
      call grswnd (umin,umax,0.001*zbot,0.001*ztop) ! set viewport
      call grstrn (1) ! linear or log
      call grstrf

      call uuslnt (1)
      call uuslni (2)

!      do ilayer=1,nlayer
!       call uulin (2,(/ ubg(ilayer), ubg(ilayer) /), 0.001*(/ zbnds(ilayer), zbnds(ilayer+1) /) )
!      end do

!       call uuslnt (3)
!      do ilayer=1,nlayer
!       call uulin (2,(/ umin+fac*(bn2(ilayer)-bn2min), umin+fac*(bn2(ilayer)-bn2min) /), 0.001*(/ zbnds(ilayer), zbnds(ilayer+1) /) )
!      end do

      call uulin  (nlev,ubg,0.001*zlev)
       call uuslnt (3)
      call uulin  (nlev,umin+fac*(bn2-bn2min),0.001*zlev)

do ilev=1,nlev
 write(*,*) 1.0e-3*zlev(ilev),bn2(ilev),umin+fac*(bn2(ilev)-bn2min)
end do


      call uzinit
      call uziset ('INDEXT2',3)
      call uziset ('INDEXL1',5)
      call uzrset ('RSIZEL1',0.014)
      call uzrset ('RSIZEC1',0.014)
      call uzrset ('RSIZET2',0.010)
      call uzrset ('RSIZET1',0.004)
      call uxsfmt ('(I2)')
      call uzlset ('LOFFSET',.false.)
      if ((umax-umin).gt.80.0)then
      call uxaxdv ('B',20.0,40.0)
      elseif ((umax-umin).gt.40.0)then
      call uxaxdv ('B',10.0,20.0)      
      else
      call uxaxdv ('B',5.0,10.0)
      end if
      call uxsttl ('B','U(m/s)',0.0)

      call uzlset ('LABELXT',.true.)
      call uzlset ('LOFFSET',.true.)
      call uzrset ('XOFFSET',(bn2min-umin/fac)/1.0e-4)
!      call uzrset ('XOFFSET',-umin/fac/1.0e-4)
      call uzrset ('XFACT',1.0/fac/1.0e-4)
      call uxaxdv ('T',0.5,1.)
      call uxsttl ('T','N|2"(10|-4"s|-2")',0.0)
!      call uysfmt ('(I2)')
!      call uzlset ('LABELYR',.FALSE.)
!      call uzlset ('LABELYL',.FALSE.)


      call sgtxzv (0.45,0.65,trim(title1),0.020,0,0,5) ! title
      call sgtxzv (0.17,0.65,trim(title2),0.016,0,-1,5) ! title
      call grcls
end do

return
end subroutine quicklook_yz

!============================================!

subroutine quicklook_xy
use setup

integer::iwork(nx_data,ny_data)

character*3::cvarc,cvarl,cvars, cvarv
character*20::title1,title2
namelist /quicklook_xy_nml/ reflev,iout, cvarc,  cvars, cvarv, title1, title2, cvmin,cvmax,cdv

open(11,file='calc.conf')
 read(11,nml=quicklook_xy_nml)
close(11)


ixs=nx/2-nx_data/2+1
ixe=nx/2+nx_data/2
iys=ny/2-ny_data/2+1
iye=ny/2+ny_data/2

izsmp=iblkge(zlev,nlev,reflev) !!! TORI AEZU

open(11,file='calc.conf')
 read(11,nml=quicklook_xy_nml)
close(11)


 write(*,*) 'max,min h',maxval(h),minval(h)
 write(*,*) 'max,min z',maxval(dispz(ixs:ixe,iys:iye,izsmp)),minval(dispz(ixs:ixe,iys:iye,izsmp))

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
      call swlset ('LSEP',.TRUE.) ! psfile numbering
      call slmgn (0.0,0.0,0.0,0.0) ! margin 
      call grfrm 
      call grswnd (-1.0e-3*0.5*vxlen,1.0e-3*0.5*vxlen,-1.0e-3*0.5*vylen,1.0e-3*0.5*vylen) ! set window
      call grsvpt (0.15,0.75,0.13,0.63) ! set viewport !!! 1:1
      call grstrn (1) ! linear or log
      call grstrf


  zasp= (ztop-zbot) / vxlen
  yasp= vylen / vxlen

! **** contour ****

      call udlset ('LMSG',.FALSE.) ! 'contour interval'
      call udrset ('RSIZEL',0.015) ! label
      call uwsgxa ((1.0e-3)*axx(ixs:ixe),nx_data) ! yaxis value
      call uwsgya ((1.0e-3)*axy(iys:iye),ny_data) ! zaxis value (km)
      call udgcla (cvmin,cvmax,cdv)
      call uddclv (0.0)
 select case(trim(cvarc))
  case('u');  call udcntz (u(ixs:ixe,iys:iye,izsmp),nx_data,nx_data,ny_data,iwork,nx_data*ny_data)
  case('v');  call udcntz (v(ixs:ixe,iys:iye,izsmp),nx_data,nx_data,ny_data,iwork,nx_data*ny_data)
  case('w');  call udcntz (w(ixs:ixe,iys:iye,izsmp),nx_data,nx_data,ny_data,iwork,nx_data*ny_data)
  case('x');  call udcntz (dispx(ixs:ixe,iys:iye,izsmp),nx_data,nx_data,ny_data,iwork,nx_data*ny_data)
  case('y');  call udcntz (dispy(ixs:ixe,iys:iye,izsmp),nx_data,nx_data,ny_data,iwork,nx_data*ny_data)
  case('z');  call udcntz (dispz(ixs:ixe,iys:iye,izsmp),nx_data,nx_data,ny_data,iwork,nx_data*ny_data)
 end select

! *** vector ***

select case (trim(cvarv))
 case ('uv'); call ugvect (u(ixs:ixe,iys:iye,izsmp),nx_data,v(ixs:ixe,iys:iye,izsmp)/yasp,nx_data,nx_data,ny_data)
 case ('xy'); call ugvect (dispx(ixs:ixe,iys:iye,izsmp),nx_data,dispy(ixs:ixe,iys:iye,izsmp)/yasp,nx_data,nx_data,ny_data)
end select


! *** topo ***

  call udiclv
  call udlset ('LMSG',.FALSE.) ! 'contour interval'
  call udrset ('RSIZEL',0.015) ! label
  call udcntz(h(ixs:ixe,iys:iye),nx_data,nx_data,ny_data,iwork,nx_data*ny_data)

! **** x ,y axis ****
      call uzinit
      call uziset ('INDEXT2',3)
      call uziset ('INDEXL1',5)
      call uzrset ('RSIZEL1',0.018)
      call uzrset ('RSIZEC1',0.018)
      call uzrset ('RSIZET2',0.010)
      call uzrset ('RSIZET1',0.004)
      call uxsfmt ('(I4)')
      call uysfmt ('(I4)')
      call uzlset ('LOFFSET',.false.)
      call uzlset ('LABELXT',.FALSE.)
      call uzlset ('LABELXB',.TRUE.)
      call uxaxdv ('B',50.0,100.0)
      call uxaxdv ('T',50.0,100.0)
      call uxsttl ('B','x (km)',0.0)
      call uzlset ('LABELYR',.FALSE.)
      call uzlset ('LABELYL',.TRUE.)
      call uyaxdv ('L',50.0,100.0)
      call uyaxdv ('R',50.0,100.0)
      call uziset ('IROTCYL',1)
      call uysttl ('L','y (km)',0.0)

      call sgtxzv (0.45,0.65,trim(title1),0.020,0,0,5) ! title
      call sgtxzv (0.17,0.65,trim(title2),0.016,0,-1,5) ! title
      call grcls
end do

return
end subroutine quicklook_xy

!============================================!
