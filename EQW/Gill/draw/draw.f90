!==================================================!
module setup
implicit real(a-h,o-z)

include 'Pcon.h'
include 'Grid.h'

real(4)::heat(nt,nx,ny,nlev)
real(4)::vb(nt,nx,ny,nlev)
real(4)::vu(nt,nx,ny,nlev)
real(4)::vv(nt,nx,ny,nlev)
real(4)::vw(nt,nx,ny,nlev)
real(4)::vp(nt,nx,ny,nlev)
real(4)::vpv(nt,nx,ny,nlev)

real(4)::plev(nlev)
real(4)::tbg(nlev)

real(4)::vz(nt,nx,ny,nlev)

 real,parameter::t00=300.0

 real,parameter::bn2=0.012**2


 character*100::cdir
 character*5  ::cmode
 namelist /draw_nml/ &
 iout,               &    
 cdir,               &    
 cmode

 character*30::cname
 namelist /output_nml/ &
 cname,      &
 itimev

end module setup
!==================================================!
program main
use setup

open(11,file='output.conf')           ; read(11,nml=draw_nml)   ; close(11)
open(11,file=trim(cdir)//'calc.conf') ; read(11,nml=output_nml) ; close(11)

call load

call prep

if (trim(cmode).eq.'latp')  call draw_latp
if (trim(cmode).eq.'lonp')  call draw_lonp
if (trim(cmode).eq.'map')   call draw_map

stop
end program main

!==================================================!

subroutine draw_map
use setup

real(4),allocatable::vc(:,:),vc2(:,:),vs(:,:),vs2(:,:),vh(:,:)
character*30::title1, title2
character*5::varc,vars,varc2,vars2
!-----
namelist /map_nml/ refp, range_xl, range_xr, range_yl, range_yr, varc, vars, varc2, vars2, title1, title2, iheat
!-----

open(11,file='output.conf')  ; read(11,nml=map_nml)  ; close(11)

ixl=iblkle(axx,nx,range_xl)
ixr=iblkge(axx,nx,range_xr)
iyl=iblkle(axy,ny,range_yl)
iyr=iblkge(axy,ny,range_yr)
lenx=ixr-ixl+1
leny=iyr-iyl+1

ireflev=iblkge(-plev,nlev,-refp)

allocate(vc(lenx,leny),vc2(lenx,leny),vs(lenx,leny),vs2(lenx,leny),vh(lenx,leny))



xstics=1.0
xmtics=5.0
ystics=1.0
ymtics=2.0

!ireflev=nlev


if (itimev.eq.0) ntime=1
!if (itimev.eq.1) ntime=nt
if (itimev.eq.1) ntime=58

do itime=1,ntime

select case(trim(varc))
 case('z'); vc=vz(itime,ixl:ixr,iyl:iyr,ireflev)
 case('t'); vc=vb(itime,ixl:ixr,iyl:iyr,ireflev)*t00/grav*(plev(ireflev)*1.0e2/ps0)**(rd/cp)
 case('pv');vc=vpv(itime,ixl:ixr,iyl:iyr,ireflev)
 case('w'); vc=vw(itime,ixl:ixr,iyl:iyr,ireflev)
end select

if (iheat.eq.1) vh=heat(itime,ixl:ixr,iyl:iyr,ireflev)

! output
iouts=1
ioute=2
if(iout.eq.1) ioute=1
if(iout.eq.2) iouts=2
do ioutl=iouts,ioute

 call ADCL_general(ioutl)

 call grswnd (1.0e-3*range_xl,1.0e-3*range_xr,1.0e-3*range_yl,1.0e-3*range_yr) ! set window
 call grsvpt (0.15,0.85,0.23,0.73) ! set viewport 
 call grstrn (1) ! linear or log
 call grstrf

 call uwsgxa ((1.0e-3)*axx(ixl:ixr),lenx) ! xaxis value
 call uwsgya ((1.0e-3)*axy(iyl:iyr),leny) ! yaxis value

! **** contour ****

 call udlset ('LMSG',.FALSE.) ! 'contour interval'
 call udrset ('RSIZEL',0.015) ! label
 call udiset ('INDXMJ',13) ! major line
 call udiset ('INDXMN',11) ! minor line

! call udgcla(-0.0005,0.001,0.00005)

 call udcntr (vc,lenx,lenx,leny)

 call udiclv
 call udiset ('INDXMJ',23) ! major line
 call udiset ('INDXMN',21) ! minor line
 call udcntr (vh,lenx,lenx,leny)


! **** y, z axis ****
 call uzinit
 call uziset ('INDEXT2',3)
 call uziset ('INDEXL1',5)
 call uzrset ('RSIZEL1',0.018)
 call uzrset ('RSIZEC1',0.018)
 call uzrset ('RSIZET2',0.010)
 call uzrset ('RSIZET1',0.004)
 !      call uzlset ('LABELXB',.FALSE.)
 call uzlset ('LOFFSET',.TRUE.)
 call uzrset ('XFACT',1.0e-3)
 call uzrset ('YFACT',1.0e-3)
 
 call uysfmt ('(I4)')
 call uxaxdv ('B',xstics,xmtics)
 call uxaxdv ('T',xstics,xmtics)
 call uzlset ('LABELYR',.FALSE.)
 call uzlset ('LABELYL',.TRUE.)
 call uyaxdv ('L',ystics,ymtics)
 call uyaxdv ('R',ystics,ymtics)
 call uxsttl ('B','x (10|3"km)',0.0)
 call uysttl ('L','y (10|3"km)',0.0)
 call sglset ('LCLIP',.FALSE.) ! clipping
 call sgtxzv (0.50,0.76,trim(title1),0.024,0,0,5) ! title
 call sgtxzv (0.17,0.75,trim(title2),0.016,0,-1,5) ! title
 call grcls
end do

end do


deallocate(vc,vc2,vs,vs2,vh)

return
end subroutine draw_map

!==================================================!

subroutine draw_latp
use setup

real(4),allocatable::vc(:,:),vc2(:,:),vs(:,:),vs2(:,:),vh(:,:)
character*30::title1, title2
character*5::varc,vars,varc2,vars2
real(4),parameter::vplab(5)=(/1000.,500.,300.,200.,100./)
character*4,parameter::cplab(5)=(/'1000',' 500',' 300',' 200',' 100'/)
!-----
namelist /latp_nml/ refx, range_yl, range_yr, range_levb, range_levt, varc, vars, varc2, vars2, title1, title2, iheat
!-----

open(11,file='output.conf')  ; read(11,nml=latp_nml)  ; close(11)

iyl=iblkle(axy,ny,range_yl)
iyr=iblkge(axy,ny,range_yr)

if (iyr.ge.ny+1) iyr=ny
if (iyl.le.0) iyl=1


if (range_levb.le.0.0.or.range_levt.le.0.0)then
 ilevb=1
 ilevt=nlev
 range_levb=plev(1)
 range_levt=plev(nlev)
else
 ilevb=iblkle(-plev,nlev,range_levb)
 ilevt=iblkge(-plev,nlev,range_levt)
end if

leny   = iyr-iyl+1
lenlev = ilevt-ilevb+1

ixref=iblkge(axx,nx,refx)

allocate(vc(leny,lenlev),vc2(leny,lenlev),vs(leny,lenlev),vs2(leny,lenlev),vh(leny,lenlev))



xstics=1.0
xmtics=5.0

if (itimev.eq.0) ntime=1
if (itimev.eq.1) ntime=nt

do itime=1,ntime

select case(trim(varc))
 case('z'); vc=vz(itime,ixref,iyl:iyr,ilevb:ilevt)
 case('t'); do ilev=1,lenlev; vc(:,ilev)=vb(itime,ixref,iyl:iyr,ilevb-1+ilev)*t00/grav*(plev(ilevb-1+ilev)*1.0e2/ps0)**(rd/cp) ; end do
 case('th'); do ilev=1,lenlev; vc(:,ilev)=vb(itime,ixref,iyl:iyr,ilevb-1+ilev)*t00/grav ; end do
 case('w'); vc=vw(itime,ixref,iyl:iyr,ilevb:ilevt)
end select

if (iheat.eq.1) vh=heat(itime,ixref,iyl:iyr,ilevb:ilevt)

! output
iouts=1
ioute=2
if(iout.eq.1) ioute=1
if(iout.eq.2) iouts=2
do ioutl=iouts,ioute

 call ADCL_general(ioutl)

 call grswnd (1.0e-3*range_yl,1.0e-3*range_yr,range_levb,range_levt) ! set window
 call grsvpt (0.15,0.85,0.23,0.73) ! set viewport 
 call grstrn (2) ! linear or log
 call grstrf


 call uwsgxa ((1.0e-3)*axy(iyl:iyr),leny) 
 call uwsgya (plev(ilevb:ilevt),lenlev)

! **** contour ****

 call udlset ('LMSG',.FALSE.) ! 'contour interval'
 call udrset ('RSIZEL',0.015) ! label
 call udiset ('INDXMJ',13) ! major line
 call udiset ('INDXMN',11) ! minor line

 call udcntr (vc,leny,leny,lenlev)

 call udiclv
 call udiset ('INDXMJ',23) ! major line
 call udiset ('INDXMN',21) ! minor line
 call udcntr (vh,leny,leny,lenlev)


! **** y, z axis ****
 call uzinit
 call uziset ('INDEXT2',3)
 call uziset ('INDEXL1',5)
 call uzrset ('RSIZEL1',0.018)
 call uzrset ('RSIZEC1',0.018)
 call uzrset ('RSIZET2',0.010)
 call uzrset ('RSIZET1',0.004)
 !      call uzlset ('LABELXB',.FALSE.)
 call uzlset ('LOFFSET',.TRUE.)
 call uzrset ('XFACT',1.0e-3)
 
 call uysfmt ('(I4)')
 call uxaxdv ('B',xstics,xmtics)
 call uxaxdv ('T',xstics,xmtics)
 call uzlset ('LABELYR',.FALSE.)
 call uzlset ('LABELYL',.TRUE.)
 call uyaxlb ('L',vplab,5,vplab,cplab,4,5)
 call uyaxlb ('R',vplab,5,vplab,cplab,4,5)
 call uxsttl ('B','y (10|3"km)',0.0)
 call uysttl ('L','Pressure (hPa)',0.0)
 call sglset ('LCLIP',.FALSE.) ! clipping
 call sgtxzv (0.50,0.76,trim(title1),0.024,0,0,5) ! title
 call sgtxzv (0.17,0.75,trim(title2),0.020,0,-1,5) ! title
 call grcls
end do

end do


deallocate(vc,vc2,vs,vs2,vh)

return
end subroutine draw_latp

!==================================================!


subroutine draw_lonp
use setup

real(4),allocatable::vc(:,:),vc2(:,:),vs(:,:),vs2(:,:),vh(:,:)
character*30::title1, title2
character*5::varc,vars,varc2,vars2
real(4),parameter::vplab(5)=(/1000.,500.,300.,200.,100./)
character*4,parameter::cplab(5)=(/'1000',' 500',' 300',' 200',' 100'/)
!-----
namelist /lonp_nml/ refy, range_xl, range_xr, range_levb, range_levt, varc, vars, varc2, vars2, title1, title2, iheat
!-----

open(11,file='output.conf')  ; read(11,nml=lonp_nml)  ; close(11)

ixl=iblkle(axx,nx,range_xl)
ixr=iblkge(axx,nx,range_xr)

if (ixr.ge.nx+1) ixr=nx
if (ixl.le.0) ixl=1


if (range_levb.le.0.0.or.range_levt.le.0.0)then
 ilevb=1
 ilevt=nlev
 range_levb=plev(1)
 range_levt=plev(nlev)
else
 ilevb=iblkle(-plev,nlev,range_levb)
 ilevt=iblkge(-plev,nlev,range_levt)
end if

lenx   = ixr-ixl+1
lenlev = ilevt-ilevb+1

iyref=iblkge(axy,ny,refy)

allocate(vc(lenx,lenlev),vc2(lenx,lenlev),vs(lenx,lenlev),vs2(lenx,lenlev),vh(lenx,lenlev))



xstics=1.0
xmtics=5.0

if (itimev.eq.0) ntime=1
if (itimev.eq.1) ntime=nt

do itime=1,ntime

select case(trim(varc))
 case('z'); vc=vz(itime,ixl:ixr,iyref,ilevb:ilevt)
 case('t'); do ilev=1,lenlev; vc(:,ilev)=vb(itime,ixl:ixr,iyref,ilevb-1+ilev)*t00/grav*(plev(ilevb-1+ilev)*1.0e2/ps0)**(rd/cp) ; end do
 case('th'); do ilev=1,lenlev; vc(:,ilev)=vb(itime,ixl:ixr,iyref,ilevb-1+ilev)*t00/grav ; end do
 case('w'); vc=vw(itime,ixl:ixr,iyref,ilevb:ilevt)
 case('pv'); vc=vpv(itime,ixl:ixr,iyref,ilevb:ilevt)
end select

write(*,*) maxval(vc),minval(vc)

if (iheat.eq.1) vh=heat(itime,ixl:ixr,iyref,ilevb:ilevt)

! output
iouts=1
ioute=2
if(iout.eq.1) ioute=1
if(iout.eq.2) iouts=2
do ioutl=iouts,ioute

 call ADCL_general(ioutl)

 call grswnd (1.0e-3*range_xl,1.0e-3*range_xr,range_levb,range_levt) ! set window
 call grsvpt (0.15,0.85,0.23,0.73) ! set viewport 
 call grstrn (2) ! linear or log
 call grstrf


 call uwsgxa ((1.0e-3)*axx(ixl:ixr),lenx) 
 call uwsgya (plev(ilevb:ilevt),lenlev)

! **** contour ****

 call udlset ('LMSG',.FALSE.) ! 'contour interval'
 call udrset ('RSIZEL',0.015) ! label
 call udiset ('INDXMJ',13) ! major line
 call udiset ('INDXMN',11) ! minor line

 call udcntr (vc,lenx,lenx,lenlev)

 call udiclv
 call udiset ('INDXMJ',23) ! major line
 call udiset ('INDXMN',21) ! minor line
 call udcntr (vh,lenx,lenx,lenlev)


! **** y, z axis ****
 call uzinit
 call uziset ('INDEXT2',3)
 call uziset ('INDEXL1',5)
 call uzrset ('RSIZEL1',0.018)
 call uzrset ('RSIZEC1',0.018)
 call uzrset ('RSIZET2',0.010)
 call uzrset ('RSIZET1',0.004)
 !      call uzlset ('LABELXB',.FALSE.)
 call uzlset ('LOFFSET',.TRUE.)
 call uzrset ('XFACT',1.0e-3)
 
 call uysfmt ('(I4)')
 call uxaxdv ('B',xstics,xmtics)
 call uxaxdv ('T',xstics,xmtics)
 call uzlset ('LABELYR',.FALSE.)
 call uzlset ('LABELYL',.TRUE.)
 call uyaxlb ('L',vplab,5,vplab,cplab,4,5)
 call uyaxlb ('R',vplab,5,vplab,cplab,4,5)
 call uxsttl ('B','x (10|3"km)',0.0)
 call uysttl ('L','Pressure (hPa)',0.0)
 call sglset ('LCLIP',.FALSE.) ! clipping
 call sgtxzv (0.50,0.76,trim(title1),0.024,0,0,5) ! title
 call sgtxzv (0.17,0.75,trim(title2),0.020,0,-1,5) ! title
 call grcls
end do

end do


deallocate(vc,vc2,vs,vs2,vh)

return
end subroutine draw_lonp

!==================================================!
subroutine load
use setup

write(*,*) 'loading...'

if (itimev.eq.1)then
 open(21,file=trim(cdir)//'heat',form='unformatted') ; read(21)   heat ; close(21) 
 open(21,file=trim(cdir)//'u',   form='unformatted') ; read(21)     vu ; close(21) 
 open(21,file=trim(cdir)//'v',   form='unformatted') ; read(21)     vv ; close(21) 
 open(21,file=trim(cdir)//'w',   form='unformatted') ; read(21)     vw ; close(21) 
 open(21,file=trim(cdir)//'p',   form='unformatted') ; read(21)     vp ; close(21) 
 open(21,file=trim(cdir)//'b',   form='unformatted') ; read(21)     vb ; close(21) 
 open(21,file=trim(cdir)//'pv',  form='unformatted') ; read(21)    vpv ; close(21) 
else
 open(21,file=trim(cdir)//'heat',form='unformatted') ; read(21) heat(1,:,:,:) ; close(21) 
 open(21,file=trim(cdir)//'u',   form='unformatted') ; read(21)   vu(1,:,:,:) ; close(21) 
 open(21,file=trim(cdir)//'v',   form='unformatted') ; read(21)   vv(1,:,:,:) ; close(21) 
 open(21,file=trim(cdir)//'w',   form='unformatted') ; read(21)   vw(1,:,:,:) ; close(21) 
 open(21,file=trim(cdir)//'p',   form='unformatted') ; read(21)   vp(1,:,:,:) ; close(21) 
 open(21,file=trim(cdir)//'b',   form='unformatted') ; read(21)   vb(1,:,:,:) ; close(21) 
 open(21,file=trim(cdir)//'pv',  form='unformatted'); read(21)  vpv(1,:,:,:) ; close(21) 
end if

vpv = vpv*1.0e6 !!! PVU

!write(*,*) maxval(vw(1,:,:,:))
!stop

return
end subroutine load
!==================================================!
subroutine prep
use setup

real(4)::thmid(nlev-1),zmid(nlev-1)

zmid(1:nlev-1)  = 0.5 * (zlev(1:nlev-1)+zlev(2:nlev))
thmid(1:nlev-1) = t00 + t00/grav * bn2 * (zmid(1:nlev-1)-zbot) !!! theta


plev(1) = 1.0 !!! (p/ps0) ** kp


do ilev=2,nlev
 plev(ilev) = plev(ilev-1) - grav/cp/thmid(ilev-1)  * (zlev(ilev)-zlev(ilev-1))
! plev(ilev) = plev(ilev-1) - grav/cp/t00 * (zlev(ilev)-zlev(ilev-1))
end do 

plev = plev ** (cp/rd) * ps0 * 1.0e-2 !!!! hPa

tbg(1:nlev) = t00 + t00/grav * bn2 * (zlev(1:nlev)-zbot) !!! theta
tbg(1:nlev) = tbg(1:nlev) * (plev(1:nlev)*1.0e2/ps0)**(rd/cp)



do ilev=1,nlev
 vz(:,:,:,ilev) = vp(:,:,:,ilev) * rd * tbg(ilev) / grav / plev(ilev) / 1.0e2
end do

return
end subroutine prep
!==================================================!

