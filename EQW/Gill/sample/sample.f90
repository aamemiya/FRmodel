!======================================!
module setup
implicit real(a-h,o-z)
include 'Pcon.h'
include 'Grid.h'

 real(4)::heat(ny,nlev)

 real(4)::var_yz(ny,nlev)
 real(4)::var_y(ny)
 real(4)::var_z(nlev)

 real,parameter::t00=250.0
 real,parameter::H=t00*rd/grav

 real,parameter::bn2=0.01**2
 real,parameter::beta=1.0e-11

 integer,parameter::iout=1

end module setup
!======================================!
program sample
use setup
use func

write(*,*) 'scale height = ',H,'(m)'


do in=0,2
 vm=pi*real(in)/ztop
 ce=sqrt(bn2/(vm**2+(0.5/H)**2))
 do ilev=1,nlev
  var_z(ilev)= exp(0.5*zlev(ilev)/H) * cos(vm*zlev(ilev))
 end do
do il=-1,2
 write(*,*) 'mode n=',in,', l=',il
 if (il.eq.-1)then 
  write(*,*) 'Kelvin wave'
  else
  write(*,*) 'Rossby wave'
 end if
 gamma=sqrt(2.0*beta/ce)
! write(*,*) 1.0/gamma,beta,ce
! stop
 do iy=1,ny
  var_y(iy)= para_cyl(gamma*axy(iy),il)
 end do
 
 !!! KIKA KUKA
 
 var_z=var_z / sqrt(sum(var_z**2)/real(nlev))
 var_y=var_y / sqrt(sum(var_y**2)/real(ny))
 
 do ilev=1,nlev
 do iy=1,ny
  var_yz(iy,ilev)=var_z(ilev)*var_y(iy)
 end do
 end do

 call quicklook
end do
end do


end program sample
!======================================!
subroutine quicklook
use setup

character*20::title1,title2

real(4)::u_snap(ny,nlev)
real(4)::v_snap(ny,nlev)
real(4)::w_snap(ny,nlev)
real(4)::h_snap(ny,nlev)
real(4)::b_snap(ny,nlev)
real(4)::z_snap(ny,nlev)

ystics=500.0
ymtics=1000.0



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
      call grswnd (1.0e-3*yl,1.0e-3*yr,1.0e-3*zbot,1.0e-3*ztop) ! set window
      call grsvpt (0.15,0.85,0.13,0.73) ! set viewport !!! 1:1
      call grstrn (1) ! linear or log
      call grstrf


  zasp= (ztop-zbot) / (xr-xl) 

! **** contour ****

      call udlset ('LMSG',.FALSE.) ! 'contour interval'
      call udrset ('RSIZEL',0.015) ! label
      call uwsgxa ((1.0e-3)*axy,ny) ! yaxis value
      call uwsgya ((1.0e-3)*zlev,nlev) ! zaxis value (km)
      call udcntr (var_yz,ny,ny,nlev)
! select case(trim(cvarc))
! case('u');  call udcntr (u_snap,nx,nx,nlev)
! case('w');  call udcntr (w_snap,nx,nx,nlev)
! end select

! *** vector ***

!select case (trim(cvarv))
! case ('uw'); call ugvect (u_snap,nx,w_snap/zasp,nx,nx,nlev)
!end select

! *** lines on xy plane ***

!  call uuslnt(1) ! line type
!  call uuslni(5) ! line index (color*10 + index*1)
!  do ilev=1,nlev
!   select case (trim(cvarl))
!   case ('xz');call uulin (nx,0.001*(axx(:)+dispx_snap(:,ilev)),0.001*(zlev(ile!v)+dispz_snap(:,ilev)))
!   case ('z'); call uulin (nx,0.001*(axx(:)),0.001*(zlev(ilev)+dispz_snap(:,ile!v)))
!   end select
!  end do

! **** y, z axis ****
      call uzinit
      call uziset ('INDEXT2',3)
      call uziset ('INDEXL1',5)
      call uzrset ('RSIZEL1',0.018)
      call uzrset ('RSIZEC1',0.018)
      call uzrset ('RSIZET2',0.010)
      call uzrset ('RSIZET1',0.004)
!      call uzlset ('LABELXB',.FALSE.)
!      call uzlset ('LOFFSET',.TRUE.)
 !     call uzrset ('XFACT',0.001)
      call uzlset ('LOFFSET',.FALSE.)
      call uxaxdv ('B',ystics,ymtics)
      call uxaxdv ('T',ystics,ymtics)
      call uxsttl ('B','y(km)',0.0)
      call uysfmt ('(I2)')
      call uzlset ('LABELYR',.FALSE.)
      call uzlset ('LABELYL',.TRUE.)
      call uzlset ('LOFFSET',.FALSE.)
      call uyaxdv ('L',5.0,10.0)
      call uyaxdv ('R',5.0,10.0)
      call uziset ('IROTCYL',1)
      call uysttl ('L','Height(km)',0.0)

      call sgtxzv (0.50,0.76,trim(title1),0.024,0,0,5) ! title
      call sgtxzv (0.17,0.76,trim(title2),0.020,0,-1,5) ! title
      call grcls
end do

return
end subroutine quicklook
!======================================!
