
!============================================!

subroutine quicklook
use setup

integer::iwork(nx,nlev)

character*3::cvarc,cvarl,cvars, cvarv
character*20::title1,title2
namelist /quicklook_nml/ iout, cvarc, cvarl, cvars, cvarv, title1, title2



real(4)::u_snap(nx,nlev)
real(4)::w_snap(nx,nlev)
real(4)::dispx_snap(nx,nlev)
real(4)::dispz_snap(nx,nlev)
real(4)::h_snap(nx)
real(4)::h_diag_snap(nx)



open(11,file='calc.conf')
 read(11,nml=quicklook_nml)
close(11)


iy=ny/2

do it=1,nt,2

u_snap=u(:,it,iy,:)
w_snap=w(:,it,iy,:)
dispx_snap=dispx(:,it,iy,:)
dispz_snap=dispz(:,it,iy,:)
h_snap=h(:,it,iy)
h_diag_snap=h_diag(:,it,iy)

 write(*,*) 'test'
 write(*,*) h(nx/2,it,ny/2-5:ny/2-3)
 write(*,*) h_diag(nx/2,it,ny/2-5:ny/2-3)
 write(*,*) dispz(nx/2,it,ny/2-5:ny/2-3,10)



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
      call grswnd (1.0e-3*xl,1.0e-3*xr,1.0e-3*zbot,1.0e-3*ztop) ! set window
      call grsvpt (0.15,0.85,0.13,0.73) ! set viewport !!! 1:1
      call grstrn (1) ! linear or log
      call grstrf


  zasp= (ztop-zbot) / (xr-xl) 

! **** contour ****

      call udlset ('LMSG',.FALSE.) ! 'contour interval'
      call udrset ('RSIZEL',0.015) ! label
      call uwsgxa ((1.0e-3)*axx,nx) ! yaxis value
      call uwsgya ((1.0e-3)*zlev,nlev) ! zaxis value (km)
 select case(trim(cvarc))
 case('u');  call udcntr (u_snap,nx,nx,nlev)
 case('w');  call udcntr (w_snap,nx,nx,nlev)
 end select

! *** vector ***

select case (trim(cvarv))
 case ('uw'); call ugvect (u_snap,nx,w_snap/zasp,nx,nx,nlev)
end select

! *** lines on xy plane ***

  call uuslnt(1) ! line type
  call uuslni(5) ! line index (color*10 + index*1)
  do ilev=1,nlev
   select case (trim(cvarl))
   case ('xz');call uulin (nx,0.001*(axx(:)+dispx_snap(:,ilev)),0.001*(zlev(ilev)+dispz_snap(:,ilev)))
   case ('z'); call uulin (nx,0.001*(axx(:)),0.001*(zlev(ilev)+dispz_snap(:,ilev)))
   end select
  end do

! *** topo ***

  call uuslnt(1) ! line type
  call uuslni(1) ! line index (color*10 + index*1)
  call uulin (nx,0.001*axx(:),0.001*h_snap(:))
  call uuslnt(1) ! line type
  call uuslni(22) ! line index (color*10 + index*1)
  call uulin (nx,0.001*axx(:),0.001*h_diag_snap(:))

! **** x ,y axis ****
      call uziset ('INDEXT2',3)
      call uziset ('INDEXL1',5)
      call uzrset ('RSIZEL1',0.018)
      call uzrset ('RSIZEC1',0.018)
      call uzrset ('RSIZET2',0.010)
      call uzrset ('RSIZET1',0.004)
!      call uzlset ('LABELXB',.FALSE.)
      call uxaxdv ('B',5000.0,10000.0)
      call uxaxdv ('T',5000.0,10000.0)
      call uxsttl ('B','x(km)',0.0)
      call uysfmt ('(I2)')
      call uzlset ('LABELYR',.FALSE.)
!     call uzlset ('LABELYL',.FALSE.)
      call uyaxdv ('L',5.0,10.0)
      call uyaxdv ('R',5.0,10.0)
      call uziset ('IROTCYL',1)
      call uysttl ('L','Height(km)',0.0)

      call sgtxzv (0.50,0.76,trim(title1),0.024,0,0,5) ! title
      call sgtxzv (0.17,0.76,trim(title2),0.020,0,-1,5) ! title
      call grcls
end do

end do

return
end subroutine quicklook

!============================================!
