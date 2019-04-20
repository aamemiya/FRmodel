
!============================================!

subroutine quicklook
use setup

character*3::cvar
character*20::title1,title2
namelist /quicklook_nml/ iout, cvar, title1, title2


open(11,file='calc.conf')
 read(11,nml=quicklook_nml)
close(11)


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
      call grsvpt (0.15,0.85,0.13,0.73) ! set viewport
      call grstrn (1) ! linear or log
      call grstrf


! **** contour ****

      call udlset ('LMSG',.FALSE.) ! 'contour interval'
      call udrset ('RSIZEL',0.015) ! label
      call uwsgxa (axx,nx) ! yaxis value
      call uwsgya (1.0e-3*zlev,nlev) ! zaxis value (km)

! *** lines on xy plane ***

 write(*,*) maxval(dispz),minval(dispz)

  call uuslnt(1) ! line type
  call uuslni(5) ! line index (color*10 + index*1)
  do ilev=1,nlev
   call uulin (nx,0.001*axx(:),0.001*(zlev(ilev)+dispz(:,ilev)))
  end do

! *** topo ***

  call uuslnt(1) ! line type
  call uuslni(1) ! line index (color*10 + index*1)
  call uulin (nx,0.001*axx(:),0.001*h(:))

! **** x ,y axis ****
      call uziset ('INDEXT2',3)
      call uziset ('INDEXL1',5)
      call uzrset ('RSIZEL1',0.018)
      call uzrset ('RSIZEC1',0.018)
      call uzrset ('RSIZET2',0.010)
      call uzrset ('RSIZET1',0.004)
!      call uzlset ('LABELXB',.FALSE.)
      call uxaxdv ('B',10.0,50.0)
      call uxaxdv ('T',10.0,50.0)
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

return
end subroutine quicklook

!============================================!
