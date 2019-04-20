
!============================================!

subroutine quicklook
use setup

integer::iwork(nx,nlev)

character*3::cvarc,cvarl,cvars, cvarv
character*20::title1,title2
namelist /quicklook_nml/ iout, cvarc, cvarl, cvars, cvarv, title1, title2


open(11,file='calc.conf')
 read(11,nml=quicklook_nml)
close(11)


 write(*,*) maxval(dispz),minval(dispz)

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
      call grswnd (1.0e-3*xl,1.0e-3*xr,1.0e-3*zbot,1.0e-3*ztop) ! set window
      call grsvpt (0.15,0.75,0.13,0.63) ! set viewport !!! 1:1
      call grstrn (1) ! linear or log
      call grstrf


  zasp= (ztop-zbot) / (xr-xl) 

! **** contour ****

      call udlset ('LMSG',.FALSE.) ! 'contour interval'
      call udrset ('RSIZEL',0.015) ! label
      call uwsgxa ((1.0e-3)*axx,nx) ! yaxis value
      call uwsgya ((1.0e-3)*zlev,nlev) ! zaxis value (km)
 select case(trim(cvarc))
 case('u');  call udcntr (u,nx,nx,nlev)
 case('w');  call udcntr (w,nx,nx,nlev)
 end select

! *** vector ***

select case (trim(cvarv))
 case ('uw'); call ugvect (u,nx,w/zasp,nx,nx,nlev)
end select

! *** lines on xy plane ***

  call uuslnt(1) ! line type
  call uuslni(1) ! line index (color*10 + index*1)
  do ilev=1,nlev,2
   select case (trim(cvarl))
   case ('xz');call uulin (nx,0.001*(axx(:)+dispx(:,ilev)),0.001*(zlev(ilev)+dispz(:,ilev)))
   case ('z'); call uulin (nx,0.001*(axx(:)),0.001*(zlev(ilev)+dispz(:,ilev)))
   end select
  end do

! *** topo ***

  call uuslnt(1) ! line type
  call uuslni(3) ! line index (color*10 + index*1)
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


! **** background U and N **** 

      umin=10.*real(int(0.1*minval(ubg)))-15.0
      umax=10.*real(int(0.1*minval(ubg)))+15.0
      bn2min=minval(bn2) - 0.1 * abs(minval(bn2))
      bn2max=maxval(bn2) + 0.1 * abs(maxval(bn2))
      
      fac = ( umax - umin ) / ( bn2max - bn2min )


      write(*,*) bn2max,bn2min
      write(*,*) umax,umin,fac

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

      call uzinit
      call uziset ('INDEXT2',3)
      call uziset ('INDEXL1',5)
      call uzrset ('RSIZEL1',0.014)
      call uzrset ('RSIZEC1',0.014)
      call uzrset ('RSIZET2',0.010)
      call uzrset ('RSIZET1',0.004)
      call uxaxdv ('B',5.0,10.0)
      call uxsttl ('B','U(m/s)',0.0)

      call uzlset ('LABELXT',.true.)
      call uzlset ('LOFFSET',.true.)
      call uzrset ('XOFFSET',bn2min-umin/fac/1.0e-4)
      call uzrset ('XFACT',1/fac/1.0e-4)
      call uxaxdv ('T',0.5,1.)
      call uxsttl ('T','N2(10|-4"s|-2")',0.0)
      call uysfmt ('(I2)')
      call uzlset ('LABELYR',.FALSE.)
      call uzlset ('LABELYL',.FALSE.)


      call sgtxzv (0.45,0.65,trim(title1),0.020,0,0,5) ! title
      call sgtxzv (0.17,0.65,trim(title2),0.016,0,-1,5) ! title
      call grcls
end do

return
end subroutine quicklook

!============================================!
