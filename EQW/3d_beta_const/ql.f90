
!============================================!

subroutine quicklook_yz
use setup

integer::iwork(nx,nlev)


real(4)::u_snap(nx,nlev)
real(4)::w_snap(nx,nlev)
real(4)::dispx_snap(nx,nlev)
real(4)::dispz_snap(nx,nlev)
real(4)::heat_snap(nx,nlev)

character*3::cvarc,cvarl,cvars, cvarv
character*20::title1,title2
logical::lsub
namelist /quicklook_nml/ iout, cvarc, cvarl, cvars, cvarv, title1, title2, lsub
namelist /range_nml/ hmin, hmax, wmin, wmax


heat_fact =  grav/300.0 / 86400.0


open(11,file='calc.conf')
 read(11,nml=quicklook_nml)
 read(11,nml=range_nml)
close(11)

iy=ny/2

!do it=1,nt,1
do it=nt/2,nt/2

write(title2,'(A,I5,A)')'time= ',int(axt(it)-axt(1)),' (s)'

heat_snap=heat(:,it,iy,:) / heat_fact
w_snap=w(:,it,iy,:)

! write(*,*) it, maxval(dispz_snap),minval(dispz_snap)
 write(*,*) 'heat',maxval(heat_snap),minval(heat_snap)
 write(*,*) 'w',maxval(w_snap),minval(w_snap)

wmin=minval(w_snap)
wmax=maxval(w_snap)

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
      call grsvpt (0.15,0.75,0.13,0.63) ! set viewport !!! 1:1
      call grstrn (1) ! linear or log
      call grstrf


  zasp= (ztop-zbot) / (xr-xl) 

! **** contour ****
      call udiclv
      call udlset ('LMSG',.FALSE.) ! 'contour interval'
      call udrset ('RSIZEL',0.015) ! label
      call udiset ('INDXMJ',13) ! major line
      call udiset ('INDXMN',11) ! minor line
      call uwsgxa ((1.0e-3)*axx,nx) ! yaxis value
      call uwsgya ((1.0e-3)*zlev,nlev) ! zaxis value (km)
 select case(trim(cvarc))
 case('u')
  call udgcla(umin,umax,-10.0)
  call uddclv(0.0)
  call udcntr (u_snap,nx,nx,nlev)
 case('w')
  call udgcla(wmin,wmax,-10.0)
  call uddclv(0.0)
  call udcntr (w_snap,nx,nx,nlev)
 end select

! *** vector ***

select case (trim(cvarv))
 case ('uw'); call ugvect (u_snap,nx,w_snap/zasp,nx,nx,nlev)
end select

! *** lines on xy plane ***

  call uuslnt(1) ! line type
  call uuslni(5) ! line index (color*10 + index*1)
  do ilev=1,nlev,2
   select case (trim(cvarl))
   case ('xz');call uulin (nx,0.001*(axx(:)+dispx_snap(:,ilev)),0.001*(zlev(ilev)+dispz_snap(:,ilev)))
   case ('z'); call uulin (nx,0.001*(axx(:)),0.001*(zlev(ilev)+dispz_snap(:,ilev)))
   end select
  end do

! *** heating ***
  call udiclv
  call udiset ('INDXMJ',23) ! major line
  call udiset ('INDXMN',21) ! minor line
  call udgcla(hmin,hmax,-10.0)
  call udcntr(heat_snap,nx,nx,nlev)

! **** x ,y axis ****
      call uzinit
      call uziset ('INDEXT2',3)
      call uziset ('INDEXL1',5)
      call uzrset ('RSIZEL1',0.018)
      call uzrset ('RSIZEC1',0.018)
      call uzrset ('RSIZET2',0.010)
      call uzrset ('RSIZET1',0.004)
      call uzlset ('LOFFSET',.FALSE.)
      call uzlset ('LABELXT',.FALSE.)
      call uzlset ('LABELXB',.TRUE.)
       if (xr-xl.gt.2000.0)then
        call uxaxdv ('B',200.0,400.0)
        call uxaxdv ('T',200.0,400.0)
       elseif (xr-xl.gt.400.0)then
        call uxaxdv ('B',50.0,100.0)
        call uxaxdv ('T',50.0,100.0)
       else
        call uxaxdv ('B',10.0,50.0)
        call uxaxdv ('T',10.0,50.0)
       end if

      call uxsttl ('B','x(km)',0.0)
      call uysfmt ('(I2)')
      call uzlset ('LABELYR',.FALSE.)
      call uzlset ('LABELYL',.TRUE.)
      call uyaxdv ('L',5.0,10.0)
      call uyaxdv ('R',5.0,10.0)
      call uziset ('IROTCYL',1)
      call uysttl ('L','Height(km)',0.0)


! **** background U and N **** 
      if (lsub) then!=====================================!

      umin=10.*real(int(0.1*minval(ubg)))+5.0
      umax=10.*real(int(0.1*minval(ubg)))-5.0
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

      do ilayer=1,nlayer
       call uulin (2,(/ ubg(ilayer), ubg(ilayer) /), 0.001*(/ zbnds(ilayer), zbnds(ilayer+1) /) )
      end do

       call uuslnt (3)
      do ilayer=1,nlayer
       call uulin (2,(/ umin+fac*(bn2(ilayer)-bn2min), umin+fac*(bn2(ilayer)-bn2min) /), 0.001*(/ zbnds(ilayer), zbnds(ilayer+1) /) )
      end do

!      call uulin  (nlev,ubgz,0.001*zlev)
!      call uulin  (nlev,bn2z,0.001*zlev)

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
      call uxaxdv ('T',1.,2.)
      call uxsttl ('T','N2(s|-2")',0.0)
      call uysfmt ('(I2)')
      call uzlset ('LABELYR',.FALSE.)
      call uzlset ('LABELYL',.FALSE.)

      end if !=====================================!

      call sgtxzv (0.45,0.65,trim(title1),0.020,0,0,5) ! title
      call sgtxzv (0.17,0.65,trim(title2),0.016,0,-1,5) ! title
      call grcls
end do

end do

return
end subroutine quicklook_yz

!============================================!

subroutine quicklook_xy
use setup

integer::iwork(nx,ny)


real(4)::u_snap(nx,ny)
real(4)::w_snap(nx,ny)
real(4)::dispx_snap(nx,ny)
real(4)::dispz_snap(nx,ny)
real(4)::heat_snap(nx,ny)
real(4)::heat_diag_snap(nx,ny)

character*3::cvarc,cvarl,cvars, cvarv
character*20::title1,title2
logical::lsub
namelist /quicklook_nml/ iout, cvarc, cvarl, cvars, cvarv, title1, title2, lsub
namelist /range_nml/ hmin, hmax, wmin, wmax


heat_fact =  grav/300.0 / 86400.0


open(11,file='calc.conf')
 read(11,nml=quicklook_nml)
 read(11,nml=range_nml)
close(11)

ilev=10 !!! sample

!do it=1,nt,1
do it=nt/2,nt/2

write(title2,'(A,I5,A)')'time= ',int(axt(it)-axt(1)),' (s)'

heat_snap=heat(:,it,:,ilev) / heat_fact
heat_diag_snap=heat_diag(:,it,:,ilev) / heat_fact

w_snap=w(:,it,:,ilev)

! write(*,*) it, maxval(dispz_snap),minval(dispz_snap)
 write(*,*) 'heat ',maxval(heat_snap),minval(heat_snap)
 write(*,*) 'w    ',maxval(w_snap),minval(w_snap)


! do iy=1,ny
! write(*,*) axy(iy),heat_snap(nx/2,iy),heat_diag_snap(nx/2,iy)
! end do
! stop

 write(*,*) 'w',maxval(w_snap),minval(w_snap)

wmin=minval(w_snap)
wmax=maxval(w_snap)

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
      call grswnd (1.0e-3*xl,1.0e-3*xr,1.0e-3*yl,1.0e-3*yr) ! set window
      call grsvpt (0.15,0.75,0.13,0.63) ! set viewport !!! 1:1
      call grstrn (1) ! linear or log
      call grstrf


! **** contour ****
      call udiclv
      call udlset ('LMSG',.FALSE.) ! 'contour interval'
      call udrset ('RSIZEL',0.015) ! label
      call udiset ('INDXMJ',13) ! major line
      call udiset ('INDXMN',11) ! minor line
      call uwsgxa ((1.0e-3)*axx,nx) ! yaxis value
      call uwsgya ((1.0e-3)*axy,ny) ! zaxis value (km)
      select case(trim(cvarc))
 case('u')
  call udgcla(umin,umax,-10.0)
  call uddclv(0.0)
  call udcntr (u_snap,nx,nx,ny)
 case('w')
  call udgcla(wmin,wmax,-10.0)
  call uddclv(0.0)
  call udcntr (w_snap,nx,nx,ny)
 end select


! *** heating ***
  call udiclv
  call uwsgxa ((1.0e-3)*axx,nx) ! xaxis value
  call uwsgya ((1.0e-3)*axy,ny) ! yaxis value
  call udiset ('INDXMJ',23) ! major line
  call udiset ('INDXMN',21) ! minor line
  call udgcla(hmin,hmax,-10.0)
  call udcntr(heat_snap,nx,nx,ny)

  call udiclv
  call udiset ('INDXMJ',43) ! major line
  call udiset ('INDXMN',41) ! minor line
  call udgcla(hmin,hmax,-10.0)
  call udcntr(heat_diag_snap,nx,nx,ny)

! **** x ,y axis ****
      call uzinit
      call uziset ('INDEXT2',3)
      call uziset ('INDEXL1',5)
      call uzrset ('RSIZEL1',0.018)
      call uzrset ('RSIZEC1',0.018)
      call uzrset ('RSIZET2',0.010)
      call uzrset ('RSIZET1',0.004)
      call uzlset ('LOFFSET',.FALSE.)
      call uzlset ('LABELXT',.FALSE.)
      call uzlset ('LABELXB',.TRUE.)
       if (xr-xl.gt.2000.0e3)then
        call uxaxdv ('B',200.0,400.0)
        call uxaxdv ('T',200.0,400.0)
       elseif (xr-xl.gt.400.0e3)then
        call uxaxdv ('B',50.0,100.0)
        call uxaxdv ('T',50.0,100.0)
       else
        call uxaxdv ('B',10.0,50.0)
        call uxaxdv ('T',10.0,50.0)
       end if

      call uxsttl ('B','x(km)',0.0)
      call uysfmt ('(I5)')
      call uzlset ('LABELYR',.FALSE.)
      call uzlset ('LABELYL',.TRUE.)

       if (yr-yl.gt.2000.0e3)then
        call uyaxdv ('L',200.0,400.0)
        call uyaxdv ('R',200.0,400.0)
       elseif (yr-yl.gt.400.0e3)then
        call uyaxdv ('L',50.0,100.0)
        call uyaxdv ('R',50.0,100.0)
       else
        call uyaxdv ('L',10.0,50.0)
        call uyaxdv ('R',10.0,50.0)
       end if

      call uziset ('IROTCYL',1)
      call uysttl ('L','y (km)',0.0)



      call sgtxzv (0.45,0.65,trim(title1),0.020,0,0,5) ! title
      call sgtxzv (0.17,0.65,trim(title2),0.016,0,-1,5) ! title
      call grcls
end do

end do

return
end subroutine quicklook_xy

!============================================!
