
!==============================================================!
program main
use param
use func
use output
use spl_psp

real(4)::ylens(nh)
real(4)::fh_eqw(nh)
real(4)::fh_comp(nh)

!!! set axy

axy=(/( yl + (yr-yl) * (real(iy)-0.5)/real(ny) , iy=1,ny)/)

do iy=1,ny
 f(iy)=exp(-(real(axy(iy))/(ylen*2.0))**2)*sin(3.1416/(ylen*0.6)*axy(iy))
end do

!ylene 

ylens=(/( ylen*5.0/sqrt(real(2*ih+1)) ,ih=1,nh)/)
!do ie=1,10
 call set_grid_hermite(nh,yl,yr,ylens,nhe,yle,yre,ylene,axy_gauss,wty_gauss)
! ylene=1.4e6 
 write(*,*) nhe,yle,yre
 write(*,*) yre/ylene* ylens(nhe)
!end do
!stop

call spl_init


!ylens=ylen*2.0

call decomp_bwd_yg(f,fh)
call decomp_bwd_yg_EQW(f,fh_eqw,ylens)

!do ih=1,nh
! write(*,*) ih, fh(ih)
!end do

call decomp_fwd(fh,f_new)
call decomp_fwd_gy_EQW(fh_eqw,f_new_eqw,ylens)

do ih=1,ny
 write(*,*) ih, f_new(ih), f_new_eqw(ih)
end do


title1='sample'
xlabel='y (m)'
iout=1

astics=5.0e6
amtics=1.0e7
bstics=1.0
bmtics=1.0
ylabel=''

call draw




end program main

!=====================================================!

subroutine draw
use output

vmax=maxval(f)+0.1*abs(maxval(f))
vmin=minval(f)-0.1*abs(minval(f))

write(*,*) 'min,max',vmin,vmax

!do iy=1,ny
! write(*,*) iy,axy(iy), f_new_eqw(iy)
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
      call grswnd (yl,yr,vmin,vmax) ! set window
      call grsvpt (0.15,0.85,0.13,0.73) ! set viewport !!! 1:1
      call grstrn (1) ! linear or log
      call grstrf

 ! **** line ****
      call uulinz(ny,axy,f,1,3)
      call uulinz(ny,axy,f_new,3,23)
      call uulinz(ny,axy,f_new_eqw,3,43)
!      call uulinz(ntmax,(/(dt*real(it-1),it=1,ntmax)/),x_model,3,43)

! **** x ,y axis ****
      call uziset ('INDEXT2',3)
      call uziset ('INDEXL1',5)
      call uzrset ('RSIZEL1',0.018)
      call uzrset ('RSIZEC1',0.018)
      call uzrset ('RSIZET2',0.010)
      call uzrset ('RSIZET1',0.004)
!      call uzlset ('LABELXB',.FALSE.)
      call uxaxdv ('B',astics,amtics)
      call uxaxdv ('T',astics,amtics)
      call uxsttl ('B',trim(xlabel),0.0)
      call uzlset ('LABELYR',.FALSE.)
!     call uzlset ('LABELYL',.FALSE.)
      call uyaxdv ('L',bstics,bmtics)
      call uyaxdv ('R',bstics,bmtics)
      call uziset ('IROTCYL',1)
      call uysttl ('L',trim(ylabel),0.0)

      call sgtxzv (0.50,0.76,trim(title1),0.024,0,0,5) ! title
      call sgtxzv (0.17,0.76,trim(title2),0.020,0,-1,5) ! title
      call grcls
end do

return
end subroutine draw
!=====================================================!
