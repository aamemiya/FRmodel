
!==============================================================!
program main
use func
use output

real(4)::ylens(nh)
real(4)::fh_eqw(nh)
real(4)::fh_comp(nh)

axy=(/(yl+(real(iy)-0.5)*(yr-yl)/real(ny) ,iy=1,ny)/)

ylens=(/( ylen*6.0/sqrt(real(2*ih+1)) ,ih=1,nh)/)
!ylens(:)=ylen

do iy=1,ny
 f(iy)=exp(-(real(axy(iy))/(ylen*1.0))**2) * sin(3.1416/(ylen*4.0)*axy(iy))
end do

call decomp_fwd(f,fh)

call decomp_fwd_EQW(f,fh_eqw,ylens)

do ih=1,nh
 write(*,*) ih, fh(ih), fh_eqw(ih),ylens(ih)
end do

!fh_eqw=0.0
!fh_eqw(1)=0.3
!fh_eqw(2)=-0.2
!fh_eqw(3)=0.1

 call decomp_bwd(fh,f_new)

 call decomp_bwd_EQW(fh_eqw,f_new_EQW,ylens)

do ih=1,nh
 fh_comp=0.0
 fh_comp(ih)=fh(ih)
 call decomp_bwd(fh_comp,f_new_comp(:,ih)) 
 fh_comp=0.0
 fh_comp(ih)=fh_eqw(ih)
 call decomp_bwd(fh_comp,f_new_comp_EQW(:,ih)) 
end do

title1='sample'
xlabel='y (m)'
iout=1

astics=1000.0e3
amtics=2000.0e3
bstics=1.0
bmtics=1.0
ylabel=''

call draw




end program main

!=====================================================!

subroutine draw
use output

vmax=maxval(f)+0.1*abs(maxval(f))
!vmin=minval(f)-0.1*abs(minval(f))
vmin=-vmax

write(*,*) 'min,max',vmin,vmax


do ih=0,10

if(ih.eq.0) then
 write(title2,'(A)') 'Total'
else
 write(title2,'(A,I2)') 'WN ',ih-1
end if

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
      if (ih.eq.0)then
       call uulinz(ny,axy,f_new,3,23)
       call uulinz(ny,axy,f_new_eqw,3,43)
      else
       call uulinz(ny,axy,f_new_comp(:,ih),3,23)
       call uulinz(ny,axy,f_new_comp_eqw(:,ih),3,43)
      end if
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

end do

return
end subroutine draw
!=====================================================!
subroutine decomp_fwd_EQW(arrayin,arrayout,ylens)
use param
use func 

real(4)::arrayin(ny)
real(4)::arrayout(nh)
real(4)::arrayin_eqw(ny)
real(4)::arrayout_eqw(nh)
real(4)::ylens(nh)

real(4)::rmat_conv(nh,nh)

do ih=0,nh-1
 do iy=1,ny
 arrayin_eqw(iy)= para_cyl(axy(iy)/ylens(ih+1),ih)
 end do
 call decomp_fwd(arrayin_eqw,arrayout_eqw)
 rmat_conv(:,ih+1)=arrayout_eqw
end do

!do ih=1,nh
! write(*,*)ih,rmat_conv(ih,:)
!end do
! stop


 call decomp_fwd(arrayin,arrayout)
 call sol_real(nh,rmat_conv,arrayout)

return
end subroutine decomp_fwd_EQW
!=====================================================!
subroutine decomp_bwd_EQW(arrayin,arrayout,ylens)
use param
use func 

real(4)::arrayin(nh)
real(4)::arrayout(ny)
real(4)::ylens(nh)


do iy=1,ny
 arrayout(iy)=0.0
 do ih=0,nh-1
  arrayout(iy) = arrayout(iy) + para_cyl(axy(iy)/ylens(ih+1),ih)*arrayin(ih+1)
 end do
end do

return
end subroutine decomp_bwd_EQW
!=====================================================!

subroutine decomp_fwd(arrayin,arrayout)
use param
use func 

real(4)::arrayin(ny)
real(4)::arrayout(nh)

do ih=0,nh-1
 arrayout(ih+1)=0.0
 do iy=1,ny
  arrayout(ih+1) = arrayout(ih+1) + (yr-yl)/real(ny) * para_cyl(axy(iy)/ylen,ih)*arrayin(iy) / ylen
 end do
end do

return
end subroutine decomp_fwd
!=====================================================!
subroutine decomp_bwd(arrayin,arrayout)
use param
use func 

real(4)::arrayin(nh)
real(4)::arrayout(ny)

do iy=1,ny
 arrayout(iy)=0.0
 do ih=0,nh-1
  arrayout(iy) = arrayout(iy) + para_cyl(axy(iy)/ylen,ih)*arrayin(ih+1)
 end do
end do

return
end subroutine decomp_bwd
!=====================================================!
