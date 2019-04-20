program newton_hermite
use output
use func

integer,parameter::nhmax=20

double precision::vx,vx_init
double precision::vxs_half(nhmax/2)
double precision::vxs_half_b(nhmax/2)

double precision::vxs(nhmax)
double precision::w(nhmax)


eps1=0.30
eps2=0.1
nitr=1000

open(11,file='vxs.dat',form='formatted')
open(12,file='vws.dat',form='formatted')
!!! sample 


do ih=2,nhmax

 nv=ih/2 

 do iv=1,nv
 !! initial
  if (ih.eq.2)then
   vx=eps1
  elseif (mod(ih,2).eq.1)then
   if (iv.eq.nv)then
    vx=vxs_half_b(iv) * (1.0+eps2)
   else
    vx=sqrt(vxs_half_b(iv) * vxs_half_b(iv+1))
   end if
  elseif (mod(ih,2).eq.0)then
   if (iv.eq.1)then
    vx=eps1
   elseif (iv.eq.nv)then
    vx=vxs_half_b(iv-1) * (1.0+eps2)
   else
    vx=sqrt(vxs_half_b(iv-1) * vxs_half_b(iv))
   end if
  end if

!!! modified newtonian method
   do itr=1,nitr
!    write(*,*)itr,vx,vf
    vf=hermite_d(vx,ih)
    dfdx=2.0*dble(ih-1)*hermite_d(vx,ih-1) 
    dvx= -vf/dfdx/2.0
    vx=vx+dvx
   end do
  vxs_half(iv)=vx
 end do
write(*,*)'ih=',ih
write(*,*)vxs_half(1:nv)
vxs_half_b=vxs_half

 if(mod(ih,2).eq.0) then
  vxs(1:nv)=-vxs_half(nv:1:-1)
  vxs(nv+1:2*nv)=vxs_half(1:nv)
 else
  vxs(1:nv)=-vxs_half(nv:1:-1)
  vxs(nv+1)=0.0
  vxs(nv+2:2*nv+1)=vxs_half(1:nv)
 end if

 !!! weight
  w=0.0
  do iv=1,nhmax
!   if (vxs(iv).ne.0.0)then
   if (vxs(iv).ne.0.0.or.iv.le.2*nv)then
   w_log = dble(ih-1)*log(2.0) + 0.5*log(3.1416) + log_rfactorial(ih-1) - log(dble(ih)) - 2.0 * log(abs(hermite_d(vxs(iv),ih-1)))
   if (ih.eq.3.and.iv.eq.2) write(*,*) vxs(iv),hermite_d(vxs(iv),ih-1)
   w(iv)=exp(w_log)
!   w(iv)=w_log
   end if
  end do

write(11,'(I3,20F20.15)') ih,vxs
write(12,'(I3,20F20.15)') ih,w
end do

close(11)
close(12)

axy=(/(yl+(real(iy)-0.5)*(yr-yl)/real(ny) ,iy=1,ny)/)
do iy=1,ny
 f(iy)=para_cyl(axy(iy),nhmax)
 f_new(iy)=0.0
end do


title1='sample'
xlabel='y (m)'
iout=1

astics=1.0
amtics=2.0
bstics=0.1
bmtics=0.5
ylabel=''

call draw(real(vxs),nhmax)



end program newton_hermite

!=====================================================!

subroutine draw(vxs,nhmax)
use output

real(4)::vxs(nhmax)
real(4)::zeros(nhmax)


iout=3
zeros=0.0

!vmax=maxval(f)+0.1*abs(maxval(f))
!vmin=minval(f)-0.1*abs(minval(f))

vmax=0.8
vmin=-0.8

write(*,*) 'min,max',vmin,vmax

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
      call uulinz(ny,axy,f_new,3,13)
!      call uulinz(ntmax,(/(dt*real(it-1),it=1,ntmax)/),x_model,3,43)
     call uumrkz (nhmax,vxs,zeros,4,23,0.002)
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
