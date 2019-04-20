program sample

include 'Grid.h'

real(4)::f(nt)
real(4)::vomg0

complex::c_fh(nt)

character*20::title1,title2
character*20::xlabel,ylabel
real(4)::vmin,vmax
real(4)::astics,amtics
real(4)::bstics,bmtics
integer::iout

vomg0=0.01

do it=1,nt
 f(it)=sin(axt(it)*vomg0)
 write(*,*) it,axt(it),f(it)
end do

c_fh=f
call fft_fwd_comp(nt,c_fh,1)
!do it=1,nt
! write(*,*) it,axt(it),c_fh(it)
!end do
!stop



call fft_bwd_comp(nt,c_fh,1)
f=real(c_fh)

title1='sample'
xlabel='t (s)'
iout=1

astics=100.0
amtics=200.0
bstics=0.1
bmtics=0.5
ylabel=''

iout=3
zeros=0.0

vmax=maxval(f)+0.1*abs(maxval(f))
vmin=minval(f)-0.1*abs(minval(f))

!vmax=0.8
!vmin=-0.8

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
      call grswnd (axt(1),axt(nt),vmin,vmax) ! set window
      call grsvpt (0.15,0.85,0.13,0.73) ! set viewport !!! 1:1
      call grstrn (1) ! linear or log
      call grstrf

 ! **** line ****
      call uulinz(nt,axt,f,1,3)
      call uulinz(2,(/axt(1),axt(nt)/),(/0.0,0.0/),3,1)
      call uulinz(2,(/0.0,0.0/),(/vmin,vmax/),3,1)
!      call uulinz(ntmax,(/(dt*real(it-1),it=1,ntmax)/),x_model,3,
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


stop
end program sample
