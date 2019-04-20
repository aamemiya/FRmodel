program test

include 'Pcon.h'

integer,parameter::nt=32,nx=32,nz=10
real(4)::v(nt,nx,nz)
complex::c_v(nt,nx,nz)

complex,parameter::c_ur=(1.0,0.0)
complex,parameter::c_ui=(0.0,1.0)

do ix=1,nx
 v(:,ix,:)=cos(2.0*pi*real(ix-1)/real(nx))+sin(2.0*pi*real(ix-1)/real(nx))
 write(*,*) ix, v(1,ix,1)
end do

c_v = c_ur * v

call fft_bwd_comp_2d(nt,nx,c_v,nz)

 write(*,*)
 do ix=1,nx
 write(*,*) ix, c_v(1,ix,1)
 end do

call fft_fwd_comp_2d(nt,nx,c_v,nz)

 write(*,*)
 do ix=1,nx
 write(*,*) ix, c_v(1,ix,1)
 end do


stop
end program test
