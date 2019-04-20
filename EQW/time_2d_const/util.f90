
!====================================================!

subroutine set_grid_hermite(ny,ylen,yl,yr,grid,weight)
use func

integer,intent(in)::ny
real(4),intent(in)::ylen
real(4),intent(out)::yl,yr
real(4),intent(out)::grid(ny)
real(4),intent(out)::weight(ny)

real(4)::weight_GH(ny)

double precision::vx,vx_init
double precision::vxs_half(ny/2)
double precision::vxs_half_b(ny/2)


call set_gauss_hermite(ny,grid,weight_GH)

do iy=1,ny
 weight(iy) = weight_GH(iy) * exp(grid(iy)**2)
end do

!! length scale
grid = grid*ylen
weight  = weight*ylen

!do iy=1,ny
! write(*,*) iy,grid(iy)*1.0e-3,weight(iy)*1.0e-3
!end do

yr =  grid(ny)*1.1
yl = -yr
write(*,*) yl,yr
yr = tdg(yr) 
yl = -yr
write(*,*) yl,yr

return
end subroutine set_grid_hermite

!====================================================!

subroutine set_gauss_hermite(ndim,point,weight)
use func

integer,intent(in)::ndim
real(4),intent(out)::point(ndim)
real(4),intent(out)::weight(ndim)

double precision::vx,vx_init
double precision::vxs_half(ndim/2)
double precision::vxs_half_b(ndim/2)

eps1=0.30
eps2=0.1
nitr=1000

do ih=2,ndim

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
 end do !!! iv

!write(*,*)'ih=',ih
!write(*,*)vxs_half(1:nv)

vxs_half_b=vxs_half

if (ih.eq.ndim)then
 if (mod(ndim,2).eq.0)then
  point(1:ndim/2) = -real(vxs_half(ndim/2:1:-1))
  point(ndim/2+1:ndim) = real(vxs_half)
 else
  point(1:ndim/2) = -real(vxs_half(ndim/2:1:-1))
  point(ndim/2+1)=0
  point(ndim/2+2:ndim) = real(vxs_half)
 end if


 !!! weight
  weight=0.0
  do iv=1,ndim
   w_log = dble(ih-1)*log(2.0) + 0.5*log(3.1416) + log_rfactorial(ih-1) - log(dble(ih)) - 2.0 * log(abs(hermite_d(dble(point(iv)),ih-1)))
   weight(iv)=exp(w_log)
  end do  

end if !!! ih.eq.nhmax

end do !!! ih


return
end subroutine set_gauss_hermite
!====================================================!
subroutine decomp_bwd_yg_EQW(arrayin,arrayout,ylens)
use setup
use func 

real(4)::arrayin(ny)
real(4)::arrayout(nh)
real(4)::arrayin_eqw(nh)
real(4)::arrayout_eqw(nh)
real(4)::ylens(nh)

real(4)::rmat_conv(nh,nh)

do ih=0,nh-1
 do iy=1,nh
 arrayin_eqw(iy)= para_cyl(axy_gauss(iy)/ylens(ih+1),ih)
 end do
 call decomp_bwd_gg(arrayin_eqw,arrayout_eqw)
 rmat_conv(:,ih+1)=arrayout_eqw
end do

!do ih=1,nh
! write(*,*)ih,rmat_conv(ih,:)
!end do
! stop


 call decomp_bwd_yg(arrayin,arrayout)
!write(*,*)maxval(arrayin),maxval(arrayout)
 call sol_real(nh,rmat_conv,arrayout)
!write(*,*)maxval(arrayout)
return
end subroutine decomp_bwd_yg_EQW
!=====================================================!
subroutine decomp_bwd_yg(arrayin,arrayout) !!! Y Grid -> Guassian Grid
use setup
use func 
use spl_psp

real(4)::arrayin(ny)
real(4)::arrayin_gauss(ny)
real(4)::arrayout(nh)

call spl_yg(arrayin,arrayin_gauss)

call decomp_bwd_gg(arrayin_gauss,arrayout)

return
end subroutine decomp_bwd_yg
!====================================================!
subroutine decomp_bwd_gg(arrayin,arrayout) !!! Guassian Grid
use setup
use func 

real(4)::arrayin(nh)
real(4)::arrayout(nh)

do ih=0,nh-1
 arrayout(ih+1)=0.0
 do iy=1,nh
  arrayout(ih+1) = arrayout(ih+1) + wty_gauss(iy) * para_cyl(axy_gauss(iy)/ylen,ih)*arrayin(iy) / ylen
 end do
end do

return
end subroutine decomp_bwd_gg
!=====================================================!
subroutine decomp_bwd(arrayin,arrayout)
use setup
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
end subroutine decomp_bwd
!=====================================================!
subroutine decomp_fwd_gy_EQW(arrayin,arrayout,ylens) ! Gauss -> Y
use setup
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
end subroutine decomp_fwd_gy_EQW
!=====================================================!
subroutine decomp_fwd_gy(arrayin,arrayout) ! Gauss -> Y
use setup
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
end subroutine decomp_fwd_gy
!=====================================================!
subroutine decomp_fwd(arrayin,arrayout) ! Gauss -> Gauss
use setup
use func 

real(4)::arrayin(nh)
real(4)::arrayout(nh)

do iy=1,nh
 arrayout(iy)=0.0
 do ih=0,nh-1
  arrayout(iy) = arrayout(iy) + para_cyl(axy_gauss(iy)/ylen,ih)*arrayin(ih+1)
 end do
end do

return
end subroutine decomp_fwd
!=====================================================!
subroutine c_decomp_bwd_yg_EQW(c_arrayin,c_arrayout,ylens)
use setup

complex::c_arrayin(ny)
complex::c_arrayout(nh)
real(4)::ylens(nh)

real(4)::r_arrayin(ny,2)
real(4)::r_arrayout(nh,2)

complex,parameter::c_ur=(1.0,0.0)
complex,parameter::c_ui=(0.0,1.0)

r_arrayin(:,1) = real(c_arrayin)
r_arrayin(:,2) =aimag(c_arrayin)

call decomp_bwd_yg_EQW(r_arrayin(:,1),r_arrayout(:,1),ylens)
call decomp_bwd_yg_EQW(r_arrayin(:,2),r_arrayout(:,2),ylens)

c_arrayout = r_arrayout(:,1) * c_ur + r_arrayout(:,2) * c_ui

return
end subroutine c_decomp_bwd_yg_EQW
!=====================================================!
subroutine c_decomp_fwd_gy_EQW(c_arrayin,c_arrayout,ylens)
use setup

complex::c_arrayin(nh)
complex::c_arrayout(ny)
real(4)::ylens(nh)

real(4)::r_arrayin(nh,2)
real(4)::r_arrayout(ny,2)

complex,parameter::c_ur=(1.0,0.0)
complex,parameter::c_ui=(0.0,1.0)

r_arrayin(:,1) = real(c_arrayin)
r_arrayin(:,2) =aimag(c_arrayin)

call decomp_fwd_gy_EQW(r_arrayin(:,1),r_arrayout(:,1),ylens)
call decomp_fwd_gy_EQW(r_arrayin(:,2),r_arrayout(:,2),ylens)

c_arrayout = r_arrayout(:,1) * c_ur + r_arrayout(:,2) * c_ui

return
end subroutine c_decomp_fwd_gy_EQW
!=====================================================!

