!===============================================!
subroutine set_grid_hermite(nh,yl,yr,ylens,nhe,yle,yre,ylene,grid,weight)
use func

integer,intent(in)::nh
real(4),intent(in)::yl,yr
real(4),intent(in)::ylens(nh)
integer,intent(out)::nhe
real(4),intent(out)::yle,yre
real(4),intent(out)::ylene
real(4),intent(out)::grid(nh)
real(4),intent(out)::weight(nh)

real(4)::grid_b(nh)
real(4)::weight_GH(nh)

ih=1
grid=0.0
weight_GH=0.0

fact_taper=1.1

do while (grid(ih)*ylens(ih)*fact_taper.lt.yr.and.ih.lt.nh)
 grid_b=grid
 ih=ih+1
 call set_gauss_hermite(ih,grid(1:ih),weight_GH(1:ih))
end do



 nhe=ih
 if (grid(ih)*ylens(ih)*fact_taper.gt.yr) then
  grid=grid_b
  nhe=nhe-1
  ylene=ylens(ih-1)
 else
  ylene=ylens(ih)
 end if

 weight = ylene * weight_GH * exp(grid**2)
 grid=grid*ylene


if (nhe.lt.nh)then
 yre = yr
 yle = yl
else
 yre =  max(grid(nh),grid(nh-1))*fact_taper
 yle = -yre
end if 

return
end subroutine set_grid_hermite

!===============================================!

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

!=====================================================!
subroutine decomp_fwd_yg_EQW(arrayin,arrayout,ylens)
use param
use func 

real(4)::arrayin(ny)
real(4)::arrayout(nh)
real(4)::arrayin_eqw(nh)
real(4)::arrayout_eqw(nh)
real(4)::ylens(nh)

real(4)::rmat_conv(nh,nh)

do ih=0,nhe-1
 do iy=1,nhe
 arrayin_eqw(iy)= para_cyl(axy_gauss(iy)/ylens(ih+1),ih)
 end do
 call decomp_fwd_gg(arrayin_eqw,arrayout_eqw)
 rmat_conv(:,ih+1)=arrayout_eqw
end do

!do ih=1,nh
! write(*,*)ih,rmat_conv(ih,:)
!end do
! stop


 call decomp_fwd_yg(arrayin,arrayout)
!write(*,*)maxval(arrayin),maxval(arrayout)
 call sol_real(nhe,rmat_conv(1:nhe,1:nhe),arrayout(1:nhe))
!write(*,*)maxval(arrayout)
return
end subroutine decomp_fwd_yg_EQW
!=====================================================!
subroutine decomp_fwd_yg(arrayin,arrayout) !!! Y Grid -> Gaussian Grid
use param
use func 
use spl_psp

real(4)::arrayin(ny)
real(4)::arrayin_gauss(nh)
real(4)::arrayout(nh)


call spl_yg(arrayin,arrayin_gauss)

call decomp_fwd_gg(arrayin_gauss,arrayout)

!do ih=1,nh
! write(*,*) ih, arrayin_gauss(ih), arrayout(ih)
!end do

return
end subroutine decomp_fwd_yg
!====================================================!
subroutine decomp_fwd_gg(arrayin,arrayout) !!! Gaussian Grid
use param
use func 

real(4)::arrayin(nh)
real(4)::arrayout(nh)

do ih=0,nhe-1
 arrayout(ih+1)=0.0
 do iy=1,nhe
  arrayout(ih+1) = arrayout(ih+1) + wty_gauss(iy) * para_cyl(axy_gauss(iy)/ylene,ih)*arrayin(iy) / ylene
 end do
end do

return
end subroutine decomp_fwd_gg

!=====================================================!
subroutine decomp_bwd(arrayin,arrayout)
use param
use func 

real(4)::arrayin(nh)
real(4)::arrayout(ny)

do iy=1,ny
 arrayout(iy)=0.0
 do ih=0,nhe-1
  arrayout(iy) = arrayout(iy) + para_cyl(axy(iy)/ylene,ih)*arrayin(ih+1)
 end do
end do

return
end subroutine decomp_bwd
!=====================================================!
subroutine decomp_bwd_gy_EQW(arrayin,arrayout,ylens) ! Gauss -> Y
use param
use func 

real(4)::arrayin(nh)
real(4)::arrayout(ny)
real(4)::ylens(nh)

do iy=1,ny
 arrayout(iy)=0.0
 do ih=0,nhe-1
  arrayout(iy) = arrayout(iy) + para_cyl(axy(iy)/ylens(ih+1),ih)*arrayin(ih+1)
 end do
end do

return
end subroutine decomp_bwd_gy_EQW
!=====================================================!
