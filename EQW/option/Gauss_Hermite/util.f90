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
