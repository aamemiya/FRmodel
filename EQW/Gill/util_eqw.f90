!==============================================================!
module eqw_decomp
use setup

integer,parameter::nh=20

real(4)::ylen

real(4)::axy_gauss(nh)
real(4)::wty_gauss(nh)

integer::nhmax

  contains !===========================!
  subroutine set_grid_hermite !!! set axy_gauss and wty

   call set_gauss_hermite(nh,axy_gauss,wty_gauss)
   do ih=1,nh
      wty_gauss(ih) = wty_gauss(ih) * exp(axy_gauss(ih)**2)
   end do

   !! length scale
   axy_gauss = axy_gauss * ylen
   wty_gauss = wty_gauss * ylen

!   write(*,*) 'check'
!   do ih=1,nh
!    write(*,*) ih,axy_gauss(ih),wty_gauss(ih)
!   end do

   !! available range 
   if (axy_gauss(nh).le.yr)then
    nhmax = nh
   else
    icut=iblkge(axy_gauss,nh,yr)
    nhmax = nh - 2*(nh-icut)
    axy_gauss(1:nhmax) = axy_gauss((nh-icut)+1:icut)
    wty_gauss(1:nhmax) = wty_gauss((nh-icut)+1:icut)
   end if   

   
  end subroutine set_grid_hermite
  !=====================================================!
  subroutine y2g(arrayin,arrayout_gauss) 

    real(4)::arrayin(ny)
    real(4)::arrayout_gauss(nh) 

    do ih=1,nhmax
     iyl=min(iblkge(axy,ny,axy_gauss(ih)),ny-1)
     fact_y= (axy_gauss(ih)-axy(iyl)) / (axy(iyl+1)-axy(iyl))
     arrayout_gauss(ih)=(1.0-fact_y)*arrayin(iyl) + fact_y * arrayin(iyl+1)
    end do


    return
  end subroutine y2g
  !=====================================================!
  subroutine decomp_fwd(arrayin,arrayout) !!! Geometric Grid

    real(4)::arrayin(ny)
    real(4)::arrayin_gauss(nh) 
    real(4)::arrayout(nh)

    call y2g(arrayin,arrayin_gauss)
    call decomp_fwd_gg(arrayin_gauss,arrayout)

    return
  end subroutine decomp_fwd
  !=====================================================!
  subroutine decomp_fwd_gg(arrayin,arrayout) !!! Guassian Grid
   use func
    real(4)::arrayin(nh)
    real(4)::arrayout(nh)

    do ih=0,nhmax-1
       arrayout(ih+1)=0.0
       do iy=1,nhmax
          arrayout(ih+1) = arrayout(ih+1) + wty_gauss(iy) * para_cyl(axy_gauss(iy)/ylen,ih)*arrayin(iy) / ylen
       end do
    end do

    return
  end subroutine decomp_fwd_gg
  !=====================================================!
subroutine decomp_bwd(arrayin,arrayout)
 use func
real(4)::arrayin(nh)
real(4)::arrayout(ny)

do iy=1,ny
 arrayout(iy)=0.0
 do ih=0,nhmax-1
  arrayout(iy) = arrayout(iy) + para_cyl(axy(iy)/ylen,ih)*arrayin(ih+1)
 end do
end do

return
end subroutine decomp_bwd
  !=====================================================!

end module eqw_decomp
!====================================!

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

!====================================!
