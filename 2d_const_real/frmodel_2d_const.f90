!=============================================!

subroutine frmodel_2d_const
use setup

real(4)::hk(nx)
real(4)::emzhk(nx,nlev)

!!! h(x) --> h(k)

hk=h
call fft_fwd_real(nx,hk,1)

!!! exp(im(k)z)h(k)

do ik=0,nx/2

 vm2=bn2(1)/(ubg(1)**2)-axk(ik)**2

 if (vm2.gt.0.0)then
  vm=sqrt(vm2)
  if (ik.eq.0)then
   emzhk(1,:)=hk(1)*cos(vm*zlev(:))
  elseif (ik.eq.nx/2)then
   emzhk(2,:)=hk(2)*cos(vm*zlev(:))
  else
   emzhk(1+2*ik,:)=hk(1+2*ik)*cos(vm*zlev(:))-hk(2+2*ik)*sin(vm*zlev(:))
   emzhk(2+2*ik,:)=hk(1+2*ik)*sin(vm*zlev(:))+hk(2+2*ik)*cos(vm*zlev(:))
  end if
 else
  vm=sqrt(-vm2)
  if (ik.eq.0)then
   emzhk(1,:)=hk(1)*exp(-vm*zlev(:))
  elseif (ik.eq.nx/2)then
   emzhk(2,:)=hk(2)*exp(-vm*zlev(:))
  else
   emzhk(1+2*ik,:)=hk(1+2*ik)*exp(-vm*zlev(:))
   emzhk(2+2*ik,:)=hk(2+2*ik)*exp(-vm*zlev(:))
 end if
 end if
end do


call fft_bwd_real(nx,emzhk,nlev)

dispz=emzhk


return
end subroutine frmodel_2d_const

!=============================================!
