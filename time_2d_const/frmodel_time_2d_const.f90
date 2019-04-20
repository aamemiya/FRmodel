!=============================================!

subroutine frmodel_time_2d_const
use setup

complex::vsgn
complex::c_vm
complex::c_hk(nx,nt)
complex::c_emzhk(nx,nt,nlev)

complex::c_vx(nx,nt,nlev)
complex::c_vu(nx,nt,nlev)
complex::c_vw(nx,nt,nlev)

complex::fc_vm

complex,parameter::c_ur=(1.0,0.0)
complex,parameter::c_ui=(0.0,1.0)

!!! h(x) --> h(k)

c_hk = h * c_ur

call fft_fwd_comp_2d(nx,nt,c_hk,1)

!!! exp(im(k)z)h(k)

do ik=1,nx
do io=1,nt
 do ilev=1,nlev

   c_vm=fc_vm(bn2(ilev),ubg(ilev),axk(ik),axo(io),vnu)

   c_emzhk(ik,io,ilev)=c_hk(ik,io)*exp( c_ui * c_vm * zlev(ilev) )

   call plrz_diag(axk(ik),axo(io),ubg(ilev),bn2(ilev),c_vm,c_emzhk(ik,io,ilev),c_vx(ik,io,ilev),c_vu(ik,io,ilev),c_vw(ik,io,ilev),vnu)
 end do
end do
end do

c_vu(1:2,:,:)=0.0
c_vx(1:2,:,:)=0.0
call fft_bwd_comp_2d(nx,nt,c_emzhk,nlev)
call fft_bwd_comp_2d(nx,nt,c_vx,nlev)
call fft_bwd_comp_2d(nx,nt,c_vu,nlev)
call fft_bwd_comp_2d(nx,nt,c_vw,nlev)


dispz = real(c_emzhk)
dispx = real(c_vx)
u     = real(c_vu)
w     = real(c_vw)


return
end subroutine frmodel_time_2d_const

!=============================================!
!!! polarization relation
subroutine plrz_diag(vk,vo,vu,vbn2,c_vm,c_vz1,c_vx1,c_vu1,c_vw1,vnu)

real(4)::vk,vu,vbn2
complex,intent(in)::c_vm,c_vz1

complex,intent(out)::c_vx1,c_vu1,c_vw1

complex,parameter::c_ur=(1.0,0.0)
complex,parameter::c_ui=(0.0,1.0)

complex::c_oh

c_oh = vo - vk * vu + c_ui * vnu

c_vx1=0.0
if (vk.ne.0.0) c_vx1 = - c_vz1 * c_vm / vk
c_vu1 = - c_ui * c_oh * c_vx1
c_vw1 = - c_ui * c_oh * c_vz1

return
end subroutine plrz_diag
!=============================================!
!!! dispersion relation and sign specification
complex function fc_vm(vbn2,vu,vk,vo,vnu)

 real(4)::vbn2,vu,vk
 complex,parameter::c_ur=(1.0,0.0)
 complex,parameter::c_ui=(0.0,1.0)
 complex::c_oh

 c_oh = vo - vk * vu + c_ui * vnu

!!! fc_vm =  ( ( vbn2 / vu**2 - vk**2 ) * c_ur ) ** 0.5

 fc_vm =  ( vk**2 * ( vbn2 * c_ur - c_oh**2 ) / c_oh**2 ) ** 0.5

 if (real(fc_vm).ne.0.0)then
 !!! propagate
  fc_vm = fc_vm *  ( sign(1.0,real(-c_oh)))
 elseif (aimag(fc_vm).ne.0.0)then
 !!! evanescent
  fc_vm = fc_vm * ( aimag(fc_vm)/abs(aimag(fc_vm)) )
 end if

end function fc_vm

!=============================================!
