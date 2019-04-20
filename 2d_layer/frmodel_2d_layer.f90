!=============================================!

subroutine frmodel_2d_const
use setup

complex::vsgn
complex::c_vm(nlayer)
complex::c_hk(nx)
complex::c_emzhk(nx,nlev)

complex::c_fact(nx,2*nlayer-1)
complex::c_mat(2*nlayer-1,2*nlayer-1)
complex::c_mat_t(2*nlayer-1,2*nlayer-1)

complex::c_vx(nx,nlev)
complex::c_vu(nx,nlev)
complex::c_vw(nx,nlev)

complex::fc_vm

complex,parameter::c_ur=(1.0,0.0)
complex,parameter::c_ui=(0.0,1.0)

!!! h(x) --> h(k)

c_hk = h * c_ur

call fft_fwd_comp(nx,c_hk,1)

!!! exp(im(k)z)h(k)

do ik=1,nx
 
 do ilayer=1,nlayer
  c_vm(ilayer)=fc_vm(bn2(ilayer),ubg(ilayer),axk(ik),vnu)
 end do
 if (axk(ik).eq.0)then

   c_emzhk(ik,:)=0.0
   c_emzhk(ik,1) = c_hk(ik)
 do ilev=1,nlev
   ilayer = inlay(ilev)
   call plrz_diag(axk(ik),ubg(ilayer),bn2(ilayer),c_vm(ilayer),c_emzhk(ik,ilev),c_vx(ik,ilev),c_vu(ik,ilev),c_vw(ik,ilev),vnu)
 end do

 else 

 !!! coefficient matrix
  c_mat=0.0
 do ilayer=1,nlayer
  if(ilayer.eq.1)then 
   c_mat(ilayer,1) = exp(  c_ui * c_vm(ilayer) * zbot)
   c_mat(ilayer,2) = exp(- c_ui * c_vm(ilayer) * zbot)
  elseif(ilayer.eq.nlayer)then
   c_mat(2*(ilayer-1),  2*(ilayer-2)+1) =   exp(  c_ui * c_vm(ilayer-1) * (zbnds(ilayer)-zbnds(ilayer-1)))
   c_mat(2*(ilayer-1),  2*(ilayer-2)+2) =   exp(- c_ui * c_vm(ilayer-1) * (zbnds(ilayer)-zbnds(ilayer-1)))   
   c_mat(2*(ilayer-1),  2*(ilayer-1)+1) = - c_ur
   c_mat(2*(ilayer-1)+1,2*(ilayer-2)+1) =   exp(  c_ui * c_vm(ilayer-1) * (zbnds(ilayer)-zbnds(ilayer-1))) * c_vm(ilayer-1)
   c_mat(2*(ilayer-1)+1,2*(ilayer-2)+2) = - exp(- c_ui * c_vm(ilayer-1) * (zbnds(ilayer)-zbnds(ilayer-1))) * c_vm(ilayer-1)  
   c_mat(2*(ilayer-1)+1,2*(ilayer-1)+1) = - c_vm(ilayer)
  else
   c_mat(2*(ilayer-1),  2*(ilayer-2)+1) =   exp(  c_ui * c_vm(ilayer-1) * (zbnds(ilayer)-zbnds(ilayer-1)))
   c_mat(2*(ilayer-1),  2*(ilayer-2)+2) =   exp(- c_ui * c_vm(ilayer-1) * (zbnds(ilayer)-zbnds(ilayer-1)))   
   c_mat(2*(ilayer-1),  2*(ilayer-1)+1) = - c_ur
   c_mat(2*(ilayer-1),  2*(ilayer-1)+2) = - c_ur
   c_mat(2*(ilayer-1)+1,2*(ilayer-2)+1) =   exp(  c_ui * c_vm(ilayer-1) * (zbnds(ilayer)-zbnds(ilayer-1))) * c_vm(ilayer-1)
   c_mat(2*(ilayer-1)+1,2*(ilayer-2)+2) = - exp(- c_ui * c_vm(ilayer-1) * (zbnds(ilayer)-zbnds(ilayer-1))) * c_vm(ilayer-1)  
   c_mat(2*(ilayer-1)+1,2*(ilayer-1)+1) = - c_vm(ilayer)
   c_mat(2*(ilayer-1)+1,2*(ilayer-1)+2) =   c_vm(ilayer)
  end if
 end do


c_fact=0.0
c_fact(ik,1)=c_hk(ik)
call sol_comp(2*nlayer-1,c_mat,c_fact(ik,:))



 do ilev=1,nlev
   ilayer = inlay(ilev)
   if (ilayer.eq.nlayer) then !!! radiation condition
    c_emzhk(ik,ilev) = c_fact(ik,2*(ilayer-1)+1) * exp (c_ui * c_vm(ilayer) * (zlev(ilev)-zbnds(ilayer))) 
   else
    c_emzhk(ik,ilev) =  c_fact(ik,2*(ilayer-1)+1) * exp (  c_ui * c_vm(ilayer) * (zlev(ilev)-zbnds(ilayer))) &
                      + c_fact(ik,2*(ilayer-1)+2) * exp (- c_ui * c_vm(ilayer) * (zlev(ilev)-zbnds(ilayer))) 
   end if
 
   call plrz_diag(axk(ik),ubg(ilayer),bn2(ilayer),c_vm(ilayer),c_emzhk(ik,ilev),c_vx(ik,ilev),c_vu(ik,ilev),c_vw(ik,ilev),vnu)
 end do

 end if

end do !!! ik




c_vu(1:2,:)=0.0
c_vx(1:2,:)=0.0
call fft_bwd_comp(nx,c_emzhk,nlev)
call fft_bwd_comp(nx,c_vx,nlev)
call fft_bwd_comp(nx,c_vu,nlev)
call fft_bwd_comp(nx,c_vw,nlev)


dispz = real(c_emzhk)
dispx = real(c_vx)
u     = real(c_vu)
w     = real(c_vw)


return
end subroutine frmodel_2d_const

!=============================================!
!!! polarization relation
subroutine plrz_diag(vk,vu,vbn2,c_vm,c_vz1,c_vx1,c_vu1,c_vw1,vnu)

real(4)::vk,vu,vbn2
complex,intent(in)::c_vm,c_vz1

complex,intent(out)::c_vx1,c_vu1,c_vw1

complex,parameter::c_ur=(1.0,0.0)
complex,parameter::c_ui=(0.0,1.0)

complex::c_oh

c_oh = - vk * vu + c_ui * vnu

c_vx1=0.0
if (vk.ne.0.0) c_vx1 = - c_vz1 * c_vm / vk
c_vu1 = - c_ui * c_oh * c_vx1
c_vw1 = - c_ui * c_oh * c_vz1

return
end subroutine plrz_diag
!=============================================!
!!! dispersion relation and sign specification
complex function fc_vm(vbn2,vu,vk,vnu)

 real(4)::vbn2,vu,vk
 complex,parameter::c_ur=(1.0,0.0)
 complex,parameter::c_ui=(0.0,1.0)
 complex::c_oh

 c_oh = - vk * vu + c_ui * vnu

fc_vm=0.0
if (c_oh.ne.0.0) fc_vm =  ( vk**2 * ( vbn2 * c_ur - c_oh**2 ) / c_oh**2 ) ** 0.5



 if (real(fc_vm).ne.0.0)then
 !!! propagate
  fc_vm = fc_vm *  ( sign(1.0,real(-c_oh)))
 elseif (aimag(fc_vm).ne.0.0)then
 !!! evanescent
  fc_vm = fc_vm * ( aimag(fc_vm)/abs(aimag(fc_vm)) )
 end if


end function fc_vm

!=============================================!
