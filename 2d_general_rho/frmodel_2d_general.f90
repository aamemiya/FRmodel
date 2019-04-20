!=============================================!

subroutine frmodel_2d_const
use setup

complex::c_vm
complex::c_vm_bnd(nlayer,2)
complex::c_phase_bnd(nlayer)
complex::c_phase(nlev)
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
  do j=1,2
   c_vm_bnd(ilayer,j)=fc_vm(bn2_bnd(ilayer,j),ubg_bnd(ilayer,j),axk(ik),vnu,tbg_bnd(ilayer,j))
  end do
 end do


 c_phase = 0.0
 do ilev=1,nlev
  if (ilev.eq.1) then
       c_phase(ilev) = fc_vm(bn2(ilev),ubg(ilev),axk(ik),vnu,tbg(ilev)) * (zlev(ilev)-zbnds(inlay(ilev))) 
       c_phase_bnd(inlay(ilev))=0.0
  elseif (inlay(ilev).eq.inlay(ilev-1))then
       c_phase(ilev) = c_phase(ilev-1) + fc_vm(bn2(ilev-1),ubg(ilev-1),axk(ik),vnu,tbg(ilev-1)) * (zlev(ilev)-zlev(ilev-1)) !! forward
  else
       c_phase(ilev) = fc_vm(bn2(ilev),ubg(ilev),axk(ik),vnu,tbg(ilev)) * (zlev(ilev)-zbnds(inlay(ilev))) 
       c_phase_bnd(inlay(ilev)) = c_phase(ilev-1) + fc_vm(bn2(ilev-1),ubg(ilev-1),axk(ik),vnu,tbg(ilev-1)) * (zbnds(inlay(ilev))-zlev(ilev-1)) !! forward
  end if
 end do

 !!! coefficient matrix
  c_mat=0.0
 if (axk(ik).eq.0)then

 do ilayer=1,nlayer
  if(ilayer.eq.1)then 
   c_mat(ilayer,1) = exp(  c_ui * c_vm_bnd(ilayer,2) * zbot)
   c_mat(ilayer,2) = exp(- c_ui * c_vm_bnd(ilayer,2) * zbot)
  elseif(ilayer.eq.nlayer)then
   c_mat(2*(ilayer-1),  2*(ilayer-2)+1) =   exp(  c_ui * c_phase_bnd(ilayer))
   c_mat(2*(ilayer-1),  2*(ilayer-2)+2) =   exp(- c_ui * c_phase_bnd(ilayer))   
   c_mat(2*(ilayer-1),  2*(ilayer-1)+1) = - c_ur
   c_mat(2*(ilayer-1)+1,2*(ilayer-2)+1) =   exp(  c_ui * c_phase_bnd(ilayer)) * c_ur
   c_mat(2*(ilayer-1)+1,2*(ilayer-2)+2) = - exp(- c_ui * c_phase_bnd(ilayer)) * c_ur
   c_mat(2*(ilayer-1)+1,2*(ilayer-1)+1) = - c_ur
  else
   c_mat(2*(ilayer-1),  2*(ilayer-2)+1) =   exp(  c_ui * c_phase_bnd(ilayer))
   c_mat(2*(ilayer-1),  2*(ilayer-2)+2) =   exp(- c_ui * c_phase_bnd(ilayer))   
   c_mat(2*(ilayer-1),  2*(ilayer-1)+1) = - c_ur
   c_mat(2*(ilayer-1),  2*(ilayer-1)+2) = - c_ur
   c_mat(2*(ilayer-1)+1,2*(ilayer-2)+1) =   exp(  c_ui * c_phase_bnd(ilayer)) * c_ur
   c_mat(2*(ilayer-1)+1,2*(ilayer-2)+2) = - exp(- c_ui * c_phase_bnd(ilayer)) * c_ur
   c_mat(2*(ilayer-1)+1,2*(ilayer-1)+1) = - c_ur
   c_mat(2*(ilayer-1)+1,2*(ilayer-1)+2) =   c_ur
  end if
 end do

 else 

 do ilayer=1,nlayer
  if(ilayer.eq.1)then 
   c_mat(ilayer,1) = exp(  c_ui * c_vm_bnd(ilayer,2) * zbot)
   c_mat(ilayer,2) = exp(- c_ui * c_vm_bnd(ilayer,2) * zbot)
  elseif(ilayer.eq.nlayer)then
   c_mat(2*(ilayer-1),  2*(ilayer-2)+1) =   exp(  c_ui * c_phase_bnd(ilayer)) * sqrt(c_vm_bnd(ilayer-1,2)/c_vm_bnd(ilayer,1))
   c_mat(2*(ilayer-1),  2*(ilayer-2)+2) =   exp(- c_ui * c_phase_bnd(ilayer)) * sqrt(c_vm_bnd(ilayer-1,2)/c_vm_bnd(ilayer,1))
   c_mat(2*(ilayer-1),  2*(ilayer-1)+1) = - c_ur
   c_mat(2*(ilayer-1)+1,2*(ilayer-2)+1) =   exp(  c_ui * c_phase_bnd(ilayer)) * c_vm_bnd(ilayer,1) * sqrt(c_vm_bnd(ilayer-1,2)/c_vm_bnd(ilayer,1))
   c_mat(2*(ilayer-1)+1,2*(ilayer-2)+2) = - exp(- c_ui * c_phase_bnd(ilayer)) * c_vm_bnd(ilayer,1) * sqrt(c_vm_bnd(ilayer-1,2)/c_vm_bnd(ilayer,1))
   c_mat(2*(ilayer-1)+1,2*(ilayer-1)+1) = - c_vm_bnd(ilayer,2)
  else
   c_mat(2*(ilayer-1),  2*(ilayer-2)+1) =   exp(  c_ui * c_phase_bnd(ilayer)) * sqrt(c_vm_bnd(ilayer-1,2)/c_vm_bnd(ilayer,1))
   c_mat(2*(ilayer-1),  2*(ilayer-2)+2) =   exp(- c_ui * c_phase_bnd(ilayer)) * sqrt(c_vm_bnd(ilayer-1,2)/c_vm_bnd(ilayer,1))
   c_mat(2*(ilayer-1),  2*(ilayer-1)+1) = - c_ur
   c_mat(2*(ilayer-1),  2*(ilayer-1)+2) = - c_ur
   c_mat(2*(ilayer-1)+1,2*(ilayer-2)+1) =   exp(  c_ui * c_phase_bnd(ilayer)) * c_vm_bnd(ilayer,1) * sqrt(c_vm_bnd(ilayer-1,2)/c_vm_bnd(ilayer,1))
   c_mat(2*(ilayer-1)+1,2*(ilayer-2)+2) = - exp(- c_ui * c_phase_bnd(ilayer)) * c_vm_bnd(ilayer,1) * sqrt(c_vm_bnd(ilayer-1,2)/c_vm_bnd(ilayer,1))
   c_mat(2*(ilayer-1)+1,2*(ilayer-1)+1) = - c_vm_bnd(ilayer,2)
   c_mat(2*(ilayer-1)+1,2*(ilayer-1)+2) =   c_vm_bnd(ilayer,2)
  end if
 end do

 end if

c_fact=0.0
c_fact(ik,1)=c_hk(ik)
call sol_comp(2*nlayer-1,c_mat,c_fact(ik,:))

 do ilev=1,nlev
   ilayer = inlay(ilev)
   c_vm=fc_vm(bn2(ilev),ubg(ilev),axk(ik),vnu,tbg(ilev))
   if (c_vm.ne.0.0) then
    amp = (c_vm_bnd(ilayer,2) / c_vm)**0.5
   else
    amp = 1.0
   end if
  
   if (ilayer.eq.nlayer) then !!! radiation condition
    c_emzhk(ik,ilev) = c_fact(ik,2*(ilayer-1)+1) * exp (c_ui * c_phase(ilev))
   else
    c_emzhk(ik,ilev) =  c_fact(ik,2*(ilayer-1)+1) * exp (  c_ui * c_phase(ilev)) * amp &
                      + c_fact(ik,2*(ilayer-1)+2) * exp (- c_ui * c_phase(ilev)) * amp
   end if
 
   call plrz_diag(axk(ik),ubg(ilev),bn2(ilev),c_vm,c_emzhk(ik,ilev),c_vx(ik,ilev),c_vu(ik,ilev),c_vw(ik,ilev),vnu,tbg(ilev))
 end do


end do


do ilev=2,nlev
 c_emzhk(:,ilev) = c_emzhk(:,ilev)  * sqrt(rhobg(1)/rhobg(ilev)) !!! density factor correction
 c_vx(:,ilev) = c_vx(:,ilev)  * sqrt(rhobg(1)/rhobg(ilev)) !!! density factor correction
 c_vu(:,ilev) = c_vu(:,ilev)  * sqrt(rhobg(1)/rhobg(ilev)) !!! density factor correction
 c_vw(:,ilev) = c_vw(:,ilev)  * sqrt(rhobg(1)/rhobg(ilev)) !!! density factor correction
end do

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
subroutine plrz_diag(vk,vu,vbn2,c_vm,c_vz1,c_vx1,c_vu1,c_vw1,vnu,vt)

real(4)::vt
real(4)::vk,vu,vbn2
complex,intent(in)::c_vm,c_vz1

complex,intent(out)::c_vx1,c_vu1,c_vw1

complex,parameter::c_ur=(1.0,0.0)
complex,parameter::c_ui=(0.0,1.0)

complex::c_oh
real(4)::alpha

include 'Pcon.h'

alpha = 0.5 * grav / rd / vt !!!  1/(2H)

c_oh = - vk * vu + c_ui * vnu

c_vx1=0.0
if (vk.ne.0.0) c_vx1 = ( c_vm + alpha * c_ui ) / vk  * c_vz1 
c_vu1 = - c_ui * c_oh * c_vx1
c_vw1 = - c_ui * c_oh * c_vz1

return
end subroutine plrz_diag
!=============================================!
!!! dispersion relation and sign specification
complex function fc_vm(vbn2,vu,vk,vnu,vt)

 real(4)::vt
 real(4)::vbn2,vu,vk
 complex,parameter::c_ur=(1.0,0.0)
 complex,parameter::c_ui=(0.0,1.0)
 complex::c_oh
 real(4)::alpha

 include 'Pcon.h'


 alpha = 0.5 * grav / rd / vt !!!  1/(2H)

 c_oh = - vk * vu + c_ui * vnu

fc_vm=0.0
if (c_oh.ne.0.0) fc_vm =  ( vk**2 * ( vbn2 * c_ur - c_oh**2 ) / c_oh**2 - alpha**2 ) ** 0.5



 if (real(fc_vm).ne.0.0)then
 !!! propagate
  fc_vm = fc_vm *  ( sign(1.0,real(-c_oh)))
 elseif (aimag(fc_vm).ne.0.0)then
 !!! evanescent
  fc_vm = fc_vm * ( aimag(fc_vm)/abs(aimag(fc_vm)) )
 end if


end function fc_vm

!=============================================!
