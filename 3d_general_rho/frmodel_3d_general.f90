!=============================================!

subroutine frmodel_3d_general
use setup

complex::c_vm
complex::c_vm_bnd(nlayer,2)
complex::c_phase_bnd(nlayer)
complex::c_phase(nlev)
complex::c_hk(nx,ny)
complex::c_emzhk(nx,ny,nlev)

complex::c_fact(nx,ny,2*nlayer-1)
complex::c_mat(2*nlayer-1,2*nlayer-1)
complex::c_mat_t(2*nlayer-1,2*nlayer-1)

complex::c_vx(nx,ny,nlev)
complex::c_vu(nx,ny,nlev)
complex::c_vy(nx,ny,nlev)
complex::c_vv(nx,ny,nlev)
complex::c_vw(nx,ny,nlev)

complex::fc_vm

complex,parameter::c_ur=(1.0,0.0)
complex,parameter::c_ui=(0.0,1.0)

!!! h(x,y) --> h(k,l)

c_hk = h * c_ur

call fft_fwd_comp_2d(nx,ny,c_hk,1)


!!! exp(im(k)z)h(k)

do ik=1,nx
do il=1,ny
 
 do ilayer=1,nlayer
  do j=1,2
   if (axk(ik).ne.0.0)then
    c_vm_bnd(ilayer,j)=fc_vm(bn2_bnd(ilayer,j),ubg_bnd(ilayer,j),vbg_bnd(ilayer,j),axk(ik),axl(il),vnu,tbg_bnd(ilayer,j),f0)
   else
    c_vm_bnd(ilayer,j)=0.0
   end if
  end do
 end do


 c_phase = 0.0
 do ilev=1,nlev
  if (ilev.eq.1) then
       c_phase(ilev) = fc_vm(bn2(ilev),ubg(ilev),vbg(ilev),axk(ik),axl(il),vnu,tbg(ilev),f0) * (zlev(ilev)-zbnds(inlay(ilev))) 
       c_phase_bnd(inlay(ilev))=0.0
  elseif (inlay(ilev).eq.inlay(ilev-1))then
       c_phase(ilev) = c_phase(ilev-1) + fc_vm(bn2(ilev-1),ubg(ilev-1),vbg(ilev-1),axk(ik),axl(il),vnu,tbg(ilev-1),f0) * (zlev(ilev)-zlev(ilev-1)) !! forward
  else
       c_phase(ilev) = fc_vm(bn2(ilev),ubg(ilev),vbg(ilev),axk(ik),axl(il),vnu,tbg(ilev),f0) * (zlev(ilev)-zbnds(inlay(ilev))) 
       c_phase_bnd(inlay(ilev)) = c_phase(ilev-1) + fc_vm(bn2(ilev-1),ubg(ilev-1),vbg(ilev-1),axk(ik),axl(il),vnu,tbg(ilev-1),f0) * (zbnds(inlay(ilev))-zlev(ilev-1)) !! forward
  end if
 end do


 !!! coefficient matrix
  c_mat=0.0
 if (axk(ik).eq.0.0.or.axl(il).eq.0.0)then

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
  else
   c_mat(2*(ilayer-1),  2*(ilayer-2)+1) =   exp(  c_ui * c_phase_bnd(ilayer)) * sqrt(c_vm_bnd(ilayer-1,2)/c_vm_bnd(ilayer,1))
   c_mat(2*(ilayer-1),  2*(ilayer-2)+2) =   exp(- c_ui * c_phase_bnd(ilayer)) * sqrt(c_vm_bnd(ilayer-1,2)/c_vm_bnd(ilayer,1))
   c_mat(2*(ilayer-1),  2*(ilayer-1)+1) = - c_ur
   c_mat(2*(ilayer-1)+1,2*(ilayer-2)+1) =   exp(  c_ui * c_phase_bnd(ilayer)) * c_vm_bnd(ilayer,1) * sqrt(c_vm_bnd(ilayer-1,2)/c_vm_bnd(ilayer,1))
   c_mat(2*(ilayer-1)+1,2*(ilayer-2)+2) = - exp(- c_ui * c_phase_bnd(ilayer)) * c_vm_bnd(ilayer,1) * sqrt(c_vm_bnd(ilayer-1,2)/c_vm_bnd(ilayer,1))
   c_mat(2*(ilayer-1)+1,2*(ilayer-1)+1) = - c_vm_bnd(ilayer,2)
    if(ilayer.ne.nlayer)then
     c_mat(2*(ilayer-1)+1,2*(ilayer-1)+2) =   c_vm_bnd(ilayer,2)
     c_mat(2*(ilayer-1),  2*(ilayer-1)+2) = - c_ur
    end if
  end if
 end do

 end if

!write(*,*)
!write(*,*) ik,il
!write(*,*)
!write(*,*) c_phase_bnd
!write(*,*) c_vm_bnd
!write(*,*) c_mat

c_fact=0.0
c_fact(ik,il,1)=c_hk(ik,il)
call sol_comp(2*nlayer-1,c_mat,c_fact(ik,il,:))

 do ilev=1,nlev
   ilayer = inlay(ilev)
   c_vm=fc_vm(bn2(ilev),ubg(ilev),vbg(ilev),axk(ik),axl(il),vnu,tbg(ilev),f0)
   if (c_vm.ne.0.0) then
    amp = (c_vm_bnd(ilayer,2) / c_vm)**0.5
   else
    amp = 1.0
   end if
  
   if (ilayer.eq.nlayer) then !!! radiation condition
    c_emzhk(ik,il,ilev) = c_fact(ik,il,2*(ilayer-1)+1) * exp (c_ui * c_phase(ilev))
   else
    c_emzhk(ik,il,ilev) =  c_fact(ik,il,2*(ilayer-1)+1) * exp (  c_ui * c_phase(ilev)) * amp &
                      + c_fact(ik,il,2*(ilayer-1)+2) * exp (- c_ui * c_phase(ilev)) * amp
   end if
 call plrz_diag(axk(ik),axl(il),ubg(ilev),vbg(ilev),bn2(ilev),c_vm,c_emzhk(ik,il,ilev),c_vx(ik,il,ilev),c_vu(ik,il,ilev),c_vy(ik,il,ilev),c_vv(ik,il,ilev),c_vw(ik,il,ilev),vnu,tbg(ilev),f0)
 end do

end do !! ik
end do !! il


do ilev=2,nlev
 c_emzhk(:,:,ilev) = c_emzhk(:,:,ilev)  * sqrt(rhobg(1)/rhobg(ilev)) !!! density factor correction
 c_vx(:,:,ilev) = c_vx(:,:,ilev)  * sqrt(rhobg(1)/rhobg(ilev)) !!! density factor correction
 c_vu(:,:,ilev) = c_vu(:,:,ilev)  * sqrt(rhobg(1)/rhobg(ilev)) !!! density factor correction
 c_vw(:,:,ilev) = c_vw(:,:,ilev)  * sqrt(rhobg(1)/rhobg(ilev)) !!! density factor correction
end do

c_vu(1:2,1:2,:)=0.0
c_vx(1:2,1:2,:)=0.0
c_vv(1:2,1:2,:)=0.0
c_vy(1:2,1:2,:)=0.0
c_vw(1:2,1:2,:)=0.0
c_emzhk(1,:,:)=0.0
!c_emzhk(nx2-1:nx2,ny2-1:ny2,:)=0.0

call fft_bwd_comp_2d(nx,ny,c_emzhk,nlev)
call fft_bwd_comp_2d(nx,ny,c_vx,nlev)
call fft_bwd_comp_2d(nx,ny,c_vu,nlev)
call fft_bwd_comp_2d(nx,ny,c_vv,nlev)
call fft_bwd_comp_2d(nx,ny,c_vy,nlev)
call fft_bwd_comp_2d(nx,ny,c_vw,nlev)

dispz = real(c_emzhk)
dispx = real(c_vx)
dispy = real(c_vy)
w     = real(c_vw)
u     = real(c_vu)
v     = real(c_vv)

return
end subroutine frmodel_3d_general

!=============================================!
!!! polarization relation
subroutine plrz_diag(vk,vl,vu,vv,vbn2,c_vm,c_vz1,c_vx1,c_vu1,c_vy1,c_vv1,c_vw1,vnu,vt,vf)

real(4)::vt
real(4)::vk,vl,vu,vv,vbn2
 real(4)::vf
complex,intent(in)::c_vm,c_vz1

complex,intent(out)::c_vx1,c_vu1,c_vy1,c_vv1,c_vw1

complex,parameter::c_ur=(1.0,0.0)
complex,parameter::c_ui=(0.0,1.0)

complex::c_oh
real(4)::alpha

include 'Pcon.h'

 !!!! under construction

alpha = 0.5 * grav / rd / vt !!!  1/(2H)

c_oh = - vk * vu - vl * vv + c_ui * vnu

c_vx1=0.0
c_vy1=0.0
if (vk.ne.0.0) c_vx1 = ( - vk - c_ui * vf/c_oh * vl ) * ( c_oh**2 / ( c_oh**2-vf**2 ) ) * ( vbn2 / c_oh**2 -1.0) * c_vz1 !!!
if (vl.ne.0.0) c_vy1 = ( - vl + c_ui * vf/c_oh * vl ) * ( c_oh**2 / ( c_oh**2-vf**2 ) ) * ( vbn2 / c_oh**2 -1.0) * c_vz1 !!!


c_vu1 = - c_ui * c_oh * c_vx1
c_vw1 = - c_ui * c_oh * c_vz1

return
end subroutine plrz_diag
!=============================================!
!!! dispersion relation and sign specification
complex function fc_vm(vbn2,vu,vv,vk,vl,vnu,vt,vf)

 real(4)::vt
 real(4)::vf
 real(4)::vbn2,vu,vv,vk,vl
 complex,parameter::c_ur=(1.0,0.0)
 complex,parameter::c_ui=(0.0,1.0)
 complex::c_oh
 real(4)::alpha

 include 'Pcon.h'

 alpha = 0.5 * grav / rd / vt !!!  1/(2H)

 c_oh = - vk * vu - vl * vv + c_ui * vnu

 if (real(c_oh).eq.0.0)then
  fc_vm=0.0
  return
 end if

fc_vm=0.0
if (c_oh.ne.0.0.and.c_oh**2 -vf**2.ne.0.0) fc_vm =  ( (vk**2 + vl**2) * ( vbn2 * c_ur - c_oh**2 ) / ( c_oh**2 -vf**2 )- alpha**2 ) ** 0.5
!!if (c_oh.ne.0.0) fc_vm =  ( (vk**2 + vl**2) * ( vbn2 * c_ur - c_oh**2 ) / c_oh**2 - alpha**2 ) ** 0.5


 !!! upper limit of abs(vm)

 fc_vm_upper=1.0e-3

 if (abs(aimag(fc_vm)).gt.fc_vm_upper)then
  fc_vm=real(fc_vm)*c_ur + real(sign(fc_vm_upper,aimag(fc_vm))) * c_ui
 end if

 if (real(fc_vm).ne.0.0)then
 !!! propagate
  fc_vm = fc_vm *  ( sign(1.0,real(-c_oh)))
 elseif (aimag(fc_vm).ne.0.0)then
 !!! evanescent
  fc_vm = fc_vm * ( aimag(fc_vm)/abs(aimag(fc_vm)) )
 end if

end function fc_vm

!=============================================!
