!=============================================!

subroutine frmodel_time_EQ_const
use setup

complex::vsgn
complex::c_vm(nh)
complex::c_vm_KV
complex::c_hk(nx,nt,ny)
complex::c_hk_h(nx,nt,nh)
complex::c_emzhk(nx,nt,nh,nlev)
complex::c_emzhk_y(nx,nt,ny,nlev)

complex::c_work_ny(ny)
complex::c_work_nh(nh)

complex::c_vx(nx,nt,nlev)
complex::c_vu(nx,nt,nlev)
complex::c_vw(nx,nt,nlev)

complex::fc_vm

real(4)::ylens(nh)

complex,parameter::c_ur=(1.0,0.0)
complex,parameter::c_ui=(0.0,1.0)

real(4),parameter::beta=2.0*eomg/er*cos(drad*10.0)

!!! h(x) --> h(k)

c_hk = h * c_ur

call fft_fwd_comp_2d(nx,nt,c_hk,ny)

do ik=1,nx
do io=1,nt

!do io=1,1
!!! Kelvin waves
!   write(*,*) 'Loop'
   c_vm(1)  = fc_vm(bn2(1),ubg(1),axk(ik),axo(io),-1,vnu,ier)
   if (ier.eq.1)then
      c_vm(1)  = fc_vm(bn2(1),ubg(1),-axk(ik),axo(io),-1,vnu,ier)
   end if
   if (ier.eq.0) then
    ylens(1)=real((sqrt(bn2(1))/beta/c_vm(1))**0.5)
    ylens=ylens(1)
    
    write(*,*) ylens(1), ylen*0.1

    if (ylens(1).ge.ylen*0.1)call c_decomp_bwd_yg_EQW(c_hk(ik,io,:),c_hk_h(ik,io,:),ylens)
    write(*,*) axk(ik),axo(io),ylens(1),c_vm(1),io
!    call c_decomp_bwd_yg_EQW(c_hk(ik,io,:),c_hk_h(ik,io,:),ylens)
    if (io.eq.1) write(*,*) 'after DCP',c_hk(ik,io,ny/2),c_hk_h(ik,io,1)
stop
  c_hk_h(ik,io,2:nh)=0.0
   do ilev=1,nlev
    c_emzhk(ik,io,1,ilev)=c_hk_h(ik,io,1)*exp( c_ui * c_vm(1) * zlev(ilev) )
!    write(*,*) c_vm(1) * zlev(ilev)
!   call plrz_diag(axk(ik),axo(io),ubg(ilev),bn2(ilev),c_vm,c_emzhk(ik,io,ilev),c_vx(ik,io,ilev),c_vu(ik,io,ilev),c_vw(ik,io,ilev),vnu)
   end do
   end if
end do
end do

write(*,*) 'after the loop'
write(*,*) maxval(abs(c_hk_h)),minval(abs(c_hk_h))
write(*,*) maxval(abs(c_emzhk(:,:,:,1))),minval(abs(c_emzhk(:,:,:,1)))
stop

 c_emzhk_y=0.0
 c_hk=0.0

do ik=1,nx
do io=1,nt
 c_vm(1) = fc_vm(bn2(1),ubg(1),axk(ik),axo(io),-1,vnu,ier)
   if (ier.eq.1)then
      c_vm(1)  = fc_vm(bn2(1),ubg(1),-axk(ik),axo(io),-1,vnu,ier)
   end if
 if (ier.eq.0) then
  ylens(1)=real((sqrt(bn2(1))/beta/c_vm(1))**0.5)
!  write(*,*) 'ylen,zlen',axk(ik),axo(io),ylens(1),2.0*pi/real(c_vm(1))
  ylens=ylens(1)
  do ilev=1,nlev
  c_work_nh=c_emzhk(ik,io,:,ilev)
  if (ylens(1).ge.ylen*0.3)   call c_decomp_fwd_gy_EQW(c_work_nh,c_work_ny,ylens)
!  if (ylens(1).ge.ylen*0.3)   call c_decomp_fwd_gy_EQW(c_hk_h(ik,io,:),c_work_ny,ylens)
  c_emzhk_y(ik,io,:,ilev)=c_work_ny
  end do
!  c_hk_h(ik,io,2:nh)=0.0
!  if (ylens(1).ge.ylen*0.5)  call c_decomp_fwd_gy_EQW(c_hk_h(ik,io,:),c_hk(ik,io,:),ylens)
   if (ylens(1).ge.ylen*0.3)  call c_decomp_fwd_gy_EQW(c_hk_h(ik,io,:),c_hk(ik,io,:),ylens)
 end if
end do
end do

write(*,*) maxval(real(c_emzhk_y(:,:,:,1))),minval(real(c_emzhk_y(:,:,:,1)))
write(*,*) maxval(real(c_hk(:,:,:))),minval(real(c_hk(:,:,:)))

!c_vu(1:2,:,:)=0.0
!c_vx(1:2,:,:)=0.0
call fft_bwd_comp_2d(nx,nt,c_emzhk_y,ny*nlev)
call fft_bwd_comp_2d(nx,nt,c_hk,ny)
!call fft_bwd_comp_2d(nx,nt,c_vx,nlev)
!call fft_bwd_comp_2d(nx,nt,c_vu,nlev)
!call fft_bwd_comp_2d(nx,nt,c_vw,nlev)



h_diag= real(c_hk)
dispz = real(c_emzhk_y)
!dispx = real(c_vx)
!u     = real(c_vu)
!w     = real(c_vw)


return
end subroutine frmodel_time_EQ_const

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
complex function fc_vm(vbn2,vu,vk,vo,ih,vnu,ier)

 integer::ih !!! WN + 1
 real(4)::vbn2,vu,vk
 complex,parameter::c_ur=(1.0,0.0)
 complex,parameter::c_ui=(0.0,1.0)
 complex::c_oh

 c_oh = vo - vk * vu + c_ui * vnu

!!! fc_vm =  ( ( vbn2 / vu**2 - vk**2 ) * c_ur ) ** 0.5

 if (ih.eq.-1 .and. real(vk/c_oh).gt.0.0 .or. & !!! Kelvin wave
     ih.ge.0  .and. real(vk/c_oh).lt.0.0 ) then !!! Rossby wave 
   fc_vm = ( vbn2 * (vk/ c_oh)**2 * abs(real(2*ih+1))) **0.5

   if (real(fc_vm).ne.0.0)then   !!! propagate
    fc_vm = fc_vm *  ( sign(1.0,real(-c_oh)))
   end if
    if (aimag(fc_vm).ne.0.0)then   !!! evanescent
    fc_vm = fc_vm * ( aimag(fc_vm)/abs(aimag(fc_vm)) )
   end if

   ier=0
!   write(*,*)'ier',ier
  else !!! Prohibited mode
   ier=1
!   write(*,*)'ier',ier
  end if

return
end function fc_vm

!=============================================!
