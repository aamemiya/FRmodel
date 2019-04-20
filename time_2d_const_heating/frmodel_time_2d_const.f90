!=============================================!

subroutine frmodel_time_2d_const
use setup

real(4)::heat_j(nx,nt,nj)
complex::c_heat_ko_j(nx,nt,nj)

complex::c_vw_hlayer(nx,nt,nj)

complex::c_nj_work(nj)
real(4)::r_2_nj_work(2,nj)
complex::c_upb_const


complex::vsgn
complex::c_vm(nlayer)

complex::c_vw(nx,nt,nlev)
complex::c_vu(nx,nt,nlev)
complex::c_vz(nx,nt,nlev)
complex::c_vx(nx,nt,nlev)

complex::c_fact(2*nlayer-2)
complex::c_mat(2*nlayer-2,2*nlayer-2)
complex::c_mat_t(2*nlayer-2,2*nlayer-2)


complex::fc_vm

complex,parameter::c_ur=(1.0,0.0)
complex,parameter::c_ui=(0.0,1.0)

!!! heat(x,t,z) --> heat(k,omega,j)


heat_j(:,:,:)=heat(:,:,1:nj)

 call fft_fwd_sin_trans(nx*nt,heat_j,nj)

c_heat_ko_j = heat_j * c_ur

 call fft_fwd_comp_2d(nx,nt,c_heat_ko_j,nj)


!!! for each set of k and omega...

do ik=1,nx
do io=1,nt
if (axo(io).ne.0.0.or.axk(ik)*ubg(1).ne.0.0)then

 do ilayer=1,nlayer
  c_vm(ilayer)=fc_vm(bn2(ilayer),ubg(ilayer),axk(ik),axo(io),vnu)
 end do

 !!! specific solution in the heating layer
 c_upb_const = 0.0
 do ij=1,nj !!! vertical WN space
!!  c_nj_work(ij) = c_heat_ko_j(ik,io,ij)/(c_vm(1)**2-(pi*real(ij)/(zbnds(2)-zbnds(1)))**2)
  if(axo(io).ne.0.0) c_nj_work(ij) = c_heat_ko_j(ik,io,ij)/(c_vm(1)**2-(pi*real(ij)/(zbnds(2)-zbnds(1)))**2) * (1.0*axk(ik)/(axo(io)-axk(ik)*ubg(1)))**2
  c_upb_const=c_upb_const + real(ij)*pi/(zbnds(2)-zbnds(1)) * c_nj_work(ij) * real(1-2*mod(ij,2))
 end do
 r_2_nj_work(1,:)=real(c_nj_work(:))
 r_2_nj_work(2,:)=aimag(c_nj_work(:))

 call fft_bwd_sin_trans(2,r_2_nj_work,nj)
 c_vw_hlayer(ik,io,:) = c_ur * r_2_nj_work(1,:) + c_ui * r_2_nj_work(2,:)

!!! heat(ik,io,:)=c_vw_hlayer(ik,io,:) !!! test


 !!! coefficient matrix
  c_mat=0.0
 do ilayer=1,nlayer-1
  if (ilayer.eq.1)then
   c_mat(2*(ilayer-1)+1,1) = sin(c_vm(ilayer) * (zbnds(2)-zbnds(1)))
   c_mat(2*(ilayer-1)+1,2*(ilayer)  )   = -c_ur   
   if(nlayer.gt.2) c_mat(2*(ilayer-1)+1,2*(ilayer)+1)   = -c_ur   
   c_mat(2*(ilayer-1)+2,1) = c_vm(ilayer)*cos(c_vm(ilayer) * (zbnds(2)-zbnds(1))) 
   c_mat(2*(ilayer-1)+2,2*(ilayer)  )   = -c_ui * c_vm(ilayer+1)
   if(nlayer.gt.2)c_mat(2*(ilayer-1)+2,2*(ilayer)+1)   =  c_ui * c_vm(ilayer+1)
  elseif(ilayer.eq.nlayer-1)then
   c_mat(2*(ilayer-1)+1,2*(ilayer-1))   = exp( c_ui*c_vm(ilayer)*(zbnds(ilayer+1)-zbnds(ilayer)))
   c_mat(2*(ilayer-1)+1,2*(ilayer-1)+1) = exp(-c_ui*c_vm(ilayer)*(zbnds(ilayer+1)-zbnds(ilayer)))
   c_mat(2*(ilayer-1)+1,2*(ilayer)  )   = -c_ur   
   c_mat(2*(ilayer-1)+2,2*(ilayer-1))   = exp( c_ui*c_vm(ilayer)*(zbnds(ilayer+1)-zbnds(ilayer))) * c_vm(ilayer)
   c_mat(2*(ilayer-1)+2,2*(ilayer-1)+1) = -exp(-c_ui*c_vm(ilayer)*(zbnds(ilayer+1)-zbnds(ilayer))) * c_vm(ilayer)
   c_mat(2*(ilayer-1)+2,2*(ilayer)  )   = -c_vm(ilayer+1)
  else
   c_mat(2*(ilayer-1)+1,2*(ilayer-1))   = exp( c_ui*c_vm(ilayer)*(zbnds(ilayer+1)-zbnds(ilayer)))
   c_mat(2*(ilayer-1)+1,2*(ilayer-1)+1) = exp(-c_ui*c_vm(ilayer)*(zbnds(ilayer+1)-zbnds(ilayer)))
   c_mat(2*(ilayer-1)+1,2*(ilayer)  )   = -c_ur   
   c_mat(2*(ilayer-1)+1,2*(ilayer)+1)   = -c_ur 
   c_mat(2*(ilayer-1)+2,2*(ilayer-1))   = exp( c_ui*c_vm(ilayer)*(zbnds(ilayer+1)-zbnds(ilayer))) * c_vm(ilayer)
   c_mat(2*(ilayer-1)+2,2*(ilayer-1)+1) = -exp(-c_ui*c_vm(ilayer)*(zbnds(ilayer+1)-zbnds(ilayer))) * c_vm(ilayer)
   c_mat(2*(ilayer-1)+2,2*(ilayer)  )   = -c_vm(ilayer+1)
   c_mat(2*(ilayer-1)+2,2*(ilayer)+1)   =  c_vm(ilayer+1)
  end if
 end do

c_fact=0.0
c_fact(2)=-c_upb_const !!! constant complex is involved in the first layer bnd as dw/dz (not w). 
call sol_comp(2*nlayer-2,c_mat,c_fact)
!c_fact=0.0
!if (ik.eq.2) c_fact(1)=1.0 !!! fixed 
!if (ik.eq.nx-1) c_fact(1)=-1.0 !!! fixed 


 do ilev=1,nlev
   ilayer = inlay(ilev)
   if (ilayer.eq.1)then !!! heating layer
    c_vw(ik,io,ilev) = c_fact(1)*sin(c_vm(1)*(zlev(ilev)-zbnds(ilayer)))
    if(ilev.le.nj) c_vw(ik,io,ilev) = c_vw(ik,io,ilev) + c_vw_hlayer(ik,io,ilev) !!! specific solution
   elseif (ilayer.eq.nlayer) then !!! radiation condition
    c_vw(ik,io,ilev) = c_fact(2*(ilayer-1)) * exp (c_ui * c_vm(ilayer) * (zlev(ilev)-zbnds(ilayer))) 
   else
    c_vw(ik,io,ilev) =  c_fact(2*(ilayer-1)) * exp (  c_ui * c_vm(ilayer) * (zlev(ilev)-zbnds(ilayer))) &
                      + c_fact(2*(ilayer-1)+1) * exp (- c_ui * c_vm(ilayer) * (zlev(ilev)-zbnds(ilayer))) 
   end if
 
  call plrz_diag(axk(ik),axo(io),ubg(ilayer),bn2(ilayer),c_vm(ilayer),c_vw(ik,io,ilev),c_vx(ik,io,ilev),c_vu(ik,io,ilev),c_vz(ik,io,ilev),vnu)
 end do

end if !!! nonzero
end do !!! io
end do !!! ik

c_vu(1:2,:,:)=0.0
c_vx(1:2,:,:)=0.0

!c_vw(1:2,:,:)=0.0
!c_vw(nx-1:nx,:,:)=0.0
!c_vw(:,1:4,:)=0.0
!c_vw(:,nt-3:nt,:)=0.0
call fft_bwd_comp_2d(nx,nt,c_vw,nlev)
call fft_bwd_comp_2d(nx,nt,c_vu,nlev)
call fft_bwd_comp_2d(nx,nt,c_vz,nlev)
call fft_bwd_comp_2d(nx,nt,c_vx,nlev)





!!!call fft_bwd_comp_2d(nx,nt,c_vw_hlayer,nj) !!! test
!!!heat(:,:,1:nj)=real(c_vw_hlayer)

w     = real(c_vw)
u     = real(c_vu)
dispz = real(c_vz)
dispx = real(c_vx)




return
end subroutine frmodel_time_2d_const

!=============================================!
!!! polarization relation
subroutine plrz_diag(vk,vo,vu,vbn2,c_vm,c_vw1,c_vx1,c_vu1,c_vz1,vnu)

real(4),intent(in)::vk,vo,vu,vbn2
complex,intent(in)::c_vm,c_vw1

complex,intent(out)::c_vx1,c_vu1,c_vz1

complex,parameter::c_ur=(1.0,0.0)
complex,parameter::c_ui=(0.0,1.0)

complex::c_oh

c_oh = vo - vk * vu + c_ui * vnu

c_vu1=0.0; c_vx1=0.0; c_vz1=0.0
if (vk.ne.0.0) c_vu1 = - c_vw1 * c_vm / vk
if (c_oh.ne.(0.0,0.0))then
 c_vx1 = c_vu1 / (- c_ui * c_oh)
 c_vz1 = c_vw1 / (- c_ui * c_oh)
end if

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

subroutine test_fftsin
use setup
real(4)::sample(nx/2-1)

!sample=sin(2.0*pi*axx(nx/2+2:nx)/(xr-xl))

!do ix=1,nx/2-1
! write(*,*) axx(nx/2+1+ix),sample(ix)
!end do

sample=0.0
sample(3)=1.0

call fft_bwd_sin_trans(1,sample,nx/2-1)

do ix=1,nx/2-1
 write(*,*) ix,sample(ix)
end do

stop
return
end subroutine test_fftsin
!=============================================!
