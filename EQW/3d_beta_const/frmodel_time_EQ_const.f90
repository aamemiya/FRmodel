!=============================================!

subroutine frmodel_time_eq_const
use setup
use func

real(4)::heat_j(nx,nt,ny,nj)
complex::c_heat_ko_j(nx,nt,ny,nj)

complex::c_heat_ko_jh(nx,nt,nh,nj)

complex::cwork_nh(nh)
real(4)::rwork_nh(nh)
real(4)::rwork_ny(ny)



complex::c_vw_hlayer(nx,nt,nh,nj)

complex::c_nj_work(nj)
real(4)::r_2_nj_work(2,nj)
complex::c_upb_const


complex::vsgn
complex::c_vm(nlayer)

complex::c_vw(nx,nt,ny,nlev)
complex::c_vu(nx,nt,ny,nlev)
complex::c_vz(nx,nt,ny,nlev)
complex::c_vx(nx,nt,ny,nlev)

complex::c_fact(2*nlayer-2)
complex::c_mat(2*nlayer-2,2*nlayer-2)
complex::c_mat_t(2*nlayer-2,2*nlayer-2)


complex::fc_vm

real(4),parameter::beta=2.0*eomg/er
complex,parameter::c_ur=(1.0,0.0)
complex,parameter::c_ui=(0.0,1.0)

!!! heat(x,t,y,z) --> heat_j(x,t,y,j)

heat_j(:,:,:,:)=heat(:,:,:,1:nj)

 call fft_bwd_sin_trans(nx*nt*ny,heat_j,nj)

!!! heat_j(x,t,y,j) --> c_heat_ko_j(k,omega,y,j)

c_heat_ko_j = heat_j * c_ur

 call fft_bwd_comp_2d(nx,nt,c_heat_ko_j,ny*nj)
do ij=1,nj
  write(*,*) 'decomp',ij,' /',nj
  vce2=bn2(1)*((zbnds(2)-zbnds(1))/(real(ij)*pi))**2
  ylen_cebr=sqrt(sqrt(vce2)/beta)
 do ik=1,nx
 do io=1,nt
  rwork_ny=real(c_heat_ko_j(ik,io,:,ij))
  call decomp_bwd(rwork_ny,rwork_nh,ylen_cebr)
  c_heat_ko_jh(ik,io,:,ij)=c_ur*rwork_nh
  rwork_ny=aimag(c_heat_ko_j(ik,io,:,ij))
  call decomp_bwd(rwork_ny,rwork_nh,ylen_cebr)
  c_heat_ko_jh(ik,io,:,ij)=c_heat_ko_jh(ik,io,:,ij)+c_ui*rwork_nh
 end do
 end do
end do

!!! for each set of k, omega and ih ...

do ik=1,nx
do io=1,nt

do ih=1,nh

if (axo(io).ne.0.0.and.axk(ik).ne.0.0)then

 do ilayer=1,nlayer
  c_vm(ilayer)=fc_vm(bn2(ilayer),ubg(ilayer),axk(ik),axo(io),ih-1,beta,vnu) !!! WN ih-1
 end do
! write(*,*) 'k,omega,ih,c_vm',axk(ik),axo(io),ih,c_vm(1)

 !!! specific solution in the heating layer
 c_upb_const = 0.0
 do ij=1,nj !!! vertical WN space
 vce2=bn2(1)*((zbnds(2)-zbnds(1))/(real(ij)*pi))**2
  if(axo(io).ne.0.0) c_nj_work(ij) = c_heat_ko_jh(ik,io,ih,ij)/(c_vm(1)**2-(pi*real(ij)/(zbnds(2)-zbnds(1)))**2) / vce2
  c_upb_const=c_upb_const + real(ij)*pi/(zbnds(2)-zbnds(1)) * c_nj_work(ij) * real(1-2*mod(ij,2))
 end do
 r_2_nj_work(1,:)=real(c_nj_work(:))
 r_2_nj_work(2,:)=aimag(c_nj_work(:))

 call fft_bwd_sin_trans(2,r_2_nj_work,nj)
 c_vw_hlayer(ik,io,ih,:) = c_ur * r_2_nj_work(1,:) + c_ui * r_2_nj_work(2,:)

!!! heat(ik,io,:)=c_vw_hlayer(ik,io,:) !!! test


 !!! coefficient matrix
  c_mat=0.0
 do ilayer=1,nlayer-1
  if (ilayer.eq.1)then 
!   write(*,*) c_vm(ilayer), (zbnds(2)-zbnds(1)), sin(real(c_vm(ilayer)) * (zbnds(2)-zbnds(1)))
   c_mat(2*(ilayer-1)+1,1) = sin(real(c_vm(ilayer)) * (zbnds(2)-zbnds(1)))
   c_mat(2*(ilayer-1)+1,2*(ilayer)  )   = -c_ur   
   if(nlayer.gt.2) c_mat(2*(ilayer-1)+1,2*(ilayer)+1)   = -c_ur   
   c_mat(2*(ilayer-1)+2,1) = c_vm(ilayer)*cos(real(c_vm(ilayer)) * (zbnds(2)-zbnds(1))) 
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
    c_vw(ik,io,ih,ilev) = c_fact(1)*sin(real(c_vm(1))*(zlev(ilev)-zbnds(ilayer)))
    if(ilev.le.nj) c_vw(ik,io,ih,ilev) = c_vw(ik,io,ih,ilev) + c_vw_hlayer(ik,io,ih,ilev) !!! specific solution
   elseif (ilayer.eq.nlayer) then !!! radiation condition
    c_vw(ik,io,ih,ilev) = c_fact(2*(ilayer-1)) * exp (c_ui * c_vm(ilayer) * (zlev(ilev)-zbnds(ilayer))) 
   else
    c_vw(ik,io,ih,ilev) =  c_fact(2*(ilayer-1)) * exp (  c_ui * c_vm(ilayer) * (zlev(ilev)-zbnds(ilayer))) &
                         + c_fact(2*(ilayer-1)+1) * exp (- c_ui * c_vm(ilayer) * (zlev(ilev)-zbnds(ilayer))) 
   end if
 
!  call plrz_diag(axk(ik),axo(io),ubg(ilayer),bn2(ilayer),c_vm(ilayer),c_vw(ik,io,ilev),c_vx(ik,io,ilev),c_vu(ik,io,ilev),c_vz(ik,io,ilev),vnu)
 end do

end if !!! nonzero
end do !!! ih
end do !!! io
end do !!! ik

 write(*,*) 'decomp'
do ik=1,nx
 write(*,*) ik,' / ',nx
do io=1,nt
do ilev=1,nlev
 ilayer = inlay(ilev)
 vce2=bn2(1)*((zbnds(2)-zbnds(1))/(real(ij)*pi))**2
 ylen_cebr=sqrt(sqrt(vce2)/beta)
! write(*,*) 'k,omega,c_vm,ylen',axk(ik),axo(io),c_vm(1),ylen_cebr
  cwork_nh=c_vw(ik,io,1:nh,ilev)
  rwork_nh=real(cwork_nh)
  call decomp_fwd(rwork_nh,rwork_ny,ylen_cebr)
  c_vw(ik,io,:,ilev)=c_ur*rwork_ny
  rwork_nh=aimag(cwork_nh)
  call decomp_fwd(rwork_nh,rwork_ny,ylen_cebr)
  c_vw(ik,io,:,ilev)=c_vw(ik,io,:,ilev)+ c_ui*rwork_ny
end do
end do
end do

call fft_bwd_comp_2d(nx,nt,c_vw,nlev)
!call fft_bwd_comp_2d(nx,nt,c_vu,nlev)
!call fft_bwd_comp_2d(nx,nt,c_vz,nlev)
!call fft_bwd_comp_2d(nx,nt,c_vx,nlev)

w     = real(c_vw)
u     = real(c_vu)
dispz = real(c_vz)
dispx = real(c_vx)




return
end subroutine frmodel_time_eq_const

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
complex function fc_vm(vbn2,vu,vk,vo,ih,beta,vnu)

 real(4)::vbn2,vu,vk
 integer::ih
 complex,parameter::c_ur=(1.0,0.0)
 complex,parameter::c_ui=(0.0,1.0)
 complex::c_oh

 complex::c_vce,c_kwb,c_kwbw

 c_oh   = vo - vk * vu + c_ui * vnu

 c_kwb  = c_oh * vk / beta
 c_kwbw = c_kwb**2 + c_kwb 

if (ih.eq.-1) then !!! Kelvin wave
 c_vce = c_oh / vk

elseif (ih.eq.0)then
 if (real(c_kwb).gt.-1.0)then !!! Rossby-Gravity wave
  c_vce =  c_oh**2/beta / (1+c_kwb)
 end if

else
 if (real(c_kwbw).ge.0.0)then !!! Inertia-Gravity wave
  c_vce =  c_oh**2/beta / ((real(ih)+0.5)+((real(ih)+0.5)**2+c_kwbw)**0.5)
 else !!! EQ-Rossby wave
  c_vce =  c_oh**2/beta / ((real(ih)+0.5)-((real(ih)+0.5)**2+c_kwbw)**0.5)
 end if

  
end if

 fc_vm = vbn2**0.5 / c_vce
 if (real(fc_vm).ne.0.0)then
 !!! propagate
  fc_vm = fc_vm *  ( sign(1.0,real(-c_oh)))
 elseif (aimag(fc_vm).ne.0.0)then
 !!! evanescent
  fc_vm = fc_vm * ( aimag(fc_vm)/abs(aimag(fc_vm)) )
 end if

! if (aimag(fc_vm).ne.0.0)then
! !!! evanescent
!  fc_vm = fc_vm * ( aimag(fc_vm)/abs(aimag(fc_vm)) )
! elseif (real(fc_vm).ne.0.0)then
! !!! propagate
!  fc_vm = fc_vm *  ( sign(1.0,real(-c_oh)))
! end if


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



subroutine decomp_bwd(arrayin,arrayout,ylen)
use setup
use func 

real(4)::arrayin(ny)
real(4)::arrayout(nh)
real(4)::ylen


do ih=0,nh-1
 arrayout(ih+1)=0.0
 do iy=1,ny
  arrayout(ih+1) = arrayout(ih+1) + (yr-yl)/real(ny) * para_cyl(axy(iy)/ylen,ih)*arrayin(iy) / ylen
 end do
end do

return
end subroutine decomp_bwd
!=====================================================!
subroutine decomp_fwd(arrayin,arrayout,ylen)
use setup
use func 

real(4)::arrayin(nh)
real(4)::arrayout(ny)
real(4)::ylen

do iy=1,ny
 arrayout(iy)=0.0
 do ih=0,nh-1
  arrayout(iy) = arrayout(iy) + para_cyl(axy(iy)/ylen,ih)*arrayin(ih+1)
 end do
end do

return
end subroutine decomp_fwd
!=====================================================!
