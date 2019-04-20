!=============================================!
subroutine decomp_test
use setup
use func
real(4)::heat_j(nx,nt,ny,nj)
complex::c_heat_ko_j(nx,nt,ny,nj)

complex::c_heat_ko_jh(nx,nt,nh,nj)

real(4)::rwork_nh(nh)
real(4)::rwork_ny(ny)

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
 write(*,*) ij
 vce2=bn2(1)*((zbnds(2)-zbnds(1))/(real(ij)*pi))**2
 ylen_cebr=sqrt(sqrt(vce2)/beta)
 do ik=1,nx
! do io=1,nt
 do io=1,1
  rwork_ny=real(c_heat_ko_j(ik,io,:,ij))
  call decomp_bwd(rwork_ny,rwork_nh,ylen_cebr)
  c_heat_ko_jh(ik,io,:,ij)=c_ur*rwork_nh
  rwork_ny=aimag(c_heat_ko_j(ik,io,:,ij))
  call decomp_bwd(rwork_ny,rwork_nh,ylen_cebr)
  c_heat_ko_jh(ik,io,:,ij)=c_heat_ko_jh(ik,io,:,ij)+c_ui*rwork_nh

! if (ij.eq.2)then
!   write(*,*) 'aliasing test'
!  do ih=1,ny
!   write(*,'(2I6,3E12.3)')ih, int(ylen_cebr*1.0e-3)/ih, (/(para_cyl(axy(iy)/yle!n_cebr,ih) ,iy=ny/2-1,ny/2+1 )/)
!  end do 
!  stop

! end if
end do
end do
end do


!c_heat_ko_jh(:,:,11:nh,:)=0.0





do ik=1,nx
!do io=1,1
do io=1,1

do ij=1,nj
 vce2=bn2(1)*((zbnds(2)-zbnds(1))/(real(ij)*pi))**2
 ylen_cebr=sqrt(sqrt(vce2)/beta)
 do ih=1,nh
  rwork_nh=real(c_heat_ko_jh(ik,io,:,ij))
  call decomp_fwd(rwork_nh,rwork_ny,ylen_cebr)
  c_heat_ko_j(ik,io,:,ij)=c_ur*rwork_ny
  rwork_nh=aimag(c_heat_ko_jh(ik,io,:,ij))
  call decomp_fwd(rwork_nh,rwork_ny,ylen_cebr)
  c_heat_ko_j(ik,io,:,ij)=c_heat_ko_j(ik,io,:,ij)+ c_ui*rwork_ny
 end do
end do

end do
end do


 call fft_fwd_comp_2d(nx,nt,c_heat_ko_j,ny*nj)

 heat_j=real(c_heat_ko_j)

 call fft_fwd_sin_trans(nx*nt*ny,heat_j,nj)

heat_diag(:,:,:,1:nj)=heat_j


 
return
end subroutine decomp_test
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
