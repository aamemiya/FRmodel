!====================================!
! ISPACK FFT wrapper
!====================================!

subroutine fft_fwd_sin_trans(ndim_t,rarray,ndim)
 integer::ndim_x,ndim_t
 real(4)::rarray(ndim_t,ndim)  !!! ndim

 integer::iwork(5)
 real(8)::ywork(ndim_t * (ndim+1))
 real(8)::twork(5*(ndim+1)/2)
 real(8)::darray(ndim_t,ndim+1)

 darray(:,1:ndim)=dble(rarray)
 darray(:,ndim+1)=0.0

 call fttsti(ndim+1,iwork,twork)
 call fttstf(ndim_t,ndim+1,darray,ywork,iwork,twork)

 rarray=real(darray(:,1:ndim))

return
end subroutine fft_fwd_sin_trans

!====================================!

subroutine fft_bwd_sin_trans(ndim_t,rarray,ndim)
 integer::ndim_x,ndim_t
 real(4)::rarray(ndim_t,ndim)  !!! ndim

 integer::iwork(5)
 real(8)::ywork(ndim_t * (ndim+1))
 real(8)::twork(5*(ndim+1)/2)
 real(8)::darray(ndim_t,ndim+1)

 darray(:,1:ndim)=dble(rarray)
 darray(:,ndim+1)=0.0

 call fttsti(ndim+1,iwork,twork)
 call fttstb(ndim_t,ndim+1,darray,ywork,iwork,twork)

 rarray=real(darray(:,1:ndim))

return
end subroutine fft_bwd_sin_trans

!====================================!

subroutine fft_fwd_comp_2d(ndim_x,ndim_y,carray,ndim_z)

 integer::ndim_x,ndim_y,ndim_z
 complex::carray(ndim_x,ndim_y,ndim_z)

 complex::carray_temp(ndim_x,ndim_y,ndim_z)
 complex::carray_temp_t(ndim_y,ndim_x,ndim_z)

 carray_temp=carray

 call fft_fwd_comp(ndim_x,carray_temp,ndim_y*ndim_z)

do iz=1,ndim_z
 carray_temp_t(:,:,iz)=transpose(carray_temp(:,:,iz))
end do

 call fft_fwd_comp(ndim_y,carray_temp_t,ndim_x*ndim_z)

do iz=1,ndim_z
 carray(:,:,iz)=transpose(carray_temp_t(:,:,iz))
end do   

return
end subroutine fft_fwd_comp_2d

!====================================!

subroutine fft_bwd_comp_2d(ndim_x,ndim_y,carray,ndim_z)

 integer::ndim_x,ndim_y,ndim_z
 complex::carray(ndim_x,ndim_y,ndim_z)

 complex::carray_temp(ndim_x,ndim_y,ndim_z)
 complex::carray_temp_t(ndim_y,ndim_x,ndim_z)


 carray_temp=carray

 call fft_bwd_comp(ndim_x,carray_temp,ndim_y*ndim_z)

do iz=1,ndim_z
 carray_temp_t(:,:,iz)=transpose(carray_temp(:,:,iz))
end do

 call fft_bwd_comp(ndim_y,carray_temp_t,ndim_x*ndim_z)

do iz=1,ndim_z
 carray(:,:,iz)=transpose(carray_temp_t(:,:,iz))
end do   


return
end subroutine fft_bwd_comp_2d

!====================================!
subroutine fft_fwd_comp(ndim,carray,ndim_t)

 integer::ndim,ndim_t
 complex::carray(ndim,ndim_t)

 integer::iwork(5)
 real(8)::ywork(ndim_t,ndim,2)
 real(8)::twork(2*ndim)

 real(8)::darray(ndim_t,ndim,2)

 complex,parameter::c_ur=(1.0,0.0)
 complex,parameter::c_ui=(0.0,1.0)

 darray(:,:,1)=transpose(real(carray))
 darray(:,:,2)=transpose(aimag(carray))

 call fttzui(ndim,iwork,twork)
 call fttzuf(ndim_t,ndim,darray,ywork,iwork,twork)

 carray=transpose(darray(:,:,1))*c_ur+transpose(darray(:,:,2))*c_ui

return
end subroutine fft_fwd_comp

!====================================!

subroutine fft_bwd_comp(ndim,carray,ndim_t)

 integer::ndim,ndim_t
 complex::carray(ndim,ndim_t)

 integer::iwork(5)
 real(8)::ywork(ndim_t,ndim,2)
 real(8)::twork(2*ndim)

 real(8)::darray(ndim_t,ndim,2)

 complex,parameter::c_ur=(1.0,0.0)
 complex,parameter::c_ui=(0.0,1.0)

 darray(:,:,1)=dble(transpose(real(carray)))
 darray(:,:,2)=dble(transpose(aimag(carray)))

 call fttzui(ndim,iwork,twork)
 call fttzub(ndim_t,ndim,darray,ywork,iwork,twork)

 carray=transpose(real(darray(:,:,1)))*c_ur+transpose(real(darray(:,:,2)))*c_ui

return
end subroutine fft_bwd_comp

!====================================!

subroutine fft_fwd_real(ndim,rarray,ndim_t)


 integer::ndim,ndim_t
 real(4)::rarray(ndim,ndim_t)

 integer::iwork(5)
 real(8)::ywork(ndim_t,ndim)
 real(8)::twork(2*ndim)

 real(8)::darray(ndim_t,ndim)

 darray=dble(transpose(rarray))

 call fttrui(ndim,iwork,twork)
 call fttruf(ndim_t,ndim,darray,ywork,iwork,twork)

 rarray=real(transpose(darray))

return
end subroutine fft_fwd_real

!====================================!

subroutine fft_bwd_real(ndim,rarray,ndim_t)


 integer::ndim,ndim_t
 real(4)::rarray(ndim,ndim_t)

 integer::iwork(5)
 real(8)::ywork(ndim_t,ndim)
 real(8)::twork(2*ndim)

 real(8)::darray(ndim_t,ndim)

 darray=dble(transpose(rarray))

 call fttrui(ndim,iwork,twork)
 call fttrub(ndim_t,ndim,darray,ywork,iwork,twork)

 rarray=real(transpose(darray))

return
end subroutine fft_bwd_real

!====================================!
