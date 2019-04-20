!====================================!
! ISPACK FFT wrapper
!====================================!

subroutine fft_fwd_comp(ndim,carray,ndim_z)

 integer::ndim,ndim_z
 complex::carray(ndim,ndim_z)

 integer::iwork(5)
 real(8)::ywork(ndim_z,ndim,2)
 real(8)::twork(2*ndim)

 real(8)::darray(ndim_z,ndim,2)

 complex,parameter::c_ur=(1.0,0.0)
 complex,parameter::c_ui=(0.0,1.0)

 darray(:,:,1)=transpose(real(carray))
 darray(:,:,2)=transpose(aimag(carray))

 call fttzui(ndim,iwork,twork)
 call fttzuf(ndim_z,ndim,darray,ywork,iwork,twork)

 carray=transpose(darray(:,:,1))*c_ur+transpose(darray(:,:,2))*c_ui

return
end subroutine fft_fwd_comp

!====================================!

subroutine fft_bwd_comp(ndim,carray,ndim_z)

 integer::ndim,ndim_z
 complex::carray(ndim,ndim_z)

 integer::iwork(5)
 real(8)::ywork(ndim_z,ndim,2)
 real(8)::twork(2*ndim)

 real(8)::darray(ndim_z,ndim,2)

 complex,parameter::c_ur=(1.0,0.0)
 complex,parameter::c_ui=(0.0,1.0)

 darray(:,:,1)=transpose(real(carray))
 darray(:,:,2)=transpose(aimag(carray))

 call fttzui(ndim,iwork,twork)
 call fttzub(ndim_z,ndim,darray,ywork,iwork,twork)

 carray=transpose(darray(:,:,1))*c_ur+transpose(darray(:,:,2))*c_ui

return
end subroutine fft_bwd_comp

!====================================!

subroutine fft_fwd_real(ndim,rarray,ndim_z)


 integer::ndim,ndim_z
 real(4)::rarray(ndim,ndim_z)

 integer::iwork(5)
 real(8)::ywork(ndim_z,ndim)
 real(8)::twork(2*ndim)

 real(8)::darray(ndim_z,ndim)

 darray=dble(transpose(rarray))

 call fttrui(ndim,iwork,twork)
 call fttruf(ndim_z,ndim,darray,ywork,iwork,twork)

 rarray=real(transpose(darray))

return
end subroutine fft_fwd_real

!====================================!

subroutine fft_bwd_real(ndim,rarray,ndim_z)


 integer::ndim,ndim_z
 real(4)::rarray(ndim,ndim_z)

 integer::iwork(5)
 real(8)::ywork(ndim_z,ndim)
 real(8)::twork(2*ndim)

 real(8)::darray(ndim_z,ndim)

 darray=dble(transpose(rarray))

 call fttrui(ndim,iwork,twork)
 call fttrub(ndim_z,ndim,darray,ywork,iwork,twork)

 rarray=real(transpose(darray))

return
end subroutine fft_bwd_real

!====================================!
