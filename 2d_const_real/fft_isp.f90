!====================================!
! ISPACK FFT wrapper
!====================================!

subroutine fft_fwd_comp(ndim,carray,ndim_z)

 integer::ndim,ndim_z
 complex::carray(ndim,ndim_z)



return
end subroutine fft_fwd_comp

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
