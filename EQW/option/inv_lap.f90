!!!========================================!!!
!!! LAPACK inverse subroutine wrapper
!!!========================================!!!

subroutine inv_comp(ndim,carray)

integer,intent(in)::ndim
complex,intent(inout)::carray(ndim,ndim)

integer::ipiv(ndim)

!!integer,parameter::nwork_64=320 !!! optimal
integer,parameter::nwork_64=192 !!! optimal
!!integer,parameter::nwork_64=1 !!! optimal
complex::cwork(nwork_64)


!nwork_test=-1

ipiv=0

!!write(*,*) 'check'

!!write(*,*) ndim,carray,ndim,ipiv,cwork,nwork_64,info


call cgetri(ndim,carray,ndim,ipiv,cwork,nwork_64,info)

!!if(info.ne.0) write(*,*) 'inv_comp info',info

!!write(*,*) 'optimal nwork', cwork
!!stop

return
end subroutine inv_comp
!!!========================================!!!

subroutine sol_comp(ndim,carray,cvec)

integer,intent(in)::ndim
complex,intent(inout)::carray(ndim,ndim)
complex,intent(inout)::cvec(ndim)
complex::cvec_1(ndim,1)

integer::ipiv(ndim)

cvec_1(:,1)=cvec

call cgesv(ndim,1,carray,ndim,ipiv,cvec_1,ndim,info)

cvec=cvec_1(:,1)

return
end subroutine sol_comp

!!!========================================!!!

subroutine sol_real(ndim,rarray,rvec)

integer,intent(in)::ndim
real(4),intent(inout)::rarray(ndim,ndim)
real(4),intent(inout)::rvec(ndim)
real(4)::rvec_1(ndim,1)

integer::ipiv(ndim)

rvec_1(:,1)=rvec

call sgesv(ndim,1,rarray,ndim,ipiv,rvec_1,ndim,info)

rvec=rvec_1(:,1)

return
end subroutine sol_real

!!!========================================!!!
