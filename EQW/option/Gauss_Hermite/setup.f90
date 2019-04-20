!==============================================================!
module param

integer,parameter::ny=32
integer,parameter::nh=15

!real(4),parameter::yl=-5000.0e3,yr=5000.0e3
!real(4),parameter::yl=-12.0,yr=12.0
real(4),parameter::ylen=2.0e6

real(4)::yl,yr
real(4)::axy(ny)
real(4)::wty(ny)
end module param
!====================================!

module func
contains !===========================!
real function hermite_r(vy,n)
 real(4)::vy
 integer::n

 if (n.eq.0)then
  hermite_r=1.0
 elseif (n.eq.1)then
  hermite_r=2.0*vy
 else
  bb=1.0
  b=2.0*vy
  do i=2,n
   a=2.0*vy*b-2.0*real(i-1)*bb
   bb=b
   b=a
  end do
  hermite_r=a
 end if

! select case(n)
! case(0); hermite_r=1.0
! case(1); hermite_r=2.0*vy
! case(2); hermite_r=4.0*vy**2-2.0
! case(3); hermite_r=8.0*vy**3-16.0*vy
! case(4); hermite_r=16.0*vy**4-64.0*vy**2+16.0
! end select
end function hermite_r
!
double precision function hermite_d(vy,n)
 double precision::vy
 double precision::a,b,bb
 integer::n

 if (n.eq.0)then
  hermite_d=1.0
 elseif (n.eq.1)then
  hermite_d=2.0*vy
 else
  bb=1.0
  b=2.0*vy
  do i=2,n
   a=2.0*vy*b-2.0*real(i-1)*bb
   bb=b
   b=a
  end do
  hermite_d=a
 end if

end function hermite_d
!
real function para_cyl(vy,n)
 real(4)::vy
 integer::n
if (n.eq.-1)then !!! Kelvin wave
 para_cyl=exp(-0.5*vy**2)
elseif (n.ge.0)then
 para_cyl=sqrt(1.0/((2**n)*sqrt(3.1416)*rfactorial(n)))*exp(-0.5*vy**2)*hermite_r(vy,n)
else
 write(*,*)'ERROR: function para_cyl'
 write(*,*)'n should be n.ge.-1'
 stop
end if


end function para_cyl
!
integer function factorial(j)
integer j
  factorial=1
 if (j.ge.2) then
  do jl=2,j
   factorial=factorial*jl
  end do
 end if
end function factorial
!
real function rfactorial(j)
integer j
  rfactorial=1.0
 if (j.ge.2) then
  do jl=2,j
   rfactorial=rfactorial*real(jl)
  end do
 end if
end function rfactorial
!
real function log_rfactorial(j)
integer j
  log_rfactorial=0.0
 if (j.ge.2) then
  do jl=2,j
   log_rfactorial=log_rfactorial+log(real(jl))
  end do
 end if
end function log_rfactorial
!
real function tdg(val) !!! 3.1415 ==> 3.1
 real(4):: val
 idig_val=int(log(val)/log(10.0))
 tdg=(int(val/10**(idig_val-1))+1)*10**(idig_val-1) 

end function tdg
end module func
!====================================!
module output
use param

real(4)::f(ny)
real(4)::f_new(ny)
real(4)::f_new_eqw(ny)
real(4)::f_new_comp(ny,nh)
real(4)::f_new_comp_eqw(ny,nh)

real(4)::fh(nh)

character*20::title1,title2
character*20::xlabel,ylabel
real(4)::vmin,vmax
real(4)::astics,amtics
real(4)::bstics,bmtics
integer::iout

end module output
!==============================================================!
