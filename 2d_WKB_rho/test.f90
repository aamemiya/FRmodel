program test

complex::c_test
complex,parameter::c_ur=(1.0,0.0)
complex,parameter::c_ui=(0.0,1.0)

a=-1.0

c_test=sqrt(a*c_ui)

write(*,*) c_test

c_test=(a*c_ui)**0.5

write(*,*) c_test

stop
end program test
