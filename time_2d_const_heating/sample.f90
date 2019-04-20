program sample

real(4)::a
complex::c_a

complex::c_f
complex::c_g

complex,parameter::c_ur=(1.0,0.0)
complex,parameter::c_ui=(0.0,1.0)

a = 3.1415 * 0.5

c_a= a * c_ur + a * c_ui 

c_f=sin(a)

write(*,*) c_f

c_f=sin(c_a)

c_g = ( exp(c_ui * c_a) - exp (- c_ui * c_a) ) * 0.5 / c_ui

write(*,*) c_f,c_g

stop

end program sample
