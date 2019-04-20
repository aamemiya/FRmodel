!!
!! Linear algebla with LAPACK libraries
!! 
!! comparison between solution of linear equations
!! obtained through inverse matrix and direct solver
!!
!! inverse matrix method is not useful and currently problematic
!! while directo solver (sol_comp) is simple and correct
!!
program sample

complex::ca(3,3)
complex::ch(3)
complex::cx(3)

ch(1:3)=(/1.0,0.0,0.0/)

ca(1:3,1)=(/1.0,1.0,0.0/)
ca(1:3,3)=(/1.0,0.0,-1.0/)
ca(1:3,2)=(/0.0,1.0,0.0/)

call inv_comp(3,ca)

do i=1,3
 write(*,*) ca(1:3,i)
end do

cx(1:3)=ca(1,1:3)*ch(1:3)

write(*,*) cx



ca(1,1:3)=(/1.0,1.0,0.0/)
ca(3,1:3)=(/1.0,0.0,-1.0/)
ca(2,1:3)=(/0.0,1.0,0.0/)
cx=ch
call sol_comp(3,ca,cx)

write(*,*) cx

stop
end program sample
