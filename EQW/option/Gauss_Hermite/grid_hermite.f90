subroutine set_grid_hermite(ny,ylen,yl,yr,grid,weight)
use func

integer,intent(in)::ny
real(4),intent(in)::ylen
real(4),intent(out)::yl,yr
real(4),intent(out)::grid(ny)
real(4),intent(out)::weight(ny)

real(4)::weight_GH(ny)

double precision::vx,vx_init
double precision::vxs_half(ny/2)
double precision::vxs_half_b(ny/2)


call set_gauss_hermite(ny,grid,weight_GH)

do iy=1,ny
 weight(iy) = weight_GH(iy) * exp(grid(iy)**2)
end do

!! length scale
grid = grid*ylen
weight  = weight*ylen

!do iy=1,ny
! write(*,*) iy,grid(iy)*1.0e-3,weight(iy)*1.0e-3
!end do

yr =  grid(ny)*1.1
yr = tdg(yr) 
yl = -yr

write(*,*) yl,yr

return
end subroutine set_grid_hermite
