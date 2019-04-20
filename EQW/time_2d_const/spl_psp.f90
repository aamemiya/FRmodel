!============================================!
module spl_psp
use EZspline_obj
use EZspline
use setup
type(EZspline1_r4)::v_spl
contains!------------------------------------!
subroutine spl_init

call EZspline_init(v_spl,ny,(/0,0/),ier)
v_spl%x1=axy

return
end subroutine spl_init
!------------------------------------!
subroutine spl_free

call EZspline_free(v_spl,ier)

return
end subroutine spl_free
!------------------------------------!
subroutine spl_yg(array_y,array_g)
real(4)::array_y(ny)
real(4)::array_g(nh)

call EZspline_setup(v_spl, array_y, ier)

call EZspline_interp(v_spl, nh, axy_gauss, array_g, ier)

return
end subroutine spl_yg
!------------------------------------!
end module spl_psp
!============================================!
