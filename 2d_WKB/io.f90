!============================================!
subroutine read_BG
use setup

select case(trim(cbgmode))
 case('test')
  !!! constant wind and stability 

   u0   = 10.0
   bn20 = 0.018**2

   u0_top=-3.0

 do ilev=1,nlev
   ubg(ilev)=u0 + ( u0_top - u0 ) * (zlev(ilev)-zlev(1)) / (ztop-zlev(1))
 end do

   bn2=bn20






end select

return
end subroutine read_BG

!=============================================!
subroutine read_init
use setup

  !!! sample

 xlim=5.0e4
 xwid=2.0e4
 height=4.0e2

 do ix=1,nx
  if (abs(axx(ix)).le.xlim)then
   h(ix)=height*exp(-((axx(ix))/xwid)**2)
  else
   h(ix)=0.0
  end if
 end do

return
end subroutine read_init

!=============================================!

subroutine output
use setup




return
end subroutine output

!=============================================!
