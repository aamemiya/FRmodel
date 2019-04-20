!============================================!
subroutine read_BG
use setup

select case(trim(cbgmode))
 case('test')
  !!! constant wind and stability 

   u0   = 10.0
   bn20 = 0.018**2

   ubg=u0
   bn2=bn20
end select

return
end subroutine read_BG

!=============================================!
subroutine read_init
use setup

  !!! sample

 xlim=2.0e6
 xwid=1.0e6
 ylim=2.0e6
 ywid=1.0e6
! tlim=1.0e6
! twid=5.0e5
 tlim=9.9e10
 twid=9.9e10
 height=4.0e2

 do ix=1,nx
 do it=1,nt
 do iy=1,ny
  if (abs(axx(ix)).le.xlim.and.abs(axt(it)).le.tlim.and.abs(axy(iy)).le.ylim)then
   h(ix,it,iy)=height*exp(-((axx(ix))/xwid)**2-((axt(it))/twid)**2-((axy(iy))/ywid)**2)
  else
   h(ix,it,iy)=0.0
  end if
 end do
 end do
 end do

return
end subroutine read_init

!=============================================!

subroutine output
use setup




return
end subroutine output

!=============================================!
