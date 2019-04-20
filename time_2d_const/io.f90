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

 xlim=5.0e4
 xwid=2.0e4
 tlim=5.0e4
 twid=2.0e4
 height=4.0e2

 do ix=1,nx
 do it=1,nt
  if (abs(axx(ix)).le.xlim.and.abs(axt(it)).le.tlim)then
   h(ix,it)=height*exp(-((axx(ix))/xwid)**2-((axt(it))/twid)**2)
  else
   h(ix,it)=0.0
  end if
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
