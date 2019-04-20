!============================================!
subroutine read_BG
use setup

select case(trim(cbgmode))
 case('test')
  !!! constant wind and stability 

!   u0   = 10.0
!   bn20 = 0.018**2

!   ubg=u0
!   bn2=bn20

!   ubg = (/ 10.0,10.0,10.0 /)
!   bn2 = (/ 0.01**2,0.01**2,0.01**2 /)
   ubg = (/ 10.0,10.0 /)
   bn2 = (/ 0.01**2,0.02**2 /)

end select

!!! layer #

do ilev=1,nlev
 do ilayer=1,nlayer
  if(zlev(ilev).ge.zbnds(ilayer).and.zlev(ilev).lt.zbnds(ilayer+1)) inlay(ilev)=ilayer  
 end do
end do

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
