!============================================!
subroutine read_BG
use setup

!!! layer #

do ilev=1,nlev
 do ilayer=1,nlayer
  if(zlev(ilev).ge.zbnds(ilayer).and.zlev(ilev).lt.zbnds(ilayer+1)) inlay(ilev)=ilayer  
 end do
end do

select case(trim(cbgmode))
 case('test')
  !!! constant wind and stability 

!   u0   = 10.0
!   bn20 = 0.018**2

!   ubg=u0
!   bn2=bn20

!   ubg = (/ 10.0,10.0,10.0 /)
!   bn2 = (/ 0.01**2,0.01**2,0.01**2 /)
!   ubg = (/ 10.0,10.0 /)
!   bn2 = (/ 0.01**2,0.02**2 /)

 do ilev=1,nlev
!  ubg(ilev)=10.0+(2-real(inlay(ilev)))*(zbnds(2)-zlev(ilev))/(zbnds(2)-zbnds(1))*10.0
  ubg(ilev)=10.0+(real(inlay(ilev))-1)*(zlev(ilev)-zbnds(2))/(zbnds(2)-zbnds(1))*10.0
  bn2(ilev)=(0.01*real(inlay(ilev)))**2
 end do

   ubg_bnd(1,1:2) = (/ 10.0,10.0 /)
   bn2_bnd(1,1:2) = (/ 0.01**2,0.01**2 /)
   ubg_bnd(2,1:2) = (/ 10.0,10.0 /)
   bn2_bnd(2,1:2) = (/ 0.01**2,0.02**2 /)





end select




return
end subroutine read_BG

!=============================================!
subroutine read_init
use setup

  !!! sample

 xlim=5.0e4
 xwid=2.0e4
 height=8.0e2

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
