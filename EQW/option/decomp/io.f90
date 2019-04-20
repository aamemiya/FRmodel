!============================================!
subroutine read_BG
use setup

select case(trim(cbgmode))
 case('test')
  !!! constant wind and stability 

   u0   = 0.0
   bn20 = 0.018**2

   ubg=u0
   bn2=bn20
!   ubg = (/ 0.0,0.0 /)
!   bn2 = (/ 0.01**2,0.01**2 /)

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
 
 xlim=1.0e6
 xwid=5.0e5
 ylim=5.0e5
 ywid=2.0e5
 tlim=8.0e4
 twid=1.0e4

 zlim=3.0e3
 zwid=2.0e3
 zcen=5.0e3


 dKdt=10.0 !!! K/day
 t00=300.0

 heat_amp = dKdt * grav/t00 / 86400.0

 do ix=1,nx
 do it=1,nt
 do iy=1,ny
 do ilev=1,nlev

  if (abs(axx(ix)).le.xlim.and.abs(axy(iy)).le.ylim.and. &
      abs(axt(it)).le.tlim.and.abs(zlev(ilev)-zcen).le.zlim.and.inlay(ilev).eq.1)then
   heat(ix,it,iy,ilev)=heat_amp*exp(-((axx(ix))/xwid)**2-((axy(iy))/ywid)**2-((axt(it))/twid)**2)*sin(pi*zlev(ilev)/(zbnds(2)-zbnds(1)))
!   heat(ix,it,iy,ilev)=heat_amp*sin(pi*2.0*axt(it)/twid)*exp(-((axx(ix))/xwid)**2-((axt(it))/twid)**2)*sin(pi*zlev(ilev)/(zbnds(2)-zbnds(1)))
  else
   heat(ix,it,iy,ilev)=0.0
  end if
 end do
 end do
 end do
 end do


heat(:,1,:,:)=sum(heat,2)/real(nt)

do it=2,nt
 heat(:,it,:,:)=heat(:,1,:,:)
end do

return
end subroutine read_init

!=============================================!

subroutine output
use setup




return
end subroutine output

!=============================================!
