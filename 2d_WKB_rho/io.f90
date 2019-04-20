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

call diag_BG


end select

return
end subroutine read_BG

!=============================================!
subroutine diag_BG
use setup

real(4),parameter::t00=300.0
real(4)::tmp(nlev)

pbg(1)=ps0
tbg(1)=t00
thbg(1)=t00
rhobg(1)=pbg(1)/rd/tbg(1)


tmp(1)=log(thbg(1))
do ilev=2,nlev
 tmp(ilev)=tmp(ilev-1)+ bn2(ilev) * (zlev(ilev)-zlev(ilev-1)) / grav
end do

thbg=exp(tmp)

tmp(1)=1.0
do ilev=2,nlev
 tmp(ilev)=tmp(ilev-1)- grav/cp/thbg(ilev) * (zlev(ilev)-zlev(ilev-1))
end do

pbg=ps0*tmp**(cp/rd)
tbg=thbg*(pbg/ps0)**(rd/cp)
rhobg=pbg/rd/tbg

return
end subroutine diag_BG
!=============================================!

subroutine read_init
use setup

  !!! sample

 xlim=5.0e4
 xwid=2.0e4
 height=2.0e2

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
