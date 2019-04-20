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
  vbg(ilev)=0.0
  bn2(ilev)=(0.01*real(inlay(ilev)))**2
 end do

   ubg_bnd(1,1:2) = (/ 10.0,10.0 /)
   vbg_bnd(1,1:2) = (/ 0.0,0.0 /)
   bn2_bnd(1,1:2) = (/ 0.01**2,0.01**2 /)
   ubg_bnd(2,1:2) = (/ 10.0,10.0 /)
   vbg_bnd(2,1:2) = (/ 0.0,0.0 /)
   bn2_bnd(2,1:2) = (/ 0.01**2,0.02**2 /)

call diag_BG
   
   tbg_bnd(1,1:2) = tbg(1)

   do ilev=1,nlev-1
    if (inlay(ilev).eq.1.and.inlay(ilev+1).eq.2) ilev_tmp=ilev
   end do
   fac = (zbnds(2)-zlev(ilev_tmp)) / (zlev(ilev_tmp+1)-zlev(ilev_tmp))
   tbg_bnd(2,1:2) = (1.0-fac) * tbg(ilev_tmp) + fac * tbg(ilev_tmp+1)

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

integer::ih(nx,ny)

  !!! sample

 xlim=5.0e4
 xwid=2.0e4
 ylim=5.0e4
 ywid=2.0e4
 height=8.0e2

 do iy=1,ny
 do ix=1,nx
  if (abs(axx(ix)).le.xlim.and.abs(axy(iy)).le.ylim)then
   h(ix,iy)=height*exp(-((axx(ix))/xwid)**2-((axy(iy))/ywid)**2)
  else
   h(ix,iy)=0.0
  end if
 end do
 end do

return


open(11,file='/work1/amemiya/Data/TOPO/ETOPO1/ascii/dir055/319')
 do iy=1,ny/2
 do ix=1,nx/2
 read(11,*)dummy1,dummy2,ih(ix,ny/2-iy+1)
 end do
 end do
close(11)
open(11,file='/work1/amemiya/Data/TOPO/ETOPO1/ascii/dir055/320')
 do iy=1,ny/2
 do ix=1,nx/2
 read(11,*)dummy1,dummy2,ih(nx/2+ix,ny/2-iy+1)
 end do
 end do
close(11)
open(11,file='/work1/amemiya/Data/TOPO/ETOPO1/ascii/dir054/319')
 do iy=1,ny/2
 do ix=1,nx/2
 read(11,*)dummy1,dummy2,ih(ix,ny/2-iy+1+ny/2)
 end do
 end do
close(11)
open(11,file='/work1/amemiya/Data/TOPO/ETOPO1/ascii/dir054/320')
 do iy=1,ny/2
 do ix=1,nx/2
 read(11,*)dummy1,dummy2,ih(nx/2+ix,ny/2-iy+1+ny/2)
 end do
 end do
close(11)

h=real(ih)
return

where(ih.lt.0) ih=0
have=sum(ih)/nx/ny

ih=ih-have

ih(:,nx-9:nx)=0
ih(:,1:10)=0
ih(ny-9:ny,:)=0
ih(1:10,:)=0

!!! smth

ismth=5
h=0.0
do ix=1+ismth,nx-ismth
 h(ix,:)= real(sum(ih(ix-ismth:ix+ismth,:),1))/real(2*ismth+1)
end do
ih=int(h)
do iy=1+ismth,ny-ismth
 h(:,iy)= real(sum(h(:,iy-ismth:iy+ismth),2))/real(2*ismth+1)
end do


!!!


return
end subroutine read_init

!=============================================!

subroutine output
use setup




return
end subroutine output

!=============================================!
