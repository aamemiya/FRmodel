!============================================!
subroutine read_BG
use setup

select case(trim(cbgmode))
 case('test')

 zbnds=(/ zbot,1.5e4,ztop /)

!!! layer #

do ilev=1,nlev
 do ilayer=1,nlayer
  if(zlev(ilev).ge.zbnds(ilayer).and.zlev(ilev).lt.zbnds(ilayer+1)) inlay(ilev)=ilayer  
 end do
end do

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
!   tbg_bnd(2,1:2) = (1.0-fac) * tbg(ilev_tmp) + fac * tbg(ilev_tmp+1)


    vudir=0.0 !!! eastward wind

case ('file')
    call profile
    call diag_BG
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
subroutine diag_BG
use setup

real(4),parameter::t00=300.0
real(4)::tmp(nlev)

pbg(1)=ps0
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
rhobg=pbg/rd/tbg

if (trim(cbgmode).eq.'test') tbg=thbg*(pbg/ps0)**(rd/cp)

return
end subroutine diag_BG
!=============================================!
subroutine read_init
use setup

integer::ih(nx,ny)


select case (trim(chmode))
 case('test')

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

 case('etopo1')

 call etopo1_init

end select


return
end subroutine read_init

!=============================================!

subroutine output
use setup




return
end subroutine output

!=============================================!
!===================================================================================!
!!!! procedure !!!!!!

!!! 1. calculate grid points vlon_grid(ix,iy),vlat_grid(ix,iy)
!!!    considering the spherisity of the earth
!!! 2. determine necessary original data range (longitude,latitude)
!!! 3. read ETOPO1 data (1/60 deg.) 
!!! 4. interpolate h to vlon_grid,vlat_grid
!!! 5. let data points below the sea surface be zero altitude 
!!! 6. convolute h with window function to smooth out small scales
!!! 7. remove data near the boundaries using proper cutoff function 
!===================================================================================!
subroutine etopo1_init
use setup
use EZspline_obj
use EZspline

namelist /etopo1_init_nml/ vlonc,vlatc,vudir_fix,hfac

real(4)::vlon_grid(nx_data,ny_data),vlat_grid(nx_data,ny_data)
real(4)::rx(nx_data,ny_data),ry(nx_data,ny_data),rz(nx_data,ny_data)
real(4)::rxx(nx_data,ny_data),ryy(nx_data,ny_data),rzz(nx_data,ny_data),rhh(nx_data,ny_data)
real(4)::h_xy(nx_data,ny_data)
real(4),allocatable::h_lonlat(:,:),axlon(:),axlat(:)
real(4),allocatable::vx_isolon(:,:),vy_isolon(:,:),vx_isolat(:,:),vy_isolat(:,:)

real(4)::h_copy(nx,ny)

type(EZspline2_r4)::h_spl


vudir_fix=-999.  !!! default not fixed

open(11,file='calc.conf')
read(11,nml=etopo1_init_nml)
close(11)

if (vudir_fix.ge.-180.0.and.vudir_fix.le.360.0) udir=vudir_fix


!!! 1. calculate grid points !!!

!!! initial location (rad)
do ix=1,nx_data
 vlon_grid(ix,:) = -0.5*vxlen/er + (real(ix)-0.5) / real(nx_data) * vxlen/er
end do
do iy=1,ny_data
 vlat_grid(:,iy) = -0.5*vylen/er + (real(iy)-0.5) / real(ny_data) * vylen/er
end do

 write(*,*) 'before'
 write(*,*) vlon_grid(1,1)  /drad,   vlat_grid(1,1)/drad
 write(*,*) vlon_grid(nx_data,1) /drad,  vlat_grid(nx_data,1)/drad
 write(*,*) vlon_grid(1,ny_data) /drad,  vlat_grid(1,ny_data)/drad
 write(*,*) vlon_grid(nx_data,ny_data)/drad, vlat_grid(nx_data,ny_data)/drad


rx=cos(vlat_grid)*cos(vlon_grid)
ry=cos(vlat_grid)*sin(vlon_grid)
rz=sin(vlat_grid)


cxy=cos(drad*vlonc); sxy=sin(drad*vlonc)
cxz=cos(drad*vlatc); sxz=sin(drad*vlatc)
cyz=cos(drad*vudir); syz=sin(drad*vudir)

rxx = cxy*cxz * rx + ( -cxy*sxz*syz - sxy*cyz ) * ry + ( -cxy*sxz*cyz + sxy*syz ) * rz
ryy = sxy*cxz * rx + ( -sxy*sxz*syz + cxy*cyz ) * ry + ( -sxy*sxz*cyz - cxy*syz ) * rz 
rzz=      sxz * rx +                cxz*syz * ry +                    cxz*cyz * rz

rhh=sqrt(rxx**2+ryy**2)

do iy=1,ny_data
do ix=1,nx_data
 call vec_to_deg(rxx(ix,iy),ryy(ix,iy),vlon_grid(ix,iy))
 call vec_to_deg(rhh(ix,iy),rzz(ix,iy),vlat_grid(ix,iy))
end do
end do

 write(*,*) 'after'
 write(*,*) rxx(1,1)  , rxx(1,1)  ,  rzz(1,1),  rhh(1,1)
 write(*,*) rxx(nx_data,1) , ryy(nx_data,1) , rzz(nx_data,1)
 write(*,*) rxx(1,ny_data) , ryy(1,ny_data) , rzz(1,ny_data)
 write(*,*) rxx(nx_data,ny_data), ryy(nx_data,ny_data), rzz(nx_data,ny_data)

 write(*,*) 
 write(*,*) 'after'
 write(*,*) vlon_grid(1,1)  ,   vlat_grid(1,1)
 write(*,*) vlon_grid(nx_data,1) ,  vlat_grid(nx_data,1)
 write(*,*) vlon_grid(1,ny_data) ,  vlat_grid(1,ny_data)
 write(*,*) vlon_grid(nx_data,ny_data), vlat_grid(nx_data,ny_data)



!!! 2. Init data range !!!


ilonl=int(minval(vlon_grid))
ilonr=int(maxval(vlon_grid))+1
ilatl=int(minval(vlat_grid))
ilatr=int(maxval(vlat_grid))+1

nlon=(ilonr-ilonl)*60
nlat=(ilatr-ilatl)*60


allocate(h_lonlat(nlon,nlat),axlon(nlon),axlat(nlat))

axlon=(/( real(ilonl)+(real(ilon)-0.5)/60.0 , ilon=1,nlon)/)
axlat=(/( real(ilatl)+(real(ilat)-0.5)/60.0 , ilat=1,nlat)/)

!!! 3. Read ETOPO1 data !!!

do ilat=ilatl,ilatr-1
do ilon=ilonl,ilonr-1
 ilons=(ilon-ilonl)*60+1
 ilone=(ilon-ilonl)*60+60
 ilats=(ilat-ilatl)*60+1
 ilate=(ilat-ilatl)*60+60
 call load(ilon,ilat,h_lonlat(ilons:ilone,ilats:ilate))
end do
end do

!!! 4. interpolate h onto vlon_grid,vlat_grid

call EZspline_init(h_spl,nlon,nlat,(/0,0/),(/0,0/),ier)
h_spl%x1=axlon
h_spl%x2=axlat
call EZspline_setup(h_spl, h_lonlat, ier)
do iy=1,ny_data
 call EZspline_interp(h_spl,nx_data,vlon_grid(:,iy),vlat_grid(:,iy),h_xy(:,iy),ier) !!! treat as 1-d array
end do
call EZspline_free(h_spl,ier)


!!! 5. let data points below the sea surface be zero altitude 

where(h_xy.lt.0.0) h_xy=0.0

write(*,*) maxval(h_lonlat),minval(h_lonlat)
write(*,*) maxval(h_xy),minval(h_xy)

!!! 6. convolute h with window function to smooth out small scales


!!! 7. remove data near the boundaries using proper cutoff function 

h=0.0

h(nx/2-nx_data/2+1:nx/2+nx_data/2,ny/2-ny_data/2+1:ny/2+ny_data/2) &
 = h_xy


h_copy=h !!! TORI AEZU

!!! smoothing TORI AEZU
ismth=1
!ismth=0
do ix=ismth+1,nx-ismth
 h(ix,:)=sum(h_copy(ix-ismth:ix+ismth,:),1)/real(2*ismth+1)
end do
h_copy=h
do iy=ismth+1,ny-ismth
 h(:,iy)=sum(h_copy(:,iy-ismth:iy+ismth),2)/real(2*ismth+1)
end do


!!! artificial amplification 

h=h*hfac


return

end subroutine etopo1_init

!===================================================================================!

subroutine vec_to_deg(vecx,vecy,rotdeg) !!! -180 to 180

include 'Pcon.h'

real(4)::vecx,vecy,rotdeg

if(vecx.eq.0.0)then
 if (vecy.gt.0.0) rotdeg=90.0
 if (vecy.lt.0.0) rotdeg=-90.0
 return
else
 rotdeg=atan(vecy/vecx)/drad
! if (vecx.gt.0.0.and.vecy.lt.0.0) rotdeg=rotdeg+180.0
 if (vecx.lt.0.0.and.vecy.gt.0.0) rotdeg=rotdeg+180.0
 if (vecx.lt.0.0.and.vecy.lt.0.0) rotdeg=180.0-rotdeg
end if
return
end subroutine vec_to_deg

!===================================================================================!

subroutine load(ilonl,ilatl,h)

integer,intent(in)::ilonl,ilatl
real(4),intent(out)::h(60,60)

character*3::clatl,clatr,i2clat
character*4::clonl,clonr,i2clon

character*60::dir_etopo1='/work1/amemiya/Data/TOPO/ETOPO1/ascii/inv/'

ilonr=ilonl+1
ilatr=ilatl+1


clonl=i2clon(ilonl)
clonr=i2clon(ilonr)
clatl=i2clat(ilatl)
clatr=i2clat(ilatr)

open(11,file=trim(dir_etopo1)//clatl//'-'//clatr//'/'//clonl//'-'//clonr)
do ilat=1,60
do ilon=1,60
 read(11,*)dumy1,dumy2,h(ilon,ilat)
end do
end do
close(11)

return
end subroutine load
!===================================================================================!

character*4 function i2clon(ilon)

character*4::temp

write(temp,'(I4)') 1000+abs(ilon)
select case(ilon.ge.0)
 case(.true.);write(i2clon,'(A3,A1)') temp(2:4),'E'
 case(.false.);write(i2clon,'(A3,A1)') temp(2:4),'W'
end select

end function i2clon

!===================================================================================!

character*3 function i2clat(ilat)

character*3::temp

write(temp,'(I3)') 100+abs(ilat)
select case(ilat.ge.0)
 case(.true.);write(i2clat,'(A2,A1)') temp(2:3),'N'
 case(.false.);write(i2clat,'(A2,A1)') temp(2:3),'S'
end select

end function i2clat

!===================================================================================!
