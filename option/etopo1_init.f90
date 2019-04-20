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

program etopo1_init
use EZspline_obj
use EZspline

implicit real(a-h,o-z)

include 'Pcon.h'

namelist /etopo1_init_nml/ vlonc,vlatc,vudir,vxlen,vylen,nx,ny
namelist /etopo1_quicklook_nml/ iout,factor,valmx,valmn,dval



real(4),allocatable::vlon_grid(:,:),vlat_grid(:,:)
real(4),allocatable::rx(:,:),ry(:,:),rz(:,:)
real(4),allocatable::rxx(:,:),ryy(:,:),rzz(:,:),rhh(:,:)
real(4),allocatable::h_lonlat(:,:),axlon(:),axlat(:)
real(4),allocatable::h_xy(:,:),axx(:),axy(:)
real(4),allocatable::vx_isolon(:,:),vy_isolon(:,:),vx_isolat(:,:),vy_isolat(:,:)

type(EZspline2_r4)::h_spl


character*20::title1='TIBET'
character*20::title2='ETOPO1'
integer::itpats(40)
real(4)::vtlevs(41)

open(11,file='calc.conf')
read(11,nml=etopo1_init_nml)
read(11,nml=etopo1_quicklook_nml)
close(11)


allocate(vlon_grid(nx,ny),vlat_grid(nx,ny),h_xy(nx,ny),axx(nx),axy(ny))
allocate(rx(nx,ny),ry(nx,ny),rz(nx,ny))
allocate(rxx(nx,ny),ryy(nx,ny),rzz(nx,ny),rhh(nx,ny))

!!! 1. calculate grid points !!!

!!! initial location (rad)
do ix=1,nx
 vlon_grid(ix,:) = -0.5*vxlen/er + (real(ix)-0.5) / real(nx) * vxlen/er
end do
do iy=1,ny
 vlat_grid(:,iy) = -0.5*vylen/er + (real(iy)-0.5) / real(ny) * vylen/er
end do

 write(*,*) 'before'
 write(*,*) vlon_grid(1,1)  /drad,   vlat_grid(1,1)/drad
 write(*,*) vlon_grid(nx,1) /drad,  vlat_grid(nx,1)/drad
 write(*,*) vlon_grid(1,ny) /drad,  vlat_grid(1,ny)/drad
 write(*,*) vlon_grid(nx,ny)/drad, vlat_grid(nx,ny)/drad


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

do iy=1,ny
do ix=1,nx
 call vec_to_deg(rxx(ix,iy),ryy(ix,iy),vlon_grid(ix,iy))
 call vec_to_deg(rhh(ix,iy),rzz(ix,iy),vlat_grid(ix,iy))
end do
end do

 write(*,*) rx(1,1)  , rx(1,1)  ,  rz(1,1)
write(*,*) -cxy*sxz, - sxy*cyz ,  -cxy*sxz*cyz, sxy*syz
 write(*,*) cxy*cxz,  ( -cxy*sxz - sxy*cyz ), ( -cxy*sxz*cyz + sxy*syz )
 write(*,*) sxy*cxz,  ( -sxy*sxz + cxy*cyz ), ( -sxy*sxz*cyz - cxy*syz )
 write(*,*) 
 write(*,*) 'after'
 write(*,*) rxx(1,1)  , rxx(1,1)  ,  rzz(1,1),  rhh(1,1)
 write(*,*) rxx(nx,1) , ryy(nx,1) , rzz(nx,1)
 write(*,*) rxx(1,ny) , ryy(1,ny) , rzz(1,ny)
 write(*,*) rxx(nx,ny), ryy(nx,ny), rzz(nx,ny)

 write(*,*) 
 write(*,*) 'after'
 write(*,*) vlon_grid(1,1)  ,   vlat_grid(1,1)
 write(*,*) vlon_grid(nx,1) ,  vlat_grid(nx,1)
 write(*,*) vlon_grid(1,ny) ,  vlat_grid(1,ny)
 write(*,*) vlon_grid(nx,ny), vlat_grid(nx,ny)



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
do iy=1,ny
 call EZspline_interp(h_spl,nx,vlon_grid(:,iy),vlat_grid(:,iy),h_xy(:,iy),ier) !!! treat as 1-d array
end do
call EZspline_free(h_spl,ier)


!!! 5. let data points below the sea surface be zero altitude 

where(h_xy.lt.0.0) h_xy=0.0

write(*,*) maxval(h_lonlat),minval(h_lonlat)
write(*,*) maxval(h_xy),minval(h_xy)

!!! 6. convolute h with window function to smooth out small scales


!!! 7. remove data near the boundaries using proper cutoff function 


!!! (optional) lat. lon. isolines

allocate(vx_isolon(nx,nlon/60-1),vy_isolon(nx,nlon/60-1),vx_isolat(nx,nlat/60-1),vy_isolat(nx,nlat/60-1))

if (nlon/60 + nlat/60 - 2 .le. ny)then
 do ilon=ilonl+1,ilonr
  ilond=ilon-ilonl
  !! reuse
  rxx(:,ilond) = (/( cos( (real(ilatl)+(real(ix)-0.5)*real(ilatr-ilatl)/real(nx) ) * drad) * cos(real(ilon)*drad) ,ix=1,nx)/)
  ryy(:,ilond) = (/( cos( (real(ilatl)+(real(ix)-0.5)*real(ilatr-ilatl)/real(nx) ) * drad) * sin(real(ilon)*drad) ,ix=1,nx)/)
  rzz(:,ilond) = (/( sin( (real(ilatl)+(real(ix)-0.5)*real(ilatr-ilatl)/real(nx) ) * drad) ,ix=1,nx)/)
 end do
 do ilat=ilatl+1,ilatr
  ilatd=ilat-ilatl
  !! reuse
  rxx(:,nlon/60-1+ilatd) = (/( cos( (real(ilonl)+(real(ix)-0.5)*real(ilonr-ilonl)/real(nx) ) * drad) * cos(real(ilat)*drad) ,ix=1,nx)/)
  ryy(:,nlon/60-1+ilatd) = (/( sin( (real(ilonl)+(real(ix)-0.5)*real(ilonr-ilonl)/real(nx) ) * drad) * cos(real(ilat)*drad) ,ix=1,nx)/)
  rzz(:,nlon/60-1+ilatd) = sin(real(ilat)*drad)
 end do

rx =                cxz*cxy   * rxx +                 cxz*sxy  * ryy +     sxz * rzz
ry = ( -cyz*sxy-syz*sxz*cxy ) * rxx + ( cyz*cxy-syz*sxz*sxy )  * ryy + syz*cxz * rzz
rz=  ( syz*sxy-cyz*sxz*cxy )  * rxx + ( -syz*cxy-cyz*sxz*sxy ) * ryy + cyz*cxz * rzz

rhh=sqrt(rx**2+ry**2)

do iy=1,ny
do ix=1,nx
 call vec_to_deg(rx(ix,iy),ry(ix,iy),vlon_grid(ix,iy))
 call vec_to_deg(rhh(ix,iy),rz(ix,iy),vlat_grid(ix,iy))
end do
end do

  vx_isolon(:,1:nlon/60-1)= er*drad*vlon_grid(:,1:nlon/60-1)
  vy_isolon(:,1:nlon/60-1)= er*drad*vlat_grid(:,1:nlon/60-1)
  vx_isolat(:,1:nlat/60-1)= er*drad*vlon_grid(:,nlon/60:nlon/60-1+nlat/60-1)
  vy_isolat(:,1:nlat/60-1)= er*drad*vlat_grid(:,nlon/60:nlon/60-1+nlat/60-1)


write(*,*) vlon_grid(:,1)
write(*,*) vx_isolon(:,1)*1.0e-3
write(*,*) vy_isolon(:,1)*1.0e-3
write(*,*) vx_isolat(:,1)*1.0e-3
write(*,*) vy_isolat(:,1)*1.0e-3


else
 write(*,*) 'error.'
 write(*,*) nlon/60,nlat/60
 write(*,*) ilonl,ilonr,ilatl,ilatr
end if

deallocate(rx,ry,rz)
deallocate(rxx,ryy,rzz,rhh)

!!! quicklook

axx=(/( -0.5*vxlen + (real(ix)-0.5)/real(nx) * vxlen ,ix=1,nx)/)
axy=(/( -0.5*vylen + (real(iy)-0.5)/real(ny) * vylen ,iy=1,ny)/)

iouts=1
ioute=2
if(iout.eq.1) ioute=1
if(iout.eq.2) iouts=2
do ioutl=iouts,ioute
! *** general settings ***
      call sgiset ('IFONT',1)
!      call swistx ('ICLRMAP',14) ! colormap blue-white-red 
      call swistx ('ICLRMAP',02) 
      call swcmll
      call swcset ('FNAME','figure')
      call swiset ('IWIDTH',1000) ! output window size
      call swiset ('IHEIGHT',800)     
      call gropn(ioutl) 
      call sglset ('LFULL',.TRUE.) ! using fullsize
      call swlset ('LSEP',.fAlSE.) ! psfile numbering
      call slmgn (0.0,0.0,0.0,0.0) ! margin 
      call grfrm 
      call grswnd (-0.5*vxlen,0.5*vxlen,-0.5*vylen,0.5*vylen) ! set window
      aratio=vxlen/vylen
      vpxl=0.20
      vpxr=0.90
      vpyl=max(0.70-0.70/aratio,0.10)
      vpyr=0.70
      if (vpyl.eq.0.10)then
         vpxl=0.55-0.30*aratio
         vpxr=0.55+0.30*aratio
      end if
      !   write(*,*)vpxl,vpxr,vpyl,vpyr
      call grsvpt (vpxl,vpxr,vpyl,vpyr) ! set viewport
      call grstrn (1) ! linear or log
      call grstrf

      call uwsgxa (axx,nx) ! xaxis value
      call uwsgya (axy,ny) ! yaxis value


! *** topo ***

  call udlset ('LMSG',.FALSE.) ! 'contour interval'
  call udrset ('RSIZEL',0.015) ! label


   ntpat=int((valmx-valmn)/dval)+1
   iinc=-(70/ntpat)
   itpats(1:ntpat-1)=(/(999+1000*(94+iinc*(i-1)),i=1,ntpat-1)/) !!! Colormap 02
   itpats(ntpat)=10999
   vtlevs(1:ntpat)=(/( valmn + real(i-1)*dval, i=1,ntpat  )/)
   vtlevs(ntpat+1)=1.0e10
   
   call ueitlv
   call uestln (vtlevs(1:ntpat+1), itpats(1:ntpat), ntpat)
   call uetone(h_xy*factor,nx,nx,ny)
   call dcbar  (vpxr+0.01,vpyl+0.01,vpyr-vpyl-0.02,vtlevs(1:ntpat+1), itpats(1:ntpat), ntpat)

 ! *** coast lines ***
   call udsclv(1.0,3,1,'',-1.0)
   
   call udcntr(h_xy,nx,nx,ny)

! *** longitude and latitude lines ***

   call sglset ('LCLIP',.TRUE.) ! Cliping

 if (ilonr-ilonl.gt.10)then
  incr=2
 else 
  incr=1
 end if
 do ilon=1,ilonr-ilonl-1,incr
  call uulinz(nx,vx_isolon(:,ilon),vy_isolon(:,ilon),3,3)  
 end do


 if (ilatr-ilatl.gt.10)then
  incr=2
 else 
  incr=1
 end if
 do ilat=1,ilatr-ilatl-1,incr
  call uulinz(nx,vx_isolat(:,ilat),vy_isolat(:,ilat),3,3)  
 end do
 

! **** x ,y axis ****
      call uziset ('INDEXT2',3)
      call uziset ('INDEXL1',5)
      call uzrset ('RSIZEL1',0.018)
      call uzrset ('RSIZEC1',0.018)
      call uzrset ('RSIZET2',0.010)
      call uzrset ('RSIZET1',0.004)
!      call uzlset ('LABELXB',.FALSE.)
      call uzlset ('LABELXT',.FALSE.)
      call uzlset ('LOFFSET',.TRUE.)
      call uzrset ('YFACT',1.0e-3)
      call uzrset ('XFACT',1.0e-3)
      call uxsfmt ('(I5)')
      if (vxlen.lt.1000.0e3)then
       call uxaxdv ('B',10.0,50.0)
       call uxaxdv ('T',10.0,50.0)
      else
       call uxaxdv ('B',100.0,500.0)
       call uxaxdv ('T',100.0,500.0)
      end if

      call uxsttl ('B','x (km)',0.0)
      call uysfmt ('(I5)')
      call uzlset ('LABELYR',.FALSE.)
!     call uzlset ('LABELYL',.FALSE.)
      if (vylen.lt.1000.0e3)then
       call uyaxdv ('L',10.0,50.0)
       call uyaxdv ('R',10.0,50.0)
      else
       call uyaxdv ('L',100.0,500.0)
       call uyaxdv ('R',100.0,500.0)
      end if
      call uziset ('IROTCYL',1)
      call uysttl ('L','y (km)',0.0)
      call sglset ('LCLIP',.FALSE.) ! Cliping
      call sgtxzv (0.5*(vpxr+vpxl),vpyr+0.03,trim(title1),0.024,0,0,5) ! title
      call sgtxzv (vpxl+0.03,vpyr+0.05,trim(title2),0.016,0,-1,5) ! title
      call grcls
      end do


end program etopo1_init

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

!--------------------------------------------------------------!

subroutine dcbar(vpxr,vpyl,dylen,vtpats,itpats,ntpat)

integer::itpats(ntpat)
real(4)::vtpats(ntpat+1)

real(4)::xm4(4),ym4(4)
real(4)::xm5(5),ym5(5)
character*4::cvar
character*10::cfact

namelist /etopo1_quicklook_nml/ iout,factor,valmx,valmn,dval

open(11,file='calc.conf')
read(11,nml=etopo1_quicklook_nml)
close(11)


call sgqvpt(vpxl,vpxr,vpyl,vpyr)


dyp=dylen/real(ntpat)

call sglset ('LCLIP',.FALSE.) ! Cliping

!!! Color bar

      do ic=1,ntpat
         xm4(1)=vpxr+0.01
         xm4(2)=xm4(1)+0.02
         xm4(3)=xm4(2)
         xm4(4)=xm4(1)
         ym4(1)=vpyl+real(ic-1)*dyp
         ym4(2)=ym4(1)
         ym4(3)=vpyl+real(ic)*dyp
         ym4(4)=ym4(3)
         itpat=itpats(ic)
         call sgtnzv(4,xm4,ym4,itpat)
      end do
!!! Waku
      xm5=(/vpxr+0.01,vpxr+0.03,vpxr+0.03,vpxr+0.01,vpxr+0.01/)
      ym5=(/vpyl,vpyl,vpyl+dylen,vpyl+dylen,vpyl/)
      call sgplzv(5,xm5,ym5,1,1) 


      if (ntpat.eq.21)then
       do ic=1,21,4
         write(cvar,'(I4)') int(vtpats(ic)) 
         call sgtxzv(xm4(2),vpyl+real(ic-1)*dyp,cvar,0.016,0,-1,3) ! 
       end do
      else !!! TORI AEZU 
      ic=2
         write(cvar,'(I4)') int(vtpats(2)) 
         call sgtxzv(xm4(2),vpyl+real(ic-1)*dyp,cvar,0.016,0,-1,3) ! 
      ic=ntpat
         write(cvar,'(I4)') int(vtpats(ntpat)) 
         call sgtxzv(xm4(2),vpyl+real(ic-1)*dyp,cvar,0.016,0,-1,3) ! 
     end if


         ixfac=-nint(log(factor)/log(10.0)) !!!
         if (ixfac.ne.0)then
         write(cfact,'(A,I2,A)') '*10|',ixfac,'"'
         call sgtxzv(xm4(1),vpyl+real(ntpat)*dyp+0.02,cfact,0.018,0,-1,3) ! 
         end if

return
end subroutine dcbar

!==============================================================!
