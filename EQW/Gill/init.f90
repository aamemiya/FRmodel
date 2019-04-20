subroutine init
use setup

!!! 
!!! test heating -- steady
efold=10.0
h0_day=10.0
zloc_h=7.0e3
zlen_h=5.0e3
yloc_h=0.0e3
ylen_h=1000.0e3
xloc_h=0.0e3
xlen_h=4000.0e3

!! namelist input 
open(11,file="calc.conf")
 read(11,nml=gill_nml)
close(11)

if (igill.eq.1)then
!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!
!!! Gill-type

write(*,*) ' Gill (1980) settings. '

 h0_day=1.0

 ylen0=2185.0e3 !!! 1st vertical mode
   ce0=47.7    !!! 1st vertical mode
 yloc_h=0.0e3
! yloc_h=ylen0*0.8
 ylen_h=ylen0
 xloc_h=0.0e3
 xlen_h=ylen0  * 2.0 !!! Gill 1980
 eps0  =  0.1 * ce0/ylen0 !!! Gill 1980 
 efold = 1.0/eps0/86400.0
 write(*,*) 'equivalent e-folding time',efold,'day'


!!! vertical 1st mode only
do ix=1,nx; do iy=1,ny; do iz=1,nlev
 if (abs(axx(ix)).le.xlen_h) heat(:,ix,iy,iz) = h0_day * exp(-0.5*((axy(iy)-yloc_h)/ylen_h)**2) * cos(0.5*pi*axx(ix)/xlen_h) * sin(pi*zlev(iz)/ztop)
end do;end do;end do

! H=t00*rd/grav
! write(*,*) 'scale height = ',H,'(m)'
 H=999.9e10
 write(*,*) ' Boussinesq approximation 1/H -> 0 '

!!! KOKOMADE
!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!
elseif (ifile.eq.1)then
 H=t00*rd/grav
 write(*,*) 'scale height = ',H,'(m)'
 call inputforcing
else

 !! gaussian heating in the troposphere

 H=t00*rd/grav

!H=999.9e10

 write(*,*) 'scale height = ',H,'(m)'

 do it=1,nt; do ix=1,nx; do iy=1,ny; do ilev=1,nlev
  heat(it,ix,iy,ilev)=h0_day * exp(-((zlev(ilev)-zloc_h)/zlen_h)**2) * exp(-((axy(iy)-yloc_h)/ylen_h)**2) * exp(-((axx(ix)-xloc_h)/xlen_h)**2) * cos(2.0*pi * axt(it)/(period0*86400.0))
 end do;end do;end do;end do


end if


!var_yz=db(nt/2,nx/2,:,:) * t00/grav * 86400.0
!call quicklook_yz('input (K/day)')

eps0=1.0/(86400.0*efold)




return
end subroutine init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! input heat(:,:,:,:) from the external file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine inputforcing
use setup

include 'include/Grid_MERRA_input.h'
include 'netcdf.inc'

!!! define MERRA 4-D grid 
real(4),allocatable::heat_merra(:,:,:,:)

real(4),allocatable::work1_merra(:,:,:,:)
real(4),allocatable::work2_merra(:,:,:,:)

real(4),allocatable::band_low(:),band_high(:)

character*100::ncdir, ncfile
character*100::cvar
namelist /input_nml/ ncdir, ncfile, cvar, mtime

real(4),parameter::rmiss_ncdf=9.9692100E+36 


open(11,file='calc.conf')
 read(11,nml=input_nml)
close(11)


allocate(heat_merra(mlon,mlat,mlev,mtime))
allocate(work1_merra(mlon,mlat,mlev,mtime))
allocate(work2_merra(mlon,mlat,mlev,mtime))

write(*,*) 'heating data from the external file.'

!!! read MERRA data

  istat=NF_OPEN(trim(ncdir)//trim(ncfile),NF_NOWRITE,idnc)
  if (istat.ne.NF_NOERR) call cabort(NF_STRERROR(istat)//trim(ncdir)//trim(ncfile))
  write(*,*) 'reading ',trim(ncfile),' ...'
  istat=NF_INQ_VARID(idnc,trim(cvar),idvar)
  istat=NF_GET_VAR_REAL(idnc,idvar,heat_merra)
  istat=NF_CLOSE(idnc)
  if (istat.ne.NF_NOERR) call cabort(NF_STRERROR(istat))

where(heat_merra.eq.rmiss_ncdf)heat_merra=0.0


!!!
heat_merra = heat_merra * 86400.0

!!! average or filtering 
write(*,*) 'averaging / filtering ...'

!!! TORI AEZU average
!heat_merra(:,:,:,1)=sum(heat_merra,4)/real(mtime)
!do itime=2,mtime
! heat_merra(:,:,:,itime)=heat_merra(:,:,:,1)
!end do


!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! KANTAN band_pass !!!!
nday_high=20
nday_low=40

!!! test
 write(*,*) maxval(heat_merra),minval(heat_merra)
! write(*,*) maxloc(heat_merra),minloc(heat_merra)

do itime=1,mtime
 itimes=max(itime-nday_high/2,1)
 itimee=min(itime+nday_high/2,mtime)
 work1_merra(:,:,:,itime)=sum(heat_merra(:,:,:,itimes:itimee),4)/real(itimee-itimes+1)
 itimes=max(itime-nday_low/2,1)
 itimee=min(itime+nday_low/2,mtime)
 work2_merra(:,:,:,itime)=sum(heat_merra(:,:,:,itimes:itimee),4)/real(itimee-itimes+1)
end do

heat_merra=work1_merra-work2_merra

!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! test
 write(*,*) maxval(heat_merra),minval(heat_merra)
! write(*,*) maxloc(heat_merra),minloc(heat_merra)
!stop

!!! interpolation
write(*,*) 'pre-processing ...'
call intp(mtime,heat_merra,heat)

!!! test
 write(*,*) maxval(heat),minval(heat)
! write(*,*) maxloc(heat),minloc(heat)



!!! 
deallocate(heat_merra)
deallocate(work1_merra)
deallocate(work2_merra)


return
end subroutine inputforcing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine cabort(cmsg)
character(*),intent(in)::cmsg


write(*,*) 'cabort'
write(*,*) trim(cmsg)
stop

return
end subroutine cabort
!==============================================================!
subroutine intp(mtime,array_merra,array)
use setup
include 'include/Grid_MERRA_input.h'


real(4),intent(in)::array_merra(mlon,mlat,mlev,mtime)
real(4),intent(out)::array(nt,nx,ny,nlev)

real(4)::array_merra_nt(mlon,mlat,mlev,nt) !!! work
real(4)::array_merra_nlevnt(mlon,mlat,nlev,nt) !!! work
real(4)::array_merra_nynlevnt(mlon,ny,nlev,nt) !!! work
real(4)::array_merra_nxnynlevnt(nx,ny,nlev,nt) !!! work

real(4)::zmlev(mlev)
real(4)::ymlat(mlat)
real(4)::xmlon(mlon)
real(4)::tmsec(mtime)

namelist /intp_nml/ axlon_base, axlon_cut_l, axlon_cut_r, axlat_cut_l, axlat_cut_r, xtaper, ytaper
namelist /intp_time_nml/ intv_sec,istart_data,icount_data
open(11,file='calc.conf')
 read(11,nml=intp_nml)
 read(11,nml=intp_time_nml)
close(11)

zmlev = -H * log(1.0e2*pmlev/ps0)
!write(*,*) "zmlev",zmlev

ymlat=er*drad*axmlatSN
!write(*,*) "ymlat",ymlat

xmlon=er*drad*(axmlon-axlon_base)
!write(*,*) "xmlon",xmlon


!!! Time intp
!!! tori aezu
!do it=1,nt
! array_merra_nt(:,:,:,it)=array_merra(:,:,:,1)
!end do


sec_datalen = real(intv_sec * (mtime-1))

tl_data = - 0.5 * sec_datalen
tr_data =   0.5 * sec_datalen

if (tr_data.ge.tr )then
 write(*,*) 'ERROR:: too large data size.  tr_data=',tr_data, 'tr=',tr 
 stop
end if

tmsec =(/( tl_data + real(intv_sec*(i-1)), i=1,mtime)/)
itldata=mtime
itrdata=1
do it=1,nt
 imtl=iblkge(tmsec,mtime,axt(it))
 if (imtl.ge.1.and.imtl.le.mtime-1)then
 if (it.lt.itldata) itldata=it
 if (it.gt.itrdata) itrdata=it
 fact=(axt(it)-tmsec(imtl)) / (tmsec(imtl+1)-tmsec(imtl))
 array_merra_nt(:,:,:,it)=(1.0-fact)*array_merra(:,:,:,imtl) &
                                +    fact *array_merra(:,:,:,imtl+1)
 else
 array_merra_nt(:,:,:,it)=0.0
 end if
end do


!!! tapering in time
ttaper=  86400 * 4.0
nttp = 0.5 * ttaper / (axt(2)-axt(1))

do it=max(itldata-nttp,1),itldata+nttp
 array_merra_nt(:,:,:,it) = array_merra_nt(:,:,:,it) * (0.5 * ( 1.0 + sin(pi*(axt(it)-axt(itldata))/ttaper) ) )
end do
do it=itrdata-nttp,min(itrdata-nttp,nt)
 array_merra_nt(:,:,:,it) = array_merra_nt(:,:,:,it) * (0.5 * ( 1.0 + sin(pi*(-axt(it)+axt(itrdata))/ttaper) ) )
end do






!!! Z intp

array_merra_nt(:,:,:,:)=array_merra_nt(:,:,mlev:1:-1,:)

do ilev=1,nlev
 imlevl=iblkge(zmlev,mlev,zlev(ilev))
 if (imlevl.ge.1.and.imlevl.le.mlev-1)then
 fact=(zlev(ilev)-zmlev(imlevl)) / (zmlev(imlevl+1)-zmlev(imlevl))
 array_merra_nlevnt(:,:,ilev,:)=(1.0-fact)*array_merra_nt(:,:,imlevl  ,:) &
                                +    fact *array_merra_nt(:,:,imlevl+1,:)
 else
 array_merra_nlevnt(:,:,ilev,:)=0.0
 end if
end do

!!! Y intp

array_merra_nlevnt(:,:,:,:)=array_merra_nlevnt(:,mlat:1:-1,:,:)

array_merra_nynlevnt=0.0
iyldata=ny
iyrdata=1
do iy=1,ny
 imlatl=iblkge(ymlat,mlat,axy(iy))
 if(imlatl.ge.1.and.imlatl.le.mlat-1)then 
 if(axmlatSN(imlatl+1).ge.axlat_cut_l.and.  &
    axmlatSN(imlatl)  .le.axlat_cut_r        ) then
 if (iy.lt.iyldata) iyldata=iy
 if (iy.gt.iyrdata) iyrdata=iy
 fact=(axy(iy)-ymlat(imlatl)) / (ymlat(imlatl+1)-ymlat(imlatl))
 array_merra_nynlevnt(:,iy,:,:)=(1.0-fact)*array_merra_nlevnt(:,imlatl  ,:,:) &
                                +    fact *array_merra_nlevnt(:,imlatl+1,:,:)
 end if
 end if
end do


!!! X intp
array_merra_nxnynlevnt=0.0
ixldata=nx
ixrdata=1
do ix=1,nx
 imlonl=iblkge(xmlon,mlon,axx(ix))
 if (imlonl.ge.1.and.imlonl.le.mlon-1)then
   if (axmlon(imlonl+1).ge.axlon_cut_l.and.  &
       axmlon(imlonl)  .le.axlon_cut_r        ) then
 if (ix.lt.ixldata) ixldata=ix
 if (ix.gt.ixrdata) ixrdata=ix
 fact=(axx(ix)-xmlon(imlonl)) / (xmlon(imlonl+1)-xmlon(imlonl))
 array_merra_nxnynlevnt(ix,:,:,:)=(1.0-fact)*array_merra_nynlevnt(imlonl  ,:,:,:) &
                                +      fact *array_merra_nynlevnt(imlonl+1,:,:,:)
    end if
 end if
end do

do it=1,nt
 array(it,:,:,:)= array_merra_nxnynlevnt(:,:,:,it)
end do

!!! tapering in X and Y

nxtp = 0.5 * xtaper / (axx(2)-axx(1))
nytp = 0.5 * ytaper / (axy(2)-axy(1))

do ix=max(ixldata-nxtp,1),ixldata+nxtp
 array(:,ix,:,:) = array(:,ix,:,:) * (0.5 * ( 1.0 + sin(pi*(axx(ix)-axx(ixldata))/xtaper) ) )
end do
do ix=ixrdata-nxtp,min(ixrdata-nxtp,nx)
 array(:,ix,:,:) = array(:,ix,:,:) * (0.5 * ( 1.0 + sin(pi*(-axx(ix)+axx(ixrdata))/xtaper) ) )
end do

do iy=max(iyldata-nytp,1),iyldata+nytp
 array(:,iy,:,:) = array(:,iy,:,:) * (0.5 * ( 1.0 + sin(pi*(axy(iy)-axy(iyldata))/ytaper) ) )
end do
do iy=iyrdata-nytp,min(iyrdata-nytp,ny)
 array(:,iy,:,:) = array(:,iy,:,:) * (0.5 * ( 1.0 + sin(pi*(-axy(iy)+axy(iyrdata))/ytaper) ) )
end do


return
end subroutine intp
!==============================================================!
