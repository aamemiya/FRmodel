
program gill
use setup
use eqw_decomp

real(4)::db(nt,nx,ny,nlev)

real(4)::vb(nt,nx,ny,nlev)
real(4)::vu(nt,nx,ny,nlev)
real(4)::vv(nt,nx,ny,nlev)
real(4)::vw(nt,nx,ny,nlev)
real(4)::vp(nt,nx,ny,nlev)

real(4)::vpv(nt,nx,ny,nlev)
real(4)::vz(nt,nx,ny,nlev)

real(4)::vb_n(nt,nx,ny,nlev)
real(4)::vu_n(nt,nx,ny,nlev)
real(4)::vv_n(nt,nx,ny,nlev)
real(4)::vw_n(nt,nx,ny,nlev)
real(4)::vp_n(nt,nx,ny,nlev)

real(4)::db_h(nh)
real(4)::db_h_tmp(nh)
real(4)::q(ny,nlev)

real(4)::r_Q_txhn(nt,nx,nh,nlev)

complex::c_dq_nh(nh,nlev)

complex::c_Q_okhn(nt,nx,nh,nlev)

complex::c_vq_okhn(nt,nx,nh,nlev)
complex::c_vr_okhn(nt,nx,nh,nlev)

complex::c_vu_okhn(nt,nx,nh,nlev)
complex::c_vv_okhn(nt,nx,nh,nlev)
complex::c_vw_okhn(nt,nx,nh,nlev)
complex::c_vb_okhn(nt,nx,nh,nlev)
complex::c_vp_okhn(nt,nx,nh,nlev)


complex,parameter::c_ur=(1.0,0.0)
complex,parameter::c_ui=(0.0,1.0)

character*40::title

character*200::command

!!! initialize -- specify heating 

 call init

nkcut=int((xr-xl)/xlen_nkcut)
nocut=int((tr-tl)/(peri_nocut*86400.0))

!if (igill.eq.1) nncut=1

write(*,*) 'cutoff:'
write(*,*) 'nkcut : ', nkcut
write(*,*) 'nocut : ', nocut
write(*,*) 'nhcut : ', nhcut
write(*,*) 'nncut : ', nncut
write(*,*)
write(*,*) ' == Here we go == '
write(*,*)

write(*,*) ' forward... '

!if (igill.eq.1)then
! db=heat
!else
 do ilev=1,nlev
  db(:,:,:,ilev)=heat(:,:,:,ilev)*exp(-0.5*zlev(ilev)/H)
 end do
 db=grav/t00*db / 86400.0 !!! K / day --> K / s and Q --> db

 do it=1,nt; do ix=1,nx
  call fft_fwd_sin_trans(ny,db(it,ix,:,1:nlev-1),nlev-1)
 end do;end do
 db(:,:,:,nlev)=0.0
!end if

!do it=1,nt
! write(*,*) 'max min db',maxval(db(it,:,:,:)),minval(db(it,:,:,:))
!end do


do in=1,nncut
 vm=pi*real(in)/ztop
 ce=sqrt(bn2/(vm**2+(0.5/H)**2))
 ylen=sqrt(ce/beta)
! write(*,*) 'ylen ',ylen
! write(*,*) 'ce ',ce
! stop


 db(:,:,:,in)=db(:,:,:,in) * vm**2/(vm**2+(0.5/H)**2) !! factor -- maybe this is not complete (need cos component)



 db(:,:,:,in)=db(:,:,:,in) / bn2 !! b--> Q


!do it=1,nt
! write(*,*) 'max min db 2',maxval(db(it,:,:,:)),minval(db(it,:,:,:))
!end do


 call set_grid_hermite
 do it=1,nt
 do ix=1,nx
       call decomp_fwd(db(it,ix,:,in),r_Q_txhn(it,ix,:,in))
 end do
 end do
end do

!do it=1,nt
! write(*,*) 'max min r_Q',maxval(r_Q_txhn(it,:,:,:))
!end do


c_Q_okhn = c_ur * r_Q_txhn

!do it=1,nt
! write(*,*) 'max min c_Q',maxval(real(c_Q_okhn(it,:,:,:)))
!end do

 call fft_fwd_comp_2d(nt,nx,c_Q_okhn,nh*nlev)

!do it=1,nt
! write(*,*) 'max min c_Q',maxval(real(c_Q_okhn(it,:,:,:)))
!end do

!stop

!
!!!do in=1,nncut


do in=3,3



 vm=pi*real(in)/ztop
 ce=sqrt(bn2/(vm**2+(0.5/H)**2))
 ylen=sqrt(ce/beta)
 call set_grid_hermite


do io=1,min(nt,2*(nocut+1))
do ik=1,min(nx,2*(nkcut+1))
 
!!! Kelvin wave
!!!   c_vq_okhn(io,ik,1,in)= - c_Q_okhn(io,ik,1,in) * c_ui / (axo(io)-ce*axk(ik)+c_ui*eps0) /ce

 do ih=3,min(nhmax-1,nhcut)
! do ih=2,2
  c_vq_okhn(io,ik,ih,in)= - ( real(ih-1)*c_Q_okhn(io,ik,ih+1,in) + sqrt(real(ih*(ih-1)))*c_Q_okhn(io,ik,ih-1,in) ) * c_ui / ( real(2*ih-1)*(axo(io)+c_ui*eps0)+ce*axk(ik)) /ce
  c_vr_okhn(io,ik,ih,in)=sqrt(1.0+1.0/real(ih))*c_vq_okhn(io,ik,ih,in)
 end do
  c_vu_okhn(io,ik,1,in)=0.5*(c_vq_okhn(io,ik,1,in)-c_vr_okhn(io,ik,2,in))
  c_vp_okhn(io,ik,1,in)=ce*0.5*(c_vq_okhn(io,ik,1,in)+c_vr_okhn(io,ik,2,in))
  c_vu_okhn(io,ik,2,in)=-0.5*c_vr_okhn(io,ik,3,in)
  c_vp_okhn(io,ik,2,in)=ce*0.5*c_vr_okhn(io,ik,3,in)
 do ih=3,min(nhmax-1,nhcut)
  c_vu_okhn(io,ik,ih,in)=0.5*(c_vq_okhn(io,ik,ih-1,in)-c_vr_okhn(io,ik,ih+1,in))
  c_vp_okhn(io,ik,ih,in)=ce*0.5*(c_vq_okhn(io,ik,ih-1,in)+c_vr_okhn(io,ik,ih+1,in))
 end do

do ih=2,min(nhmax-1,nhcut)
  c_vv_okhn(io,ik,ih,in)= (- c_ui * (axo(io)-ce*axk(ik)+c_ui*eps0) * c_vq_okhn(io,ik,ih,in) + c_Q_okhn(io,ik,ih+1,in)/ce) / sqrt(beta*ce) / sqrt(2.0*real(ih-1))
end do

!!! CORRECTED 16/11/09
! c_vb_okhn(io,ik,:,in) =  vm * c_vp_okhn(io,ik,:,in)
 c_vb_okhn(io,ik,:,in) = - vm * c_vp_okhn(io,ik,:,in)

 c_vw_okhn(io,ik,:,in)=   ( c_Q_okhn(io,ik,:,in)*bn2 + c_ui * (axo(io) + c_ui*eps0 ) * c_vb_okhn(io,ik,:,in) ) / bn2
!!! c_vw_okhn(io,ik,:,in)= c_Q_okhn(io,ik,:,in) - ( c_ui * (axo(io) + c_ui*eps0 ) * c_vb_okhn(io,ik,:,in) ) / bn2
end do
end do
end do
 
!do it=1,nt
! write(*,*) 'max c_vb',maxval(real(c_vb_okhn(it,:,:,:)))
!end do


write(*,*) ' backward... '

!!! backword transform
 call fft_bwd_comp_2d(nt,nx,c_vq_okhn,nh*nlev)
 call fft_bwd_comp_2d(nt,nx,c_vu_okhn,nh*nlev)
 call fft_bwd_comp_2d(nt,nx,c_vv_okhn,nh*nlev)
 call fft_bwd_comp_2d(nt,nx,c_vw_okhn,nh*nlev)
 call fft_bwd_comp_2d(nt,nx,c_vp_okhn,nh*nlev)
 call fft_bwd_comp_2d(nt,nx,c_vb_okhn,nh*nlev)


!do it=1,nt
! write(*,*) 'max c_vb',maxval(real(c_vb_okhn(it,:,:,:)))
!end do


!stop

vu=0.0; vv=0.0; vw=0.0; vb=0.0; vp=0.0

do in=1,nncut

 vm=pi*real(in)/ztop
 ce=sqrt(bn2/(vm**2+(0.5/H)**2))
 ylen=sqrt(ce/beta)
 alpha=asin((0.5/H)/sqrt(vm**2+(0.5/H)**2))

 call set_grid_hermite

 do it=1,nt
 do ix=1,nx
  call decomp_bwd(real(c_vu_okhn(it,ix,:,in)),vu_n(it,ix,:,in))
  call decomp_bwd(real(c_vv_okhn(it,ix,:,in)),vv_n(it,ix,:,in))
  call decomp_bwd(real(c_vw_okhn(it,ix,:,in)),vw_n(it,ix,:,in))
  call decomp_bwd(real(c_vp_okhn(it,ix,:,in)),vp_n(it,ix,:,in))
  call decomp_bwd(real(c_vb_okhn(it,ix,:,in)),vb_n(it,ix,:,in))
 end do
 end do

 do ilev=1,nlev
  vu(:,:,:,ilev)=vu(:,:,:,ilev) + exp( zlev(ilev)*0.5/H) * cos(vm*zlev(ilev)+alpha) * vu_n(:,:,:,in)
  vv(:,:,:,ilev)=vv(:,:,:,ilev) + exp( zlev(ilev)*0.5/H) * cos(vm*zlev(ilev)+alpha) * vv_n(:,:,:,in)
  vp(:,:,:,ilev)=vp(:,:,:,ilev) + exp(-zlev(ilev)*0.5/H) * cos(vm*zlev(ilev)+alpha) * vp_n(:,:,:,in)
  vw(:,:,:,ilev)=vw(:,:,:,ilev) + exp( zlev(ilev)*0.5/H) * sin(vm*zlev(ilev)) * vw_n(:,:,:,in)
  vb(:,:,:,ilev)=vb(:,:,:,ilev) + exp( zlev(ilev)*0.5/H) * sin(vm*zlev(ilev)) * vb_n(:,:,:,in)
 end do
end do


!!! *** CORRECTED 16/11/5 ***

  vz = - vb / bn2

! do ilev=1,nlev
!  vz(:,:,:,ilev)= t00*rd/grav * vp(:,:,:,ilev) / ( ps0 * exp(-zlev(ilev)/H))
! end do

  call diagnose_pv(vu,vv,vb,vpv)

  if (icuth.eq.1) call heat_cutoff

do ix=nx/2,nx/2
! var_yz=vpv(nt/2,ix-10,:,:) * 1.0e6
! var_yz=vp(nt/2,ix-10,:,:) *1.0e-2
! var_yz=vz(nt/2,ix-10,:,:) 
! var_yz=vb(nt/2,ix-10,:,:)
 var_yz=vw(1,ix+10,:,:)
 var_yz_2=heat(nt/2,ix,:,:) 
 write(title,'(A,I5,A)')'x= ',int((1.0e-3)*axx(ix)),' (km)'
 call quicklook_yz(title)
end do



!do it=11,nt-10
do it=nt/2,nt/2
!do ilev=1,nlev-5
do ilev=nlev-8,nlev-8
! var_xy=vp(it,:,:,ilev) *1.0e-2
! var_xy=vz(nt/2,:,:,ilev) 
 var_xy=vpv(nt/2,:,:,ilev) * 1.0e6
! var_xy=vu(nt/2,:,:,ilev) 
! var_xy=vv(nt/2,:,:,ilev) 
! var_xy=vz(nt/2,:,:,ilev) 
! var_xy_2=heat(nt/2,:,:,ilev) 
 var_xy_2=heat(it,:,:,ilev) 
! write(title,'(A,I5,A)')'z= ',int((1.0e-3)*zlev(ilev)),' (km)'
 if(igill.eq.1)then
 write(title,'(A,I5,A,I4,A)')'z= ',int((1.0e-3)*zlev(ilev)),' (km)'
 else
 write(title,'(A,I5,A,I4,A)')'z= ',int((1.0e-3)*zlev(ilev)),' (km)  t=',int((axt(it)-axt(1))/86400.0),' (day)'
 end if
! write(*,*)maxval(var_xy)
 call quicklook_xy(title)
end do

end do



!!! Output
open(11,file='calc.conf')
 read(11,nml=output_nml)
close(11)

call mkdir_f90('./save/'//trim(cname))
write(command,'(a)') 'cp Grid.h ./save/'//trim(cname)//'/'
lencom=len(trim(command))
call system(command(1:lencom))

write(command,'(a)') 'cp calc.conf ./save/'//trim(cname)//'/'
lencom=len(trim(command))
call system(command(1:lencom))

write(command,'(a)') 'cp Pcon.h ./save/'//trim(cname)//'/'
lencom=len(trim(command))
call system(command(1:lencom))


open(21,file='./save/'//trim(cname)//'/heat',form='unformatted')
open(22,file='./save/'//trim(cname)//'/u',form='unformatted')
open(23,file='./save/'//trim(cname)//'/v',form='unformatted')
open(24,file='./save/'//trim(cname)//'/w',form='unformatted')
open(25,file='./save/'//trim(cname)//'/p',form='unformatted')
open(26,file='./save/'//trim(cname)//'/b',form='unformatted')
open(27,file='./save/'//trim(cname)//'/pv',form='unformatted')
if(itimev.eq.1)then
 write(21) heat
 write(22) vu
 write(23) vv
 write(24) vw
 write(25) vp
 write(26) vb
 write(27) vpv
else
 write(21) heat(1,:,:,:)
 write(22) vu(1,:,:,:)
 write(23) vv(1,:,:,:)
 write(24) vw(1,:,:,:)
 write(25) vp(1,:,:,:)
 write(26) vb(1,:,:,:)
 write(27) vpv(1,:,:,:)
end if
close(21)
close(22)
close(23)
close(24)
close(25)
close(26)
close(27)
!!!


end program gill
!======================================!
subroutine diagnose_pv(arrayu,arrayv,arrayb,arraypv)
use setup

real(4)::arrayu(nt,nx,ny,nlev)
real(4)::arrayv(nt,nx,ny,nlev)
real(4)::arrayb(nt,nx,ny,nlev)
real(4)::arraypv(nt,nx,ny,nlev)

real(4)::dudz(nt,nx,ny,nlev)
real(4)::dudy(nt,nx,ny,nlev)
real(4)::dvdx(nt,nx,ny,nlev)
real(4)::dbdz(nt,nx,ny,nlev)
real(4)::dbdy(nt,nx,ny,nlev)

real(4)::rho(nlev)

 rho = ps0 * exp(-(1.0-rd/cp)*zlev/H) /rd /t00

do ix=1,nx
 ixr=min(ix+1,nx)
 ixl=max(ix-1,1)
 dvdx(:,ix,:,:)=(arrayv(:,ixr,:,:)-arrayv(:,ixl,:,:)) / (axx(ixr)-axx(ixl))
end do

do iy=1,ny
 iyr=min(iy+1,ny)
 iyl=max(iy-1,1)
 dudy(:,:,iy,:)=(arrayu(:,:,iyr,:)-arrayu(:,:,iyl,:)) / (axy(iyr)-axy(iyl))
 dbdy(:,:,iy,:)=(arrayb(:,:,iyr,:)-arrayb(:,:,iyl,:)) / (axy(iyr)-axy(iyl))
end do

do ilev=1,nlev
 ilevt=min(ilev+1,nlev)
 ilevb=max(ilev-1,1)
 dudz(:,:,:,ilev)=(arrayu(:,:,:,ilevt)-arrayu(:,:,:,ilevb)) / (zlev(ilevt)-zlev(ilevb))
 dbdz(:,:,:,ilev)=(arrayb(:,:,:,ilevt)-arrayb(:,:,:,ilevb)) / (zlev(ilevt)-zlev(ilevb))
end do




do iy=1,ny
do ilev=1,nlev
! arraypv(:,:,iy,ilev)= t00/grav / rho(ilev) * ( (bn2+dbdz(:,:,iy,ilev)) * (beta *axy(iy) + dvdx(:,:,iy,ilev) -dudy(:,:,iy,ilev))  + dbdy(:,:,iy,ilev) * dudz(:,:,iy,ilev))
! arraypv(:,:,iy,ilev)= t00/grav / rho(ilev) * ( (bn2) * (beta *axy(iy) + dvdx(:,:,iy,ilev) -dudy(:,:,iy,ilev)) )

!!!
! arraypv(:,:,iy,ilev)= t00/grav / rho(ilev) * ( bn2 * (beta *axy(iy) -dudy(:,:,iy,ilev)) +dbdz(:,:,iy,ilev) * (beta *axy(iy) ))

 arraypv(:,:,iy,ilev)= t00/grav / rho(ilev) * ( bn2 * ( -dudy(:,:,iy,ilev)) +dbdz(:,:,iy,ilev) * (beta *axy(iy) ))
end do
end do



return
end subroutine diagnose_pv
!======================================!
subroutine heat_cutoff
use setup
use eqw_decomp

include 'include/Grid_MERRA_input.h'

real(4)::heat_n(nt,nx,ny,nlev)
real(4)::r_work(nt,nx,nh,nlev)
complex::c_work(nt,nx,nh,nlev)
real(4)::ymlat(mlat)
real(4)::xmlon(mlon)
namelist /intp_nml/ axlon_base, axlon_cut_l, axlon_cut_r, axlat_cut_l, axlat_cut_r, xtaper, ytaper

 open(11,file='calc.conf')
  read(11,nml=intp_nml)
 close(11)

ymlat=er*drad*axmlatSN
xmlon=er*drad*(axmlon-axlon_base)

do ilev=1,nlev
 heat(:,:,:,ilev)=heat(:,:,:,ilev) + exp( -zlev(ilev)*0.5/H) 
end do

 do it=1,nt; do ix=1,nx
  call fft_fwd_sin_trans(ny,heat(it,ix,:,1:nlev-1),nlev-1)
 end do;end do
 heat(:,:,:,nlev)=0.0
 
r_work=0.0
do in=1,nncut
 vm=pi*real(in)/ztop
 ce=sqrt(bn2/(vm**2+(0.5/H)**2))
 ylen=sqrt(ce/beta)

 call set_grid_hermite
 do it=1,nt
 do ix=1,nx
   call decomp_fwd(heat(it,ix,:,in),r_work(it,ix,:,in))
 end do
 end do
if (nhcut.lt.nhmax) r_work(:,:,nhcut+1:nhmax,in)=0.0
end do

c_work = (1.0,0.0) * r_work

 call fft_fwd_comp_2d(nt,nx,c_work,nh*nlev)

!!! cutoff
if (2*(nocut+1).lt.nt) c_work(2*(nocut+1)+1:nt,:,:,:)=0.0
if (2*(nkcut+1).lt.nx) c_work(:,2*(nkcut+1)+1:nx,:,:)=0.0

!!! for display
!c_work(:,1,:,:)=0.0


!!! backword transform
 call fft_bwd_comp_2d(nt,nx,c_work,nh*nlev)

r_work=real(c_work)

do in=1,nncut
 vm=pi*real(in)/ztop
 ce=sqrt(bn2/(vm**2+(0.5/H)**2))
 ylen=sqrt(ce/beta)

 call set_grid_hermite

 do it=1,nt
 do ix=1,nx
  call decomp_bwd(r_work(it,ix,:,in),heat_n(it,ix,:,in))
 end do
 end do

 do ilev=1,nlev
  heat(:,:,:,ilev)=heat(:,:,:,ilev) + exp( zlev(ilev)*0.5/H) * sin(vm*zlev(ilev)) * heat_n(:,:,:,in)
 end do
end do


!!! Y intp
iyldata=ny
iyrdata=1
do iy=1,ny
 imlatl=iblkge(ymlat,mlat,axy(iy))
 if(imlatl.ge.1.and.imlatl.le.mlat-1)then 
 if(axmlatSN(imlatl+1).ge.axlat_cut_l.and.  &
    axmlatSN(imlatl)  .le.axlat_cut_r        ) then
 if (iy.lt.iyldata) iyldata=iy
 if (iy.gt.iyrdata) iyrdata=iy
 else
 heat(:,:,iy,:)=0.0
 end if
 else
 heat(:,ix,:,:)=0.0
 end if
end do

!!! X intp
ixldata=nx
ixrdata=1
do ix=1,nx
 imlonl=iblkge(xmlon,mlon,axx(ix))
 if (imlonl.ge.1.and.imlonl.le.mlon-1)then
   if (axmlon(imlonl+1).ge.axlon_cut_l.and.  &
       axmlon(imlonl)  .le.axlon_cut_r        ) then
 if (ix.lt.ixldata) ixldata=ix
 if (ix.gt.ixrdata) ixrdata=ix
 else
 heat(:,ix,:,:)=0.0
    end if
 else
 heat(:,ix,:,:)=0.0
 end if
end do

!!! tapering in X and Y

nxtp = 0.5 * xtaper / (axx(2)-axx(1))
nytp = 0.5 * xytaper / (axy(2)-axy(1))

do ix=max(ixldata-nxtp,1),ixldata+nxtp
 heat(:,ix,:,:) = heat(:,ix,:,:) * (0.5 * ( 1.0 + sin(pi*(axx(ix)-axx(ixldata))/xtaper) ) )
end do
do ix=ixrdata-nxtp,min(ixrdata-nxtp,nx)
 heat(:,ix,:,:) = heat(:,ix,:,:) * (0.5 * ( 1.0 + sin(pi*(axx(ix)-axx(ixrdata))/xtaper) ) )
end do

do iy=max(iyldata-nytp,1),iyldata+nytp
 heat(:,iy,:,:) = heat(:,iy,:,:) * (0.5 * ( 1.0 + sin(pi*(axy(iy)-axy(iyldata))/ytaper) ) )
end do
do iy=iyrdata-nytp,min(iyrdata-nytp,ny)
 heat(:,iy,:,:) = heat(:,iy,:,:) * (0.5 * ( 1.0 + sin(pi*(axy(iy)-axy(iyrdata))/ytaper) ) )
end do





return
end subroutine heat_cutoff
!======================================!
subroutine quicklook_yz(title)
use setup

character(*)::title
character*20::title1,title2

real(4)::u_snap(ny,nlev)
real(4)::v_snap(ny,nlev)
real(4)::w_snap(ny,nlev)
real(4)::h_snap(ny,nlev)
real(4)::b_snap(ny,nlev)
real(4)::z_snap(ny,nlev)


title1=title

ystics=1000.0
ymtics=2000.0

! output
iouts=1
ioute=2
if(iout.eq.1) ioute=1
if(iout.eq.2) iouts=2
do ioutl=iouts,ioute
! *** general settings ***
      call sgiset ('IFONT',1)
      call swistx ('ICLRMAP',14) ! colormap blue-white-red 
      call swcmll
      call swcset ('FNAME','figure')
      call swiset ('IWIDTH',1000) ! output window size
      call swiset ('IHEIGHT',800)     
      call gropn(ioutl) 
      call sglset ('LFULL',.TRUE.) ! using fullsize
      call swlset ('LSEP',.TRUE.) ! psfile numbering
      call slmgn (0.0,0.0,0.0,0.0) ! margin 
      call grfrm 
      call grswnd (1.0e-3*yl,1.0e-3*yr,1.0e-3*zbot,1.0e-3*ztop) ! set window
      call grsvpt (0.15,0.85,0.13,0.73) ! set viewport !!! 1:1
      call grstrn (1) ! linear or log
      call grstrf


  zasp= (ztop-zbot) / (xr-xl) 

! **** contour ****

      call udlset ('LMSG',.FALSE.) ! 'contour interval'
      call udrset ('RSIZEL',0.015) ! label
      call udiset ('INDXMJ',13) ! major line
      call udiset ('INDXMN',11) ! minor line
      call uwsgxa ((1.0e-3)*axy ,  ny) ! yaxis value
      call uwsgya ((1.0e-3)*zlev,nlev) ! zaxis value (km)

!      call udgcla (-0.3, 0.3, 0.03)
      call udcntr (var_yz,ny,ny,nlev)

      call udiclv
      call udiset ('INDXMJ',23) ! major line
      call udiset ('INDXMN',21) ! minor line
      call udcntr (var_yz_2,ny,ny,nlev)
! select case(trim(cvarc))
! case('u');  call udcntr (u_snap,nx,nx,nlev)
! case('w');  call udcntr (w_snap,nx,nx,nlev)
! end select

! *** vector ***

!select case (trim(cvarv))
! case ('uw'); call ugvect (u_snap,nx,w_snap/zasp,nx,nx,nlev)
!end select

! *** lines on xy plane ***

!  call uuslnt(1) ! line type
!  call uuslni(5) ! line index (color*10 + index*1)
!  do ilev=1,nlev
!   select case (trim(cvarl))
!   case ('xz');call uulin (nx,0.001*(axx(:)+dispx_snap(:,ilev)),0.001*(zlev(ile!v)+dispz_snap(:,ilev)))
!   case ('z'); call uulin (nx,0.001*(axx(:)),0.001*(zlev(ilev)+dispz_snap(:,ile!v)))
!   end select
!  end do

! **** y, z axis ****
      call uzinit
      call uziset ('INDEXT2',3)
      call uziset ('INDEXL1',5)
      call uzrset ('RSIZEL1',0.018)
      call uzrset ('RSIZEC1',0.018)
      call uzrset ('RSIZET2',0.010)
      call uzrset ('RSIZET1',0.004)
!      call uzlset ('LABELXB',.FALSE.)
!      call uzlset ('LOFFSET',.TRUE.)
 !     call uzrset ('XFACT',0.001)
      call uzlset ('LOFFSET',.FALSE.)
      call uxaxdv ('B',ystics,ymtics)
      call uxaxdv ('T',ystics,ymtics)
      call uxsttl ('B','y(km)',0.0)
      call uysfmt ('(I2)')
      call uzlset ('LABELYR',.FALSE.)
      call uzlset ('LABELYL',.TRUE.)
      call uyaxdv ('L',1.0,5.0)
      call uyaxdv ('R',1.0,5.0)
      call uziset ('IROTCYL',1)
      call uysttl ('L','Height(km)',0.0)

      call sgtxzv (0.50,0.76,trim(title1),0.024,0,0,5) ! title
      call sgtxzv (0.17,0.76,trim(title2),0.020,0,-1,5) ! title
      call grcls
end do

return
end subroutine quicklook_yz
!======================================!

subroutine quicklook_xy(title)
use setup

character(*)::title
character*40::title1,title2

namelist /quicklook_nml/ ylq, yrq, xlq, xrq, ylatlq, ylatrq, xlonlq, xlonrq
namelist /intp_nml/ axlon_base, axlon_cut_l, axlon_cut_r, axlat_cut_l, axlat_cut_r, xtaper, ytaper
namelist /contour_nml/ ifree, vmax1, vmin1, vmax2, vmin2

ylq=yl
yrq=yr
xlq=xl
xrq=xr

open(11,file='calc.conf')
 read(11,nml=quicklook_nml)
 read(11,nml=contour_nml)
close(11)

if (ifile.eq.1)then
 open(11,file='calc.conf')
  read(11,nml=intp_nml)
 close(11)
 ylq=er*drad*ylatlq
 yrq=er*drad*ylatrq
 xlq=er*drad*(xlonlq-axlon_base)
 xrq=er*drad*(xlonrq-axlon_base)
end if

iylq=max(iblkge(axy,ny,ylq),1)
iyrq=min(iblkle(axy,ny,yrq),ny)
ixlq=max(iblkge(axx,nx,xlq),1)
ixrq=min(iblkle(axx,nx,xrq),nx)

!write(*,*)iylq,iyrq,ixlq,ixrq

ixlen=ixrq-ixlq+1
iylen=iyrq-iylq+1


title1=title

xstics=5.0
xmtics=10.0
ystics=2.5
ymtics=5.0

! output
iouts=1
ioute=2
if(iout.eq.1) ioute=1
if(iout.eq.2) iouts=2
do ioutl=iouts,ioute
! *** general settings ***
      call sgiset ('IFONT',1)
      call swistx ('ICLRMAP',14) ! colormap blue-white-red 
      call swcmll
      call swcset ('FNAME','figure')
      call swiset ('IWIDTH',1000) ! output window size
      call swiset ('IHEIGHT',800)     
      call gropn(ioutl) 
      call sglset ('LFULL',.TRUE.) ! using fullsize
      call swlset ('LSEP',.TRUE.) ! psfile numbering
      call slmgn (0.0,0.0,0.0,0.0) ! margin 
      call grfrm 
      call grswnd (1.0e-3*xlq,1.0e-3*xrq,1.0e-3*ylq,1.0e-3*yrq) ! set window
      call grsvpt (0.15,0.85,0.23,0.63) ! set viewport !!! not 1:1
      call grstrn (1) ! linear or log
      call grstrf


  yasp= (yr-yl) / (xr-xl) 

! **** contour ****

      call sglset ('LCLIP',.TRUE.) ! clipping
      call udlset ('LMSG',.FALSE.) ! 'contour interval'
      call udrset ('RSIZEL',0.015) ! label
      call udiset ('INDXMJ',13) ! major line
      call udiset ('INDXMN',11) ! minor line
      call uwsgxa ((1.0e-3)*axx(ixlq:ixrq),ixlen) ! xaxis value
      call uwsgya ((1.0e-3)*axy(iylq:iyrq),iylen) ! yaxis value

      if (ifree.eq.0) call udgcla (vmin1, vmax1, -20.)
      if (ifree.eq.0) call uddclv (0.0)
      call udcntr (var_xy(ixlq:ixrq,iylq:iyrq),ixlen,ixlen,iylen)

      call udiset ('INDXMJ',23) ! major line
      call udiset ('INDXMN',21) ! minor line
      if (ifree.eq.0) call udgcla (vmin2, vmax2, -20.)
      if (ifree.eq.0) call uddclv (0.0)
      call udcntr (var_xy_2(ixlq:ixrq,iylq:iyrq),ixlen,ixlen,iylen)

! select case(trim(cvarc))
! case('u');  call udcntr (u_snap,nx,nx,nlev)
! case('w');  call udcntr (w_snap,nx,nx,nlev)
! end select

! *** vector ***

!select case (trim(cvarv))
! case ('uw'); call ugvect (u_snap,nx,w_snap/zasp,nx,nx,nlev)
!end select

! *** lines on xy plane ***

!  call uuslnt(1) ! line type
!  call uuslni(5) ! line index (color*10 + index*1)
!  do ilev=1,nlev
!   select case (trim(cvarl))
!   case ('xz');call uulin (nx,0.001*(axx(:)+dispx_snap(:,ilev)),0.001*(zlev(ile!v)+dispz_snap(:,ilev)))
!   case ('z'); call uulin (nx,0.001*(axx(:)),0.001*(zlev(ilev)+dispz_snap(:,ile!v)))
!   end select
!  end do

! **** y, z axis ****
      call uzinit
      call uziset ('INDEXT2',3)
      call uziset ('INDEXL1',5)
      call uzrset ('RSIZEL1',0.018)
      call uzrset ('RSIZEC1',0.018)
      call uzrset ('RSIZET2',0.010)
      call uzrset ('RSIZET1',0.004)
!      call uzlset ('LABELXB',.FALSE.)

if (ifile.eq.1)then
   call uzlset ('LOFFSET',.TRUE.)
   call uzrset ('XFACT',1.0e3/er/drad)
   call uzrset ('XOFFSET',axlon_base)
   call uzrset ('YFACT',1.0e3/er/drad)
   xstics=10.0
   xmtics=30.0
   ystics=5.0
   ymtics=10.0
else
    call uzlset ('LOFFSET',.TRUE.)
    call uzrset ('XFACT',1.0e-3)
    call uzrset ('YFACT',1.0e-3)
end if
      call uysfmt ('(I4)')
      call uxaxdv ('B',xstics,xmtics)
      call uxaxdv ('T',xstics,xmtics)
      call uzlset ('LABELYR',.FALSE.)
      call uzlset ('LABELYL',.TRUE.)
      call uyaxdv ('L',ystics,ymtics)
      call uyaxdv ('R',ystics,ymtics)


if (ifile.eq.1)then
   call uxsttl ('B','Longitude',0.0)
   call uysttl ('L','Latitude',0.0)
else
   call uxsttl ('B','x (10|3"km)',0.0)
   call uysttl ('L','y(10|3"km)',0.0)
end if


      call sglset ('LCLIP',.FALSE.) ! clipping
      call sgtxzv (0.50,0.66,trim(title1),0.024,0,0,5) ! title
      call sgtxzv (0.17,0.66,trim(title2),0.020,0,-1,5) ! title
      call grcls
end do

return
end subroutine quicklook_xy
!======================================!
