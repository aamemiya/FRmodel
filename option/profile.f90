program profile
use EZspline_obj
use EZspline


include '../alpha/Pcon.h'
include '../alpha/Grid.h'

integer,parameter::nrmax=300
character*30::fname='uvt_profile.dat'

real(4),parameter::p0=1.0e5 !!! surface pressure

real(4)::pread(nrmax)
real(4)::uread(nrmax)
real(4)::vread(nrmax)
real(4)::tread(nrmax)

real(4)::vm2smp(nrmax)

real(4)::zread(nrmax)
real(4)::bn2read(nrmax)

real(4)::ubg(nlev)
real(4)::vbg(nlev)
real(4)::bn2(nlev)
real(4)::tbg(nlev)

real(4),parameter::dmiss=-9.99e10

type(EZspline1_r4)::u_spl,v_spl,t_spl,bn2_spl

pread=dmiss
uread=dmiss
vread=dmiss
tread=dmiss

open(11,file=trim(fname),form='formatted')
do ir=1,nrmax
 read(11,*,end=100) pread(ir),uread(ir),vread(ir),tread(ir)
end do
100 close(11)

nr=count(pread.ne.dmiss)

pread=pread*1.0e5 !!! hPa ==> Pa

!!! ( pread(1:iupb-1) > p0, pread(iupb:nr) <= p0 )
iupb=nr-iblkge(pread(nr:1:-1),nr,p0)+1

!!! calculate z
zread(iupb) = 0.0 + rd/grav*0.5*(tread(iupb)) * log(p0/pread(iupb))
do iz=iupb+1,nr
 zread(iz)=zread(iz-1)+rd/grav*0.5*(tread(iz)+tread(iz-1)) * log(pread(iz-1)/pread(iz))
end do

if (iupb.gt.2)then
zread(iupb-1) = 0.0 + rd/grav*0.5*(tread(iupb-1)) * log(p0/pread(iupb))
do iz=iupb-2,1
 zread(iz) = zread(iz+1) - rd/grav*0.5*(tread(iz)+tread(iz+1)) * log(pread(iz)/pread(iz+1))
end do
end if

bn2read(1) = grav/tread(1) * (grav/cp - (tread(2)-tread(1)) / (zread(2)-zread(1)) )
bn2read(nr) = grav/tread(nr) * (grav/cp - (tread(nr)-tread(nr-1)) / (zread(nr)-zread(nr-1)) )

bn2read(2:nr-1)=2.0*grav/(tread(1:nr-2)+tread(3:nr)) *(grav/cp - (tread(3:nr)-tread(1:nr-2)) / (zread(3:nr)-zread(1:nr-2)) )



!uread=uread*3.0

!do lenx=50,500,50

!vk0=2.0*pi/ ( 30.0* 1.0e3)
!vk0=2.0*pi/ ( real(lenx)* 1.0e3)
!f0=2.0*eomg*sin(-60.0*drad)

!vm2smp= vk0**2 * (bn2read - (uread*vk0)**2) / ((uread*vk0)**2-f0**2) !!! disp. relation

!write(*,*) 'lenx=',lenx
!do ir=1,nr
! write(*,*) ir, vm2smp(ir), uread(ir), pread(ir)
!end do
!end  do
!stop

!!! intp

if (ztop .gt. zread(nr))then
 write(*,*) 'error:: ztop > max(zread)'
 stop
end if

call EZspline_init(u_spl,nr,(/0,0/),ier)
call EZspline_init(v_spl,nr,(/0,0/),ier)
call EZspline_init(t_spl,nr,(/0,0/),ier)
call EZspline_init(bn2_spl,nr,(/0,0/),ier)
u_spl%x1=zread(1:nr)
v_spl%x1=zread(1:nr)
t_spl%x1=zread(1:nr)
bn2_spl%x1=zread(1:nr)
call EZspline_setup(u_spl, uread(1:nr), ier)
call EZspline_setup(v_spl, vread(1:nr), ier)
call EZspline_setup(t_spl, tread(1:nr), ier)
call EZspline_setup(bn2_spl, bn2read(1:nr), ier)
call EZspline_interp(u_spl,nlev,zlev,ubg,ier)
call EZspline_interp(v_spl,nlev,zlev,vbg,ier)
call EZspline_interp(t_spl,nlev,zlev,tbg,ier)
call EZspline_interp(bn2_spl,nlev,zlev,bn2,ier)
call EZspline_free(u_spl,ier)
call EZspline_free(v_spl,ier)
call EZspline_free(t_spl,ier)
call EZspline_free(bn2_spl,ier)




 !!! surface wind direction

 if (ubg(1)**2+vbg(1)**2.eq.0.0)then 
  vudir=0.0
 else
  vue   = ubg(1)/sqrt(ubg(1)**2+vbg(1)**2)
  vve   = vbg(1)/sqrt(ubg(1)**2+vbg(1)**2)
  vudir = acos(vue)

  if (vve.lt.0.0)vudir = 360.0-vudir
 end  if

 ubg=cos(vudir*drad)*ubg+sin(vudir*drad)*vbg
 vbg=cos(vudir*drad)*vbg-sin(vudir*drad)*ubg

 !!! 


ibnd=maxval(maxloc(bn2(2:nlev)-bn2(1:nlev-1)))

zbnd=0.5*(zlev(ibnd)+zlev(ibnd+1))

!vmin=minval(tread,tread.ne.dmiss)
!vmax=maxval(tread,tread.ne.dmiss)

!vmin=minval(uread,uread.ne.dmiss)
!vmax=maxval(uread,uread.ne.dmiss)

vmin=minval(bn2read,bn2read.ne.dmiss)
vmax=maxval(bn2read,bn2read.ne.dmiss)



!vmin=minval(vm2smp,vm2smp.gt.-1.0e-6)
!vmax=maxval(vm2smp,vm2smp.lt.1.0e-5)

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
      call swlset ('LSEP',.FALSE.) ! psfile numbering
      call slmgn (0.0,0.0,0.0,0.0) ! margin 
      call grfrm 
      call grsvpt (0.17,0.77,0.13,0.73) ! set new window
      call grswnd (vmin,vmax,zbot,ztop) ! set viewport
      call grstrn (1) ! linear or log
      call grstrf

      call uuslnt (1)
      call uuslni (2)

      call sglset ('LCLIP',.TRUE.) ! using fullsize
!      call uulin  (nr,tread(1:nr),zread(1:nr))
!      call uulin  (nr,uread(1:nr),zread(1:nr))
      call uulin  (nr,bn2read(1:nr),zread(1:nr))
      call uuslnt (3)
      call uulin  (nlev,bn2,zlev)
!      call uulin  (nr,vm2smp(1:nr),zread(1:nr))

      call uulin  (2,(/vmin,vmax/),(/zbnd,zbnd/))

      call uzinit
      call uziset ('INDEXT2',3)
      call uziset ('INDEXL1',5)
      call uzrset ('RSIZEL1',0.014)
      call uzrset ('RSIZEC1',0.014)
      call uzrset ('RSIZET2',0.010)
      call uzrset ('RSIZET1',0.004)
      call uxaxdv ('B',10.0,50.0)
      call uxsttl ('B','T(K)',0.0)

      call uzlset ('LOFFSET',.true.)
      call uzrset ('YOFFSET',0.0)
      call uzrset ('YFACT',1.e-3)
      call uyaxdv ('L',10.0,20.0)
      call uysttl ('L','z (km)',0.0)
      call uysfmt ('(I2)')
      call uzlset ('LABELYR',.FALSE.)
      call uzlset ('LABELYL',.TRUE.)

      call sglset ('LCLIP',.FALSE.) ! using fullsize
      call sgtxzv (0.45,0.75,'sample profile',0.020,0,0,5) ! title
      call sgtxzv (0.17,0.75,'',0.016,0,-1,5) ! title
      call grcls

end do


stop
end program profile
