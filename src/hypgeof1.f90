!-------------------------------------------------------------------
! function  : f1 (complex(8))
!
! package   : F1
!
! Language  : Fortran 90
!
! author    : F. Colavecchia (flavioc@lanl.gov)
!                            (flavioc@cab.cnea.gov.ar)
!             Daniel Sabanes Bove (daniel.sabanesbove@ifspm.uzh.ch)
!             
!
! date      : 3/26/97      version: 0.1
! revision  : 6/25/02      version: 1.0
! modified  : 3/28/2012    DSB: use R warning message
!
! purpose   :  Computes the F1 function following the ideas in
!              the paper CPC 138 (2001) 29, eq. (33).
!         		 Calculates the Appell's hypergeometric function 
!			         f1(a,b,b',c,x,y) in the convergence region
!			         through the numerical integration of ODE 
!			         equation representing the PDE system of 
!			         F1.
!
! input     :  a  -> complex parameter of Appell's F1  
!              b  -> complex parameter of Appell's F1  
!              bp -> complex parameter of Appell's F1  
!              c  -> complex parameter of Appell's F1  
!              u  -> complex variable                  
!              v  -> complex variable                  
!              val -> complex function value
!              iflag -> integer integration flag
!              hyp2f1 -> 1 for Forrey and 2 for Michel & Stoitsov
!                       hypergeometric function code
!
! modules   :    uses Felhberg fifth-four Runge-Kutta method for
!                integration
!
!
!-------------------------------------------------------------------
subroutine hypgeof1(a,b,bp,c,u,v,val,iflag,hyp2f1)
  implicit none
  complex(8) val,a,b,bp,c,u,v,f21
  double precision eps,fact,zero,abserr
  parameter (eps=1.d-8)
  parameter (fact=.1)
  parameter (zero=1e-12)
  parameter (abserr=1d-40)
  integer kmax,neqn
  integer iwork(5)
  integer iflag
  integer hyp2f1
  external bsstep,hypdrvf1,rwarn
  real(8) au,av,t0,one,su,sv
  real(8) start,finish
  real(8) work(50)
  complex(8) z0,dz,aa,bb,bbp,cc,y(3),xx,yy,z,c1
  common /hypgf1/ aa,bb,bbp,cc,xx,yy,z0,dz
  common /path/ kmax
  t0=0d0
  kmax=0
  one = 1d0
  c1 = (1d0,0d0)
  au = cdabs(u)
  av = cdabs(v)
  if(au.lt.zero.and.av.lt.zero) then
     val = c1
     return
  end if
  !	
  !	  There still exists a divergence behavior when u=v 
  !     even for |u,v| < 1, so we skip this case
  !
  if(cdabs(u-v).lt.zero) then
     val = f21(a,b+bp,c,u,hyp2f1)
     return
  endif

  su = dreal(u)/au
  sv = dreal(v)/av
  !
  !	  Select z0, starting integration point
  !
  if(su.lt.0d0.and.dreal(v).lt.one.or.sv.lt.0d0.and.dreal(u).lt.one) then
     t0 = 1d0/(16d0*max(au,av))
  elseif(au.lt.one.and.av.lt.one) then
     t0 = 1d0/(5d0*max(au,av))
  elseif(au.gt.one.and.av.lt.one) then
     t0 = (1d0+fact)/au
  elseif(av.gt.one.and.au.lt.one) then
     t0 = (1d0+fact)/av
  endif
  z0 = dcmplx(t0,0d0)
  z  = dcmplx(1.d0,0.d0)
  aa=a
  bb=b
  bbp=bp
  cc=c
  xx=u
  yy=v
  dz=z-z0
  !
  !	  Initial values of the function f1
  !
  call hypstartf1(a,b,bp,c,u,v,t0*c1,y,hyp2f1)

  neqn = 6
  iflag = 1
  start = 0d0
  finish = 1d0
  !
  !   Proceed with the integration
  !
  call rkf45(hypdrvf1, neqn, y, start,finish, EPS, abserr, iflag, work, iwork )
  val=y(1)

  if(iflag.ne.2) then
     ! write(*,*) 'iflag =',iflag 
     call rwarn('hypgeof1: Problems with the integration')
  end if
  return
end subroutine hypgeof1
