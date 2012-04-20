!-------------------------------------------------------------------
! subroutine: hypstartf1
!
! package   : F1
!
! Language  : Fortran 90
!
! author    : F. Colavecchia (flavioc@lanl.gov)
!                            (flavioc@cab.cnea.gov.ar)
!
! date      : 3/26/97      version: 0.1
! revision  : 6/25/02      version: 1.0
!
! purpose   :  Computation of initial values of the function
!	          	 and its derivatives to be used in the ODE integration of
!              the system (33), see paper CPC 138 (2001) 29.
!              We make use of the series F1 as in eq. (31)
!
! input     :    a  -> complex parameter of Appell's F1
!                b  -> complex parameter of Appell's F1
!                bp -> complex parameter of Appell's F1
!                c  -> complex parameter of Appell's F1
!                u  -> complex variable
!                v  -> complex variable
!                t0  -> initial step value
!                hyp2f1 -> 1 for Forrey and 2 for Michel & Stoitsov
!                       hypergeometric function code
!
!
! output    :    y(3)  -> F1 Appell's hypergeometric function
!                         first and second derivatives at (u*t0,v*t0) 
!
!
!-------------------------------------------------------------------

subroutine hypstartf1(a,b,bp,c,u,v,t0,y,hyp2f1)
  implicit none
  complex(8) au,av,t0
  complex(8) a,b,bp,c,u,v,y(3)
  complex(8) caux21,caux22,caux31,caux32,caux33,c1
  complex(8) cf1bnl
  integer hyp2f1

  c1 = (1d0,0d0)
  au = u*t0
  av = v*t0
  caux21 = cf1bnl(1 + a,b,1 + bp,1 + c,au,av,hyp2f1)
  caux22 = cf1bnl(1 + a,1 + b,bp,1 + c,au,av,hyp2f1)
  caux31 = cf1bnl(2 + a,b,2 + bp,2 + c,au,av,hyp2f1)
  caux32 = cf1bnl(2 + a,1 + b,1 + bp,2 + c,au,av,hyp2f1)
  caux33 = cf1bnl(2 + a,2 + b,bp,2 + c,au,av,hyp2f1)
  y(1) = cf1bnl(a,b,bp,c,au,av,hyp2f1)
  y(2) = (a*(bp*v*caux21+ b*u*caux22))/c
  y(3) = (a*(1+a)*(bp*(1+bp)*v**2*caux31 + b*u*(2*bp*v*caux32+(1+b)*u*caux33)))/(c*(1+c))
  return
end subroutine hypstartf1

