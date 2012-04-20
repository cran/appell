!-------------------------------------------------------------------
! function  : g2 (complex(8))
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
! purpose   :  Computes the G2 function, Eq. (34) of paper CPC 138 (2001) 29.
!
!
! input     :    a1  -> complex parameter of Appell's F2
!                a2 -> complex parameter of Appell's F2
!                b1 -> complex parameter of Appell's F2
!                b2 -> complex parameter of Appell's F2
!                x  -> real variable
!                y  -> real variable
!                hyp2f1 -> 1 for Forrey and 2 for Michel & Stoitsov
!                       hypergeometric function code
!
! output    :    g2  ->G2 function
!
!-------------------------------------------------------------------
complex(8) function g2(a1,a2,b1,b2,x,y,hyp2f1)

  real(8) x,y  
  complex(8) a1,a2,b1,b2,f2
  complex(8) cx,cy
  logical debug
  integer hyp2f1

  debug = .false.
  cx = x*(1d0,0d0)
  cy = y*(1d0,0d0)
  if (debug) write(*,*) " Computing g2"
  if (debug) write(*,*) " x=",x
  if (debug) write(*,*) " y=",y

  if (debug) write(*,*) a1," ",a2," ",b1," ",b2
  g2 =(1+cx)**(-a1)*(1+cy)**(-a2)*f2(1-b1-b2,a1,a2,1-b1,1-b2,cx/(cx+1),cy/(cy+1),hyp2f1)
  if (debug) write(*,*) g2
  if (debug) write(*,*) " End computing g2"
  return
end function g2

!-------------------------------------------------------------------
! function  : f2 (complex(8))
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
! purpose   :  Computes the  F2 series 
!              in the convergence region.
!
!
! input     :    a  -> complex parameter of Appell's F2
!                b1 -> complex parameter of Appell's F2
!                b2 -> complex parameter of Appell's F2
!                c1 -> complex parameter of Appell's F2
!                c2 -> complex parameter of Appell's F2
!                x  -> real variable
!                y  -> real variable
!                hyp2f1 -> 1 for Forrey and 2 for Michel & Stoitsov
!                       hypergeometric function code
!
! output    :    f2  -> nth sum of F2
!
!-------------------------------------------------------------------
complex(8) function f2(a,b1,b2,c1,c2,cx,cy,hyp2f1)
  implicit none 
  logical debug, ispossible
  complex(8) a,b1,b2,c1,c2,f21,s,f2s,c0,cx,cy
  real(8) tmax1,tmax2
  integer flag
  integer hyp2f1

  ispossible =.false.
  debug = .false.

  c0 = dcmplx(0.,0.)
  if(b1.eq.c1) then
     if (debug) write(*,*) "b1=c1 " 
     f2 = (1-cx)**(-a)*f21(a,b2,c2,cy/(1-cx),hyp2f1)
     return
  else if(b2.eq.c2) then  
     if (debug) write(*,*) "b2=c2 " 
     f2 =  (1-cy)**(-a)*f21(a,b1,c1,cx/(1-cy),hyp2f1)
     return
  else if(b1.eq.0) then
     if (debug) write(*,*) "b1=0  " 
     f2 = f21(a,b2,c2,cy,hyp2f1)
     return
  else if(b2.eq.0) then 
     if (debug) write(*,*) "b2=0  " 
     f2 =  f21(a,b1,c1,cx,hyp2f1)
     return
  else if(cx.eq.0) then 
     f2 =  f21(a,b2,c2,cy,hyp2f1)
     return
  else if(cy.eq.0) then
     f2 = f21(a,b1,c1,cx,hyp2f1)
     return
  end if

  flag=  0 
  tmax1=  1.0 
  tmax2=real(cdsqrt(cx**2+cy**2))
  if (debug) write(*,*) flag," ", tmax1,"   ",tmax2
  if(tmax2.lt.tmax1) then
     flag=  1 
     ispossible = .true. 
     tmax1=tmax2
  end if
  if (debug) write(*,*) flag," ", tmax1,"   ",tmax2

  tmax2=  real(cdsqrt((cx/(1-cy))**2+(cy/(cy-1))**2)) 
  if(tmax2.lt.tmax1) then
     flag=  22
     ispossible = .true. 
     tmax1=tmax2
  end if
  if (debug) write(*,*) flag," ", tmax1,"   ",tmax2 

  tmax2=  real(cdsqrt((cx/(cx-1))**2+(cy/(1-cx))**2)) 
  if(tmax2.lt.tmax1) then
     flag=  21 
     ispossible = .true. 
     tmax1=tmax2
  end if

  if (debug) write(*,*) flag," ", tmax1,"   ",tmax2 
  tmax2=  real(cdsqrt((cx/(cx+cy-1))**2+(cy/(cx+cy-1))**2))
  if(tmax2.lt.tmax1) then
     flag=  3
     ispossible = .true. 
     tmax1=tmax2
  end if

  if (debug) write(*,*) flag," ", tmax1,"   ",tmax2
  if(flag.eq.1) then
     if (debug) write(*,*) "Series in f21(cx)     :",cx," ",cy
     s=f2s(a,b1,b2,c1,c2,cx,cy,hyp2f1)
  else if(flag.eq.21) then
     if (debug) write(*,*) "Transformation 2x:",cx/(cx-1)," ",cy/(1-cx)
     s=(1-cx)**(-a)*f2s(a,c1-b1,b2,c1,c2,cx/(cx-1),cy/(1-cx),hyp2f1)
  else if(flag.eq.22) then  
     if (debug) write(*,*) "Transformation 2y:",cx/(1-cy)," ",cx/(1-cy) 
     s=(1-cy)**(-a)*f2s(a,b1,c2-b2,c1,c2,cx/(1-cy),cy/(cy-1),hyp2f1)
  else if(flag.eq.3) then
     if (debug) write(*,*) "Transformation 3: ",cx/(cx+cy-1)," ",cy/(cx+cy-1) 
     s=(1-cx-cy)**(-a)*f2s(a,c1-b1,c2-b2,c1,c2,cx/(cx+cy-1),cy/(cx+cy-1),hyp2f1)
  else  
     if (debug) write(*,*) "Not Possible"
     f2 = c0
     return
  end if
  f2 = s
  return
end function f2
!-------------------------------------------------------------------
! function  : f2s (complex(8))
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
! purpose   :  Computes the  F2 series 
!              in the convergence region.
!
!
! input     :    a  -> complex parameter of Appell's F2
!                b1 -> complex parameter of Appell's F2
!                b2 -> complex parameter of Appell's F2
!                c1 -> complex parameter of Appell's F2
!                c2 -> complex parameter of Appell's F2
!                cx  -> complex variable
!                cy  -> complex variable
!                hyp2f1 -> 1 for Forrey and 2 for Michel & Stoitsov
!                       hypergeometric function code
!
! output    :    f2s  -> sum of F2 series
!
!-------------------------------------------------------------------
complex(8) function f2s(a,b1,b2,c1,c2,cx,cy,hyp2f1)
  implicit none 
  complex(8) a,b1,b2,c1,c2,cx,cy,f2sx
  integer hyp2f1

  if(cdabs(cy).lt.cdabs(cx)) then
     f2s = f2sx(a,b2,b1,c2,c1,cy,cx,hyp2f1)
  else
     f2s = f2sx(a,b1,b2,c1,c2,cx,cy,hyp2f1)
  end if
  return
end function f2s
!-------------------------------------------------------------------
! function  : f2sx (complex(8))
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
! purpose   :  Computes the  F2 series 
!              in the convergence region, single index,
!              single Gauss function.
!
!
! input     :    a  -> complex parameter of Appell's F2
!                b1 -> complex parameter of Appell's F2
!                b2 -> complex parameter of Appell's F2
!                c1 -> complex parameter of Appell's F2
!                c2 -> complex parameter of Appell's F2
!                x  -> real variable
!                y  -> real variable
!                hyp2f1 -> 1 for Forrey and 2 for Michel & Stoitsov
!                       hypergeometric function code
!
! output    :    f2sx  -> converged sum F2
!
! notes     :  most of the code is f77, few things come from f90.
!
!-------------------------------------------------------------------
complex(8) function f2sx(a,b1,b2,c1,c2,cx,cy,hyp2f1)
  implicit none 
  integer(4) i
  complex(8) a,b1,b2,c1,c2,cx,cy
  complex(8) coef,f21,suma,tmp,c
  real(8) rtest
  integer hyp2f1

  ! if (debug) write(*,*) a," ",b1," ",b2," ",c1," ",c2
  !if (debug) write(*,*) "Calculando f2 con serie en f21 :",x," ",y

  suma = f21(a,b1,c1,cx,hyp2f1) 
  tmp = suma
  coef = (1d0,0d0)
  c = (0d0,0d0)
  i=0 
  rtest =real(tmp/suma)
  do while(i.lt.300 .and. (rtest.gt.1e-5))
     i=i+1 
     coef = coef*(a+i-1)*(b2+i-1)/((c2+i-1)*i)
     tmp  = coef*f21(a+i,b1,c1,cx,hyp2f1)*cy**(i)
     suma=suma+tmp
     rtest =real(cdabs(tmp/suma))
  end do
  if(i.gt.200) then 
     !	if (debug) write(*,*) 'Max. number of terms obtained'
  end if
  f2sx = suma
  !	if (debug) write(*,*) i
  return
end function f2sx


!-------------------------------------------------------------------	
!
! Next functions are included for debugging purposes.
!
!-------------------------------------------------------------------
! function  : coef2 (complex(8))
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
! purpose   :  Computes the n-th coefficient of the F2 series 
!              of Burchnall and Chaundy, Eq. (36) of paper CPC 138 (2001) 29.
!
!
! input     :    a  -> complex parameter of Appell's F2
!                b1 -> complex parameter of Appell's F2
!                b2 -> complex parameter of Appell's F2
!                c1 -> complex parameter of Appell's F2
!                c2 -> complex parameter of Appell's F2
!                n  -> order
!
! output    :    coef2  -> coefficient
!
! notes     :  most of the code is f77, few things come from f90.
!
!-------------------------------------------------------------------
complex(8) function coef2(a,b1,b2,c1,c2,n)
  implicit none 
  complex(8) a,b1,b2,c1,c2,pochhammer
  real(8) fact
  integer(4) n
  coef2 = pochhammer(a,n)*pochhammer(b1,n)*pochhammer(b2,n)/(pochhammer(c1,n)*pochhammer(c2,n)*fact(n))
  return
end function coef2
!-------------------------------------------------------------------
! function  : f2term (complex(8))
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
! purpose   :  Computes the n-th term of the F2 series 
!              of Burchnall and Chaundy, Eq. (36) of paper CPC 138 (2001) 29.
!
!
! input     :    a  -> complex parameter of Appell's F2
!                b1 -> complex parameter of Appell's F2
!                b2 -> complex parameter of Appell's F2
!                c1 -> complex parameter of Appell's F2
!                c2 -> complex parameter of Appell's F2
!                x  -> real variable
!                y  -> real variable
!                n  -> order
!                hyp2f1 -> 1 for Forrey and 2 for Michel & Stoitsov
!                       hypergeometric function code
!
! output    :    coef2  -> nth term of F2
!
! notes     :  most of the code is f77, few things come from f90.
!
!-------------------------------------------------------------------
complex(8) function f2term(a,b1,b2,c1,c2,x,y,n,hyp2f1)
  implicit none 
  complex(8) a,b1,b2,c1,c2,coef2,f21
  real(8) x,y 
  complex(8) cx,cy
  integer(4) n
  integer hyp2f1

  cx = x*(1d0,0d0)
  cy = y*(1d0,0d0)
  f2term = coef2(a,b1,b2,c1,c2,n)*(x*y)**n*f21(a+n,b1+n,c1+n,cx,hyp2f1)*f21(a+n,b2+n,c2+n,cy,hyp2f1)
  return
end function f2term
!-------------------------------------------------------------------
! function  : f2n (complex(8))
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
! purpose   :  Computes the n-th sum of the F2 series 
!              of Burchnall and Chaundy, Eq. (36) of paper CPC 138 (2001) 29.
!
!
! input     :    a  -> complex parameter of Appell's F2
!                b1 -> complex parameter of Appell's F2
!                b2 -> complex parameter of Appell's F2
!                c1 -> complex parameter of Appell's F2
!                c2 -> complex parameter of Appell's F2
!                x  -> real variable
!                y  -> real variable
!                n  -> order
!                hyp2f1 -> 1 for Forrey and 2 for Michel & Stoitsov
!                       hypergeometric function code
!
! output    :    f2n  -> nth sum of F2
!
! notes     :  most of the code is f77, few things come from f90.
!
!-------------------------------------------------------------------
complex(8) function f2n(a,b1,b2,c1,c2,x,y,N,hyp2f1)
  implicit none 
  complex(8) a,b1,b2,c1,c2,f2term
  real(8) x,y 
  integer(4) i,n
  integer hyp2f1
  f2n = cmplx(0d0,0d0)
  do i=0,n
     f2n = f2n + f2term(a,b1,b2,c1,c2,x,y,i,hyp2f1)
  end do
  return 
end function f2n
