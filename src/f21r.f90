!----------------------------------------------------------------------
!     function f21						3/2/98
!     ____________
!
!     by: f. colavecchia  (flavioc@cab.cnea.gov.ar)
!    
!     additions by: Daniel Sabanes Bove (daniel.sabanesbove@ifspm.uzh.ch)  
!     (use R error message and additional input parameter)
!
!     wrapper for the hypergeometric function 
!     2f1
!
!     input:	
!     complex(8) a,b,c,z
!     integer hyp2f1 -> 1 for Forrey and 2 for Michel & Stoitsov
!                       hypergeometric function code
!     
!     output: complex(8) 2f1(a,b,c,z)
!     
!
!---------------------------------------------------------------------
!

complex(8) function f21(a,b,c,z,hyp2f1)
  implicit none

  integer(4) isneg
  integer hyp2f1
  real(8) er,zero
  complex(8) cf,a,b,c,z,cgamma,cgammar
  complex(8) ci,chyp
  external rexit

  ci = (0d0,1d0)
  f21 = (1d0,0d0)
  er = 1.0e-9
  zero = 1e-5*er
  !
  if(cdabs(z).lt.zero) return
  if(cdabs(a-c).lt.zero) then
     f21 = (1-z)**(-b)
     return
  end if
  if(cdabs(b-c).lt.zero) then
     f21 = (1-z)**(-a)
     return
  end if
  !	discard the case z=1
  if(cdabs(z-1d0).lt.zero) then
     f21 = cgamma(c)*cgamma(c-a-b)*cgammar(c-a)*cgammar(c-b)
     return
  end if
  if(cdabs(a).le.zero.or.cdabs(b).le.zero) then
     f21 = (1d0,0d0)
     return
  else if(isneg(c).eq.-1) then
     call rexit('f21: c is a negative integer')
  end if

  !	Decide which algorithm to use for computing the
  !       the Gaussian hypergeometric function
  if(hyp2f1.eq.1) then
     !	Forrey's F21
     cf = chyp(a,b,c,z)
  else if(hyp2f1.eq.2) then
     !  Michel & Stoitsov
     call HYP_2F1_SUBROUTINE(a,b,c,z,cf)
  else 
     !  unknown choice
     call rexit('f21: hyp2f1 must be either 1 or 2')
  end if

  f21 = cf

  return
end function f21


