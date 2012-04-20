!-------------------------------------------------------------------
! subroutine: f1
!
! package   : F1
!
! Language  : Fortran 90
!
! author    : F. Colavecchia (flavioc@lanl.gov)
!                            (flavioc@cab.cnea.gov.ar)
!             Daniel Sabanes Bove (daniel.sabanesbove@ifspm.uzh.ch)               
!
! date      : 3/26/97      version: 0.1
! revision  : 6/25/02      version: 1.0
! revision  : 7/25/02      version: 1.0.1  Patch eq. (23) and (24)
!                                          check parameters included.
!             3/26/2012    modified into a subroutine by DSB,
!                          several other modifications
!                          
!             
! purpose   :  Computes the F1 function following the ideas in 
!              the paper CPC 138 (2001) 29.
!
! input     :    a  -> complex parameter of Appell's F1
!                b1 -> complex parameter of Appell's F1
!                b2 -> complex parameter of Appell's F1
!                c  -> complex parameter of Appell's F1
!                x  -> real variable
!                y  -> real variable
!                algoflag -> integer flag variable, determined
!                            by the algorithm
!                userflag -> the user can supply a flag.
!                            negative value means no user flag.
!                debug -> debug logical
!                val -> complex return value
!                hyp2f1 -> 1 for Forrey and 2 for Michel & Stoitsov
!                       hypergeometric function code
!
! output    :    f1  -> F1 Appell's hypergeometric function
!
!
! notes     :  most of the code is f77, few things come from f90.
!
!-------------------------------------------------------------------

subroutine f1(a,b1,b2,c,x,y,algoflag,userflag,debug,val,hyp2f1)
  implicit none 
  logical ispossible
  logical debug
  integer(4) algoflag,isneg,pflag,userflag,flag,errorflag
  !integer get_transformation
  integer get_transformation2  
  integer check_parms
  integer hyp2f1
  real(8) x,y
  real(8) zero
  parameter(zero = 1d-18)
  complex(8) a,b1,b2,c,cgamma,f21,s,g2,cx,cy,cgammar,val
  complex(8) f1conv,sa,sb,sc,f1aux,g2aux,f1bnl  
  complex(8) g2gam

  cx=x*(1.0,0.0)
  cy=y*(1.0,0.0)

  algoflag=-1
  ispossible=.false.
  !
  !   Check particular cases
  !
  if (dabs(x-1d0).lt.zero .and. dabs(y-1d0).lt.zero) then 
     val = cgamma(c)*cgamma(c-a-b1-b2)*(cgammar(c-a)*cgammar(c-b1-b2))
     ispossible=.true.
  else if(dabs(x-y).lt.zero) then
     val = f21(a,b1+b2,c,cx,hyp2f1)
     ispossible=.true.
  else if(cdabs(a-c).lt.zero) then
     val = (1-x)**(-b1)*(1-y)**(-b2)
     ispossible=.true.
  else if(dabs(y-1d0).lt.zero) then
     val = cgamma(c)*cgamma(c-a-b2)*f21(a,b1,c-b2,cx,hyp2f1)*(cgammar(c-a)*cgammar(c-b2))
     ispossible=.true.
  else if(dabs(x-1d0).lt.zero) then
     val = cgamma(c)*cgamma(c-a-b1)*f21(a,b2,c-b1,cy,hyp2f1)*(cgammar(c-a)*cgammar(c-b1))
     ispossible=.true.
  else if(x.eq.0d0.or.dabs(x).lt.zero) then
     val = f21(a,b2,c,cy,hyp2f1)
     ispossible=.true.
  else if(y.eq.0d0.or.dabs(y).lt.zero) then
     val = f21(a,b1,c,cx,hyp2f1)
     ispossible=.true.
  else if (cdabs(b1).lt.zero) then
     val = f21(a,b2,c,cy,hyp2f1)
     ispossible=.true.
  else if (cdabs(b2).lt.zero) then
     val = f21(a,b1,c,cx,hyp2f1)
     ispossible=.true.
  else if (isneg(b1).eq.-1.or.isneg(b2).eq.-1.or.isneg(a).eq.-1) then
     ispossible=.true.
     val = f1conv(a,b1,b2,c,x,y,hyp2f1)
  endif

  if(ispossible) then
     algoflag = 5

     if(debug) call printline("Simple transformations possible: ")
     if(debug) call writef1(a,b1,b2,c,cx,cy,val)

     ! If the user did not supply a flag, return
     if(userflag.lt.0) then
        return
     else
        if(debug) call printline("But not returned due to user-specific flag")
     end if
  end if

  if(.not.ispossible) then
     !
     !   Check polynomial cases
     !
     if(isneg(a).le.0.or.isneg(b1).le.0.or.isneg(b2).le.0.or.isneg(c-a).le.0) then
        ispossible = .true.
        val = f1bnl(a,b1,b2,c,x,y,hyp2f1)
     end if
     !if(isneg(c-b1-b2).le.0) then
     !    ispossible = .true.    
     !    f1 = f1bnl(a,c-b1-b2,b2,c,x/(-1 + x),(x-y)/(-1 + x),hyp2f1)/(1-x)**a            
     !end if

     if(ispossible) then
        algoflag = 6

        if(debug) call printline("Polynomial transformations possible: ")
        if(debug) call writef1(a,b1,b2,c,cx,cy,val)

        ! If the user did not supply a flag, return
        if(userflag.lt.0) then
           return
        else
           if(debug) call printline("But not returned due to user-specific flag")
        end if
     end if
  end if

  if(.not.ispossible) then
     !
     ! Check regions
     !  
     algoflag = get_transformation2(x,y)
  end if

  ! If the user supplied a flag, than use this one
  if(userflag.ge.0) then
     flag = userflag
  else
     ! otherwise use the region algorithm flag
     flag = algoflag
  end if
  !
  ! Check parameters
  !
  pflag = check_parms(b1,b2,c,flag)

  if(flag.eq.0) then            
     if (debug) then 
        call printoutformat("Convergence r.:  x= ", x," y= ",y)
     end if
     s= f1conv(a,b1,b2,c,x,y,hyp2f1)
  else if(flag.eq.1) then  
     if (debug) then 
        call printoutformat("ODE Integr.:     x= ",x," y= ",y)
     end if
     call hypgeof1(a,b1,b2,c,cx,cy,s,errorflag,hyp2f1)
  else if(flag.eq.2) then
     if (debug) then
        call printoutformat("B&Ch Series:     x= ",x," y= ",y)
     end if
     s= f1bnl(a,b1,b2,c,x,y,hyp2f1)
  else if(flag.eq.15) then
     !
     !   Eq. (15)      
     !
     if (debug) call printoutformat("                 x= ",x," y= ",y)
     if (debug) call printoutformat("Transf. (15):    u= ",x/(x-1)," w= ",y/(y-1))
     s= f1conv(-a+c,b1,b2,c,x/(-1+x),y/(-1+y),hyp2f1)/((1-x)**b1*(1-y)**b2)
  else if(flag.eq.16) then
     !
     !   Eq. (16)      
     !
     if (debug) call printoutformat("                 x= ",x," y= ",y)
     if (debug) call printoutformat("Transf. (16):    u= ",x/(x-1)," w= ",(x-y)/(x-1))
     s= f1conv(a,-b1-b2 + c,b2,c,x/(-1 + x),(x-y)/(-1 + x),hyp2f1)/(1-x)**a
  else if(flag.eq.17) then
     !
     !   Eq. (17)      
     !
     if (debug) call printoutformat("                 x= ",x," y= ",y)
     if (debug) call printoutformat("Transf. (17):    u= ",y/(y-1)," w= ",(y-x)/(y-1))
     s= f1conv(a,b1,-b1-b2 + c,c,(-x + y)/(-1 + y),y/(-1 + y),hyp2f1)/(1-y)**a
  else if(flag.eq.21) then
     !
     !   Eq. (21)       Transformation in (1,1)x 
     !
     if (debug) call printoutformat("                 x= ",x," y= ",y)    
     if (debug) call printoutformat("Transf. (21):    u= ",1-x," w= ",1-y)    
     s=cgamma(c)*cgamma(c-a-b1-b2)*cgammar(c-a)*cgammar(c-b1-b2)                           &
          *f1conv(a,b1,b2,1+a+b1+b2-c,1-x,1-y,hyp2f1)+                                               &
          cgamma(c)*cgamma(a+b2-c)*(cgammar(a)*cgammar(b2))                                   &
          *(1-x)**(-b1)*(1-y)**(c-a-b2)*f1conv(c-a,b1,c-b1-b2,c-a-b2+1,(1-y)/(1-x),1-y,hyp2f1)+     &  
          cgamma(c)*cgamma(c-a-b2)*cgamma(a+b1+b2-c)*(cgammar(a)*cgammar(b1)*cgammar(c-a))   &
          *(1-x)**(c-a-b1-b2)*g2(c-b1-b2,b2,a+b1+b2-c,c-a-b2,x-1,(1-y)/(x-1),hyp2f1)
  else if(flag.eq.22)  then
     !
     !   Eq. (22)      Transformation in (1,1)y
     !
     if (debug) call printoutformat("                 x= ",x," y= ",y)    
     if (debug) call printoutformat("Transf. (22):    u= ",1-x," w= ",1-y)    
     sa = cgamma(c)*cgamma(c-a-b1-b2)*cgammar(c-a)*cgammar(c-b1-b2)*f1conv(a,b1,b2,1+a+b1+b2-c, 1-x,1-y,hyp2f1)
     sb = cgamma(c)*cgamma(a+b1-c)*(1-y)**(-b2)*(1-x)**(c-a-b1)                    &
          *f1conv(c-a,c-b2-b1,b2,c-a-b1+1, 1-x,(1-x)/(1-y),hyp2f1)*cgammar(a)*cgammar(b1)
     sc = cgamma(c)*cgamma(c-a-b1)*cgamma(a+b1+b2-c)*(1-y)**(c-a-b1-b2)              &
          *g2(b1,c-b1-b2,c-a-b1,a+b1+b2-c,(1-x)/(y-1),y-1,hyp2f1)*cgammar(a)*cgammar(b2)*cgammar(c-a)
     s = sa+sb+sc
  else if(flag.eq.23) then 
     !
     !   Eq. (23)      Transformation in (0,Inf.)
     !
     if (debug) call printoutformat("                 x= ",x," y= ",y)    
     if (debug) call printoutformat("Transf. (23):    u= ",x/y," w= ",1/y)
     if(pflag.eq.0) then 
        ! Normal case
        s= cgamma(c)*cgamma(b2-a)*(cgammar(b2)*cgammar(c-a))        &
             *(-y)**(-a)*f1conv(a,b1,1+a-c,a-b2+1,x/y,1/y,hyp2f1)+           &
             (cgamma(c)*cgamma(a-b2)*(cgammar(a)*cgammar(c-b2)))      &
             *(-y)**(-b2)*g2(b1,b2,1+b2-c,a-b2,-x,-1/y,hyp2f1)
     else
        ! c-b2=neg. int.
        s= cgamma(c)*cgamma(b2-a)*(cgammar(b2)*cgammar(c-a))        &
             *(-y)**(-a)*f1conv(a,b1,1+a-c,a-b2+1,x/y,1/y,hyp2f1)+           &
             (cgamma(c)*cgamma(a-b2)*cgammar(a))      &
             *(-y)**(-b2)*g2gam(b1,b2,1+b2-c,a-b2,-x,-1/y,hyp2f1)
     end if
  else if(flag.eq.24) then
     !
     !   Eq. (24)      Transformation in (Inf.0)
     !
     if (debug) call printoutformat("                 x= ",x," y= ",y)
     if (debug) call printoutformat("Transf. (24):    u= ",y/x," w= ",1/x)
     if(pflag.eq.0) then
        ! Normal case
        s=cgamma(c)*cgamma(b1-a)*cgammar(b1)*cgammar(c-a)     &
             *(-x)**(-a)*f1conv(a,1+a-c,b2,a-b1+1,1/x,y/x,hyp2f1)+        &  
             cgamma(c)*cgamma(a-b1)*(cgammar(a)* cgammar(c-b1))    & 
             *(-x)**(-b1)*g2(b1,b2,a-b1,1+b1-c,-1/x,-y,hyp2f1)      
     else   
        ! c-b1=neg. int.
        s=cgamma(c)*cgamma(b1-a)*cgammar(b1)*cgammar(c-a)     &
             *(-x)**(-a)*f1conv(a,1+a-c,b2,a-b1+1,1/x,y/x,hyp2f1)+        &  
             cgamma(c)*cgamma(a-b1)*cgammar(a)                   & 
             *(-x)**(-b1)*g2gam(b1,b2,a-b1,1+b1-c,-1/x,-y,hyp2f1)      
     end if
  else if(flag.eq.25) then
     !
     !   Eq. (25)      Transformation in (1,Inf.)
     !
     if (debug) call printoutformat("                 x= ",x," y= ",y)
     if (debug) call printoutformat("Transf. (25):    u= ",(1-x)/(1-y)," w= ",1/(1-y))    
     s=cgamma(c)*cgamma(b2-a)*(cgammar(c-a)*cgammar(b2))                                  &
          *(1-y)**(-a)*f1conv(a,b1,c-b1-b2,1+a-b2,(1-x)/(1-y),1/(1-y),hyp2f1)+                      &
          cgamma(c)*cgamma(a+b1-c)*(cgammar(a)*cgammar(b1))                                  &
          *(1-x)**(c-a-b1)*(1-y)**(-b2)*f1conv(c-a,b2,c-b2-b1,c-a-b1+ 1,(1-x)/(1-y),1-x,hyp2f1)+    &
          cgamma(c)*cgamma(a-b2)*cgamma(c-a-b1)*(cgammar(a)*cgammar(c-b1-b2)*cgammar(c-a))  &
          *(1-y)**(-b2)*g2(b1,b2,c-a-b1,a-b2,x-1,1/(y-1),hyp2f1)
  else if(flag.eq.26) then
     !
     !   Eq. (26)      Transformation in (Inf.,1)
     !
     if (debug) call printoutformat("                 x= ",x," y= ",y)
     if (debug) call printoutformat("Transf. (26):    u= ",(1-y)/(1-x)," w= ",1/(1-x))        
     s=cgamma(c)*cgamma(b1-a)*(cgammar(c-a)*cgammar(b1))                                    &
          *(1-x)**(-a)*f1conv(a,c-b1-b2,b2,1+a-b1,1/(1-x),(1-y)/(1-x),hyp2f1)+                        &
          cgamma(c)*cgamma(a+b2-c)*(cgammar(a)*cgammar(b2))                                    &
          *(1-y)**(c-a-b2)*(1-x)**(-b1)*f1conv(c-a,b1,c-b1-b2,c-a-b2+1,(1-y)/(1-x),1-y,hyp2f1)+      &
          cgamma(c)*cgamma(a-b1)*cgamma(c-a-b2)*(cgammar(a)*cgammar(c-b1-b2)*cgammar(c-a))    &
          *(1-x)**(-b1)*g2(b1,b2,a-b1,c-a-b2,1/(x-1),y-1,hyp2f1)
  else if(flag.eq.27) then
     !
     !   Eq. (27)      Transformation in (Inf.Inf.)x
     !
     if (debug) call printoutformat("                 x= ",x," y= ",y)    
     if (debug) call printoutformat("Transf. (27):    u= ",x/y," w= ",1/x) 
     s=cgamma(c)*cgamma(b2-a)*(cgammar(c-a)*cgammar(b2))                                &
          *(-y)**(-a)*f1conv(a,b1,1+a-c,1+a-b2,x/y,1/y,hyp2f1)+                                  &
          cgamma(c)*cgamma(a-b1-b2)*(cgammar(a)*cgammar(c-b1-b2))                          &
          *(-x)**(-b1)*(-y)**(-b2)*f1conv(1+b1+b2-c,b1,b2,1+b1+b2-a,1/x,1/y,hyp2f1)+              &
          cgamma(c)*cgamma(a-b2)*cgamma(b1+b2-a)*(cgammar(a)*cgammar(b1)*cgammar(c-a))    &
          *(-x)**(b2-a)*(-y)**(-b2)*g2(1+a-c,b2,b1+b2-a,a-b2,-1/x,-x/y,hyp2f1)
  else if(flag.eq.28) then
     !
     !   Eq. (28)      Transformation in (Inf.Inf.)y
     !
     if (debug) call printoutformat("                 x= ",x," y= ",y)    
     if (debug) call printoutformat("Transf. (28):    u= ",1/x," w= ",1/y) 
     s= cgamma(c)*cgamma(b1-a)*(cgammar(c-a)*cgammar(b1))                            &
          *(-x)**(-a)*f1conv(a,1+a-c,b2,1+a-b1,1/x,y/x,hyp2f1)+                                &
          cgamma(c)*cgamma(a-b1-b2)*(cgammar(a)*cgammar(c-b1-b2))                        &
          *(-x)**(-b1)*(-y)**(-b2)*f1conv(1+b1+b2-c,b1,b2,1+b1+b2-a,1/x,1/y,hyp2f1)+            &
          cgamma(c)*cgamma(a-b1)*cgamma(b1+b2-a)*(cgammar(a)*cgammar(b2)*cgammar(c-a))  &
          *(-x)**(-b1)*(-y)**(b1-a)*g2(b1,1+a-c,a-b1,b1+b2-a,-y/x,-1/y,hyp2f1)

  else if(flag.eq.29) then
     !
     !   Eq. (29)      Transformation in (Inf.Inf.)xx
     !
     if (debug) call printoutformat("                 x= ",x," y= ",y)    
     if (debug) call printoutformat("Transf. (29):    u= ",(x-y)/(y*(x-1))," w= ",1/y) 
     f1aux = f1conv(c-a,b1,1-a,1+b1+b2-a,(x-y)/(y*(x-1)),1/y,hyp2f1)
     g2aux = g2(b1,c-b1-b2,1-b1-b2,b1+b2-a,(x-y)/(1-x),-1/y,hyp2f1)
     s= cgamma(c)*cgamma(a-b1-b2)*(cgammar(a)*cgammar(c-b1-b2))                                      &
          *(-y)**(a-c)*(1-x)**(-b1)*(1-y)**(c-a-b2)*f1conv(c-a,b1,1-a,1+b1+b2-a,(x-y)/(y*(x-1)),1/y,hyp2f1)+    &
          cgamma(c)*cgamma(b1+b2-a)*(cgammar(c-a)*cgammar(b1+b2))                                        &
          *(-y)**(b1+b2-c)*(1-x)**(-b1)*(1-y)**(c-a-b2)*g2(b1,c-b1-b2,1-b1-b2,b1+b2-a,(x-y)/(1-x),-1/y,hyp2f1)  
  else if(flag.eq.30) then
     !
     !   Eq. (30)      Transformation in (Inf.Inf.)xx
     !
     if (debug) call printoutformat("                 x= ",x," y= ",y)
     if (debug) call printoutformat("Transf. (29):    u= ",(y-x)/(x*(y-1))," w= ",1/x)
     s=cgamma(c)*cgamma(a-b1-b2)*(cgammar(a)*cgammar(c-b1-b2))                                        &
          *(-x)**(a-c)*(1-y)**(-b2)*(1-x)**(c-a-b1)*f1conv(c-a,1-a,b2,1+b1+b2-a,1/x,(y-x)/(x*(y-1)),hyp2f1)+    &
          cgamma(c)*cgamma(b1+b2-a)*(cgammar(c-a)*cgammar(b1+b2))                                        &
          *(-x)**(b1+b2-c)*(1-y)**(-b2)*(1-x)**(c-a-b1)*g2(c-b1-b2,b2,b1+b2-a,1-b1-b2,-1/x,(y-x)/(1-y),hyp2f1)
  end if
  val = s

  if(flag.ge.0) then
     if(debug) call writef1(a,b1,b2,c,cx,cy,val)
  else
     call printline("This is not Possible")
     call writef1(a,b1,b2,c,cx,cy,val)
  end if

  return
end subroutine f1


