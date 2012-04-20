! subroutine printoutformat
! 20/04/2012
! modified by Daniel Sabanes Bove
! to include in R-package appell
subroutine writef1(ca,cb,cbp,cc,cx,cy,f1)
  implicit none
  complex(8) ca,cb,cc,cbp,cx,cy,f1

  ! debug = .false.
  ! if(debug) then 
  !     write(*,*) 
  !           write(*,*)	'a  =',ca
  !           write(*,*)	'b  =',cb
  !           write(*,*)	'bp =',cbp
  !           write(*,*)	'c  =',cc
  !           write(*,*)	'x  =',cx  
  !           write(*,*)	'y  =',cy  
  !           write(*,*)	'f1 =',f1  
  ! end if
  call writecomplex('a = ',ca)
  call writecomplex('b1= ',cb)
  call writecomplex('b2= ',cbp)
  call writecomplex('c = ',cc)
  call writecomplex('x = ',cx)
  call writecomplex('y = ',cy)
  call printline('---------------------')
  call writecomplex('f1= ',f1)          
  call printline('---------------------')
  return
end subroutine writef1



! subroutine printoutformat
! 20/04/2012
! modified by Daniel Sabanes Bove
! to include in R-package appell
subroutine writecomplex(str, ca)
  implicit none
  character(4) str
  complex(8) ca
  real(8) ar,ai

  call printstring(str)

  ar  = dreal(ca)
  ai  = dimag(ca)

  call printdouble(ar)

  if(dabs(ai).gt.0d0) then    
     if(ai.lt.0d0) then
        call printstring(" - ")
     else
        call printstring(" + ")
     end if

     ai = dabs(ai)
     call printdouble(ai)
     call printstring(" i")
  end if

  call printline("")
  return
end subroutine writecomplex















