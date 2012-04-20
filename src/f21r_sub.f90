!-------------------------------------------------------------------
!	subroutine: f21_sub
!
!       author: Daniel Sabanes Bove (daniel.sabanesbove@ifspm.uzh.ch)
!
!       date: 19/04/2012
!
!       purpose: This subroutine just calls the corresponding function
!		 f21 and is only necessary as an interface to R.
!
!       
!-------------------------------------------------------------------       



subroutine f21_sub(a,b,c,z,hyp2f1,val)
  implicit none

  complex(8) a,b,c,z,val
  complex(8) f21
  integer hyp2f1

  val = f21(a,b,c,z,hyp2f1)

  return
end subroutine f21_sub


