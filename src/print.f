C     subroutine printoutformat
C     20/04/2012
C     by Daniel Sabanes Bove
C     
C     This routine replaces:
C     
C     character*30 outformat
C     outformat = '(a20,f17.8,a4,f17.8)'
C     write(*,outformat) 'B&Ch Series:     x= ',x,' y= ',y
C     
C     with
C     
C     call printoutformat("B&Ch Series:     x= ",x," y= ",y)
C     
      subroutine printoutformat(label1, val1, label2, val2)
      implicit none
C     
      character(20) label1
      character(4) label2
      real(8) val1,val2
C     
      call printstring(label1)
      call printdouble(val1)
      call printstring(label2)
      call printdouble(val2)
      call printline("")
C     
      return
      end subroutine printoutformat
