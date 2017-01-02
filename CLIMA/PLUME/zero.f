      subroutine zero (arr, nar)
c
c        zeros nar words starting at arr(1)
c     
      dimension arr(nar)
      do 10 j=1,nar
        arr(j) = 0.
   10 continue
      return
      end
