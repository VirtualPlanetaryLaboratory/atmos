      subroutine scopy (n, arr, inca, brr, incb)
      dimension arr(n), brr(n)

      j = 1
      do 10 i=1,n,inca
        brr(j) = arr(i)
        j = j + incb
   10 continue
      return
      end

c************************************************

      function cvmgt (vala, valb, test)
      logical test
      
      if (test) then
        cvmgt = vala
      else
        cvmgt = valb
      endif
      return
      end
