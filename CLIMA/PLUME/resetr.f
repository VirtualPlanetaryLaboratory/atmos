      subroutine resetr (a, n, val)
   
c        Sets real array a(n) to val
        
      dimension a(*)
            
      do 10 i=1,n 
        a(i) = val
   10 continue
  
      return
      end   

