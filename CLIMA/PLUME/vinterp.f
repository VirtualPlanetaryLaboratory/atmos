      subroutine vinterp (a, sig, aprof, sigprof, nlon,norec,nlev,nprof)
      
c     Vertically interpolates basic state a(nlev) from aprof(nprof)
          
      dimension a(nlon,norec,nlev),sig(nlev),aprof(nprof),sigprof(nprof)
            
      do 100 jk=1,nlev
   
        if (sig(jk).le.sigprof(1)) then
          jkpa = 1
          jkpb = 1
          weia = 1.
        else if (sig(jk).ge.sigprof(nprof)) then
          jkpa = nprof   
          jkpb = nprof   
          weia = 1.
        else
          do 10 jkp=2,nprof
            if (sig(jk).le.sigprof(jkp)) then
              jkpa= jkp-1
              jkpb= jkp
              weia= (sigprof(jkp)-sig(jk))/(sigprof(jkp)-sigprof(jkp-1))
              goto 12
            endif
   10     continue
        endif
   12   continue

        do 20 jj=1,norec
          do 22 ji=1,nlon
            a(ji,jj,jk) = weia*aprof(jkpa) + (1.-weia)*aprof(jkpb)
   22     continue
   20   continue

  100 continue

      return
      end
