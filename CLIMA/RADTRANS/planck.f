      FUNCTION PLANCK(VA,T,HP,C,HK)
        ARG = HK*VA/T
        ARG = AMIN1(ARG,80.)
        EB = EXP(ARG)
        PLANCK = 2*HP*VA*VA*VA/C/C/(EB - 1.)
        RETURN
      END

      FUNCTION DPLANK(VA,T,BCON,HK)
        EB = EXP(HK*VA/T)
        VTB = VA*VA/(T*(EB - 1.))
        DPLANK = BCON*HK*EB*VTB*VTB
        RETURN
      END
