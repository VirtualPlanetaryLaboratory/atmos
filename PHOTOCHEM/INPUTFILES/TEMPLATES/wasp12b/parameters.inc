      PARAMETER(NZ=100, NQ=18, NP=0,NQ1=NQ)   !particles in main loop (for tri-diag NQ1=NQ+NP)
      PARAMETER(NEQ=NQ*NZ,LDA=3*NQ+1)
      PARAMETER(NR=179, NSP=21, NSP1=NSP+1, NSP2=NSP+2, NMAX=70)
      PARAMETER(kw=200,kj=10)    !kw = maximum number of wl intervals, kj=number of photolysis reactions
      PARAMETER(ks=5) ! ks= # of photo. species
      PARAMETER(NAQ=0) !number of aqueous phase species
      PARAMETER(NT=34)   !number of temperatures in sulfate/H2O vapor pressure file (DATA/aerosol.table)
      PARAMETER(NF=34)   !NT=number of pressures per temperature in DATA/aerosol.table
      PARAMETER(ML=12,ML1=ML+1,ML2=2*ML)
      PARAMETER(kin=33)      !file unit for cross sections (opened and closed each usage)
      