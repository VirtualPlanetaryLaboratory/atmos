      PARAMETER(NZ=200, NQ=39, NP=2, NQ1=NQ)   !particles in TRI_DIAG
      PARAMETER(NEQ=NQ*NZ,LDA=3*NQ+1)
      PARAMETER(NR=233, NSP=50, NSP1=NSP+1, NSP2=NSP+2, NMAX=65)
      PARAMETER(kw=250,kj=40)    !kw = maximum number of wl intervals, kj=number of photolysis reactions
      PARAMETER(ks=26) ! ks= # of photo. species
      PARAMETER(NAQ=10) !number of aqueous phase species
      PARAMETER(NT=50)   !number of temperatures in sulfate/H2O vapor pressure file (DATA/aerosol.table)
      PARAMETER(NF=50)   !NT=number of pressures per temperature in DATA/aerosol.table
      PARAMETER(ML=12,ML1=ML+1,ML2=2*ML)
      PARAMETER(kin=33)      !file unit for cross sections (opened and closed each usage)
