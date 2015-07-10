      PARAMETER(NZ=160, NQ=13, NP=2,NQ1=NQ+NP)   !particles in TRI_DIAG
      PARAMETER(NEQ=NQ*NZ,LDA=3*NQ+1)
      PARAMETER(NR=197, NSP=65, NSP1=NSP+1, NSP2=NSP+2, NMAX=65)
      PARAMETER(kw=250,kj=45)    !kw = maximum number of wl intervals, kj=number of photolysis reactions
      PARAMETER(ks=32) ! ks= # of photo. species
      PARAMETER(NAQ=10) !number of aqueous phase species
      PARAMETER(NT=50)   !number of temperatures in sulfate/H2O vapor pressure file (DATA/aerosol.table)
      PARAMETER(NF=50)   !NT=number of pressures per temperature in DATA/aerosol.table
      PARAMETER(ML=12,ML1=ML+1,ML2=2*ML)
      PARAMETER(kin=33)      !file unit for cross sections (opened and closed each usage)

!used in full up model - testing for now - SO2 and S*O2 in species list, S*O2 and SO2 reactions, but no mixing