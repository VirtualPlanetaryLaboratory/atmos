      PARAMETER(kw=250,kj=62)    !kw = maximum number of wl intervals, kj=number of photolysis reactions
      PARAMETER(ks=47) ! ks= # of photo. species
      PARAMETER(NAQ=10) !number of aqueous phase species
      PARAMETER(NT=50)   !number of temperatures in sulfate/H2O vapor pressure file (DATA/aerosol.table)
      PARAMETER(NFP=50)   !NT=number of pressures per temperature in DATA/aerosol.table
      PARAMETER(ML=12,ML1=ML+1,ML2=2*ML)
      PARAMETER(kin=269)      !file unit for cross sections (opened and closed each usage)
