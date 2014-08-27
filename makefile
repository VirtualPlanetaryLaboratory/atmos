# Make file for coupled radiative/convective climate and photochemical models
# needs to be updated anytime either of the models makefile is changed

###############################
# Set variables
###############################

PLATFORM = osx

# COMPILATION FLAGS:
#FC = ifort
#FCFLAG = -i8 -r8 -align all -zero
# For debugging, use this flag set instead:
#FCFLAG = -i8 -r8 -debug -fp-stack-check -check bounds -g -traceback -heap-arrays -gen-interfaces -warn interfaces -check arg_temp_created -align all -Bstatic
# more debugging options
#FCFLAG = -r8 -i8 -diag-enable sv-include

#optimized
FC = gfortran
#FCFLAG = -O3 -I. -fno-automatic -ftree-vectorize -fdefault-integer-8 -fdefault-real-8
#debugging versions
#FCFLAG = -g -fbounds-check -Wall -fbacktrace -finit-real=nan
#FCFLAG = -g -fno-automatic -ftrapv -fbounds-check -O
#FCFLAG = -g -fno-automatic -Wuninitialized -fbounds-check -O -ftrapv
FCFLAG = -g -fno-automatic -Wuninitialized -fbounds-check -O -ftrapv -Wall -I. -fdefault-integer-8 -fdefault-real-8

#optimized
#FCFLAG = -O3 -fno-automatic -ftree-vectorize
# The above causes errors with intrinsic functions expecting kind=4 numbers. Instead, use:
#FCFLAG = -O3 -I. -fno-automatic -ftree-vectorize
# for debugging, use these flags (recommended by fortran90.org):
#FCFLAG = -Wall -Wextra -Wimplicit-interface -fPIC -fmax-errors=1 -g -fbacktrace -fbounds-check -fcheck-array-temporaries -fdefault-integer-8 -fdefault-real-8
# for production runs, use these flags (recommended by fortran90.org):
#FCFLAG = -Wall -Wextra -Wimplicit-interface -fPIC -Werror -fmax-errors=1 -O3 -march=native -ffast-math -funroll-loops -fdefault-integer-8 -fdefault-real-8

# try new compiler and flags with HAO linux machines
# FC = pgf90
# FCFLAG = -DLINUX -DMSS -r8

# excecutable file
EXECUTABLE = couple.exe

# Main files
PHOTMAIN = Photochem
COUPLED = couple
CSUR = Clima

# Main subdirectories
CDIRCL = CLIMA/
CDIRPH = PHOTOCHEM/
CSUB = SUBROUTINES/

# Subroutines for coupling
CPATH0 = COUPLE
COUTP = output_photo
CINPI = input_interp

# Files +_o
PHOTOBJ  = OBJECTFILES
CLIMOBJ = OBJECT_CLIMA
COUPOBJ = OBJECT_COUPLE

# Subdirectories (CAPS) and subroutines (lower case)

# Photochemical code files
# folder locations and matrix solvers
SUBPATH = PHOTOCHEM/SUBROUTINES
LINPACK = PHOTOCHEM/OBJECT_PHOT
LINPATH = PHOTOCHEM/LINPACK
BLOK = PHOTOCHEM/DATA/INCLUDE
INC = PHOTOCHEM/INPUTFILES

# Subroutines
PHOTGRID = Photgrid
RATES = Rates
RAIN = Rainout
AQUO = Aqueous
PHOT = Photo
INIT = Initphoto
DENS = Densty
DIFCO = Difco
SED = Sedmnt
SAT = PhotSatrat
AERT = Aertab
OUT = Output
AERC= Aercon
DOCHEM = Dochem
CHEMPL = Chempl
LTNG = Ltning
MSCAT = Mscat
XSEC = Xsections
YSN = Youngsun
MIE = Initmie
FRACMIE = Initmiefrac
LNUM = Lnum
RESET = ireset
TWOSTR= Twostr
SAXPY = saxpy
SGBFA = sgbfa
SGEFA = sgefa
SGBSL = sgbsl
SGESL = sgesl
SSCAL = sscal
ISMAX = isamax
SGBCO = sgbco
SGTSL = sgtsl
SDOT = sdot
SASUM = sasum
DAXPY = daxpy
DGBFA = dgbfa
DGBCO = dgbco
DGBSL = dgbsl
DSCAL = dscal
DSUM = dasum
DDOT = ddot
IDMAX = idamax
RAYL = Rayleigh
SPLI = Spline

# These are the climate code files
CPATH1 = SETUP
CPROF = profile
CREAD = readsol
CSTAR = pickstar
CGRID = grid
CIRES = irexpsums
CINTP = interp
CINOZ = interpozone
COZON = ozone

CPATH2 = CONVEC
CCONV = convec
CSATH = satrat
CRELH = relhum
CSATC = satco2

CPATH3 = RADTRANS
CGASC = gascon
CSOL  = solar
CINFR = ir
CRAYL = rayley
CDSTI = delta2strir
CDSTS = delta2str
CMTRX = tridag
CPLAN = planck

CPATH4 = PRTCL
CGRDA = gridaer
CADAT = aerabsdata
CINPA = interpar1

CPATHR = RRTM
CRTM = rrtm
CREG = rtreg
CRTR = rtr
CATM = rrtatm
CSET = setcoef
CTAU = taumol
CRGC = rtregcld
CRTC = rtrcld
CUTL = util_$(PLATFORM)
CEXT = extra
CRTX = rtrcldmr
CRGX = rtregcldmr
CKGS = k_g
CCLD = cldprop

INTCO2 = interpco2cia
INTH2H2 = interph2h2cia
INTH2N2 = interph2n2cia
INTIR = interpir
INTO2 = interpo2cia
INTSOL = interpsolar
IRM = irm
SOLARM = solarm
SOLMOX = solarmox
SOLOX = solarox

OBPATH = $(CDIRPH)$(PHOTMAIN).o \
	$(CDIRPH)$(PHOTOBJ)/$(PHOTGRID).o \
	$(CDIRPH)$(PHOTOBJ)/$(RATES).o \
	$(CDIRPH)$(PHOTOBJ)/$(RAIN).o \
	$(CDIRPH)$(PHOTOBJ)/$(AQUO).o \
	$(CDIRPH)$(PHOTOBJ)/$(PHOT).o \
	$(CDIRPH)$(PHOTOBJ)/$(DENS).o \
	$(CDIRPH)$(PHOTOBJ)/$(DIFCO).o \
	$(CDIRPH)$(PHOTOBJ)/$(SED).o \
	$(CDIRPH)$(PHOTOBJ)/$(SAT).o \
	$(CDIRPH)$(PHOTOBJ)/$(AERT).o \
	$(CDIRPH)$(PHOTOBJ)/$(OUT).o \
	$(CDIRPH)$(PHOTOBJ)/$(AERC).o \
	$(CDIRPH)$(PHOTOBJ)/$(DOCHEM).o \
	$(CDIRPH)$(PHOTOBJ)/$(CHEMPL).o \
	$(CDIRPH)$(PHOTOBJ)/$(LTNG).o \
	$(CDIRPH)$(PHOTOBJ)/$(MSCAT).o \
	$(CDIRPH)$(PHOTOBJ)/$(INIT).o \
	$(CDIRPH)$(PHOTOBJ)/$(XSEC).o \
	$(CDIRPH)$(PHOTOBJ)/$(YSN).o \
	$(CDIRPH)$(PHOTOBJ)/$(MIE).o \
	$(CDIRPH)$(PHOTOBJ)/$(FRACMIE).o \
	$(CDIRPH)$(PHOTOBJ)/$(LNUM).o \
        $(CDIRPH)$(PHOTOBJ)/$(RESET).o \
	$(CDIRPH)$(PHOTOBJ)/$(TWOSTR).o \
        $(CDIRPH)$(PHOTOBJ)/$(RAYL).o \
	$(CDIRPH)$(PHOTOBJ)/$(SPLI).o \
	$(LINPACK)/$(SAXPY).o \
	$(LINPACK)/$(DAXPY).o \
	$(LINPACK)/$(SGBFA).o \
	$(LINPACK)/$(DGBFA).o \
	$(LINPACK)/$(SGEFA).o \
	$(LINPACK)/$(SGBCO).o \
	$(LINPACK)/$(SDOT).o \
	$(LINPACK)/$(DGBCO).o \
	$(LINPACK)/$(SGBSL).o \
	$(LINPACK)/$(DGBSL).o \
	$(LINPACK)/$(SGESL).o \
	$(LINPACK)/$(SSCAL).o \
	$(LINPACK)/$(DSCAL).o \
	$(LINPACK)/$(SASUM).o \
	$(LINPACK)/$(DSUM).o \
	$(LINPACK)/$(DDOT).o \
	$(LINPACK)/$(ISMAX).o \
	$(LINPACK)/$(IDMAX).o \
	$(LINPACK)/$(SGTSL).o \
	$(COUPOBJ)/$(COUPLED).o \
	$(COUPOBJ)/$(COUTP).o \
	$(COUPOBJ)/$(CINPI).o \
	$(CDIRCL)$(CLIMOBJ)/$(CSUR).o \
	$(CDIRCL)$(CLIMOBJ)/$(CPROF).o \
	$(CDIRCL)$(CLIMOBJ)/$(CREAD).o \
	$(CDIRCL)$(CLIMOBJ)/$(CSTAR).o \
	$(CDIRCL)$(CLIMOBJ)/$(CGRID).o \
	$(CDIRCL)$(CLIMOBJ)/$(CIRES).o \
	$(CDIRCL)$(CLIMOBJ)/$(CINTP).o \
	$(CDIRCL)$(CLIMOBJ)/$(CINOZ).o \
	$(CDIRCL)$(CLIMOBJ)/$(COZON).o \
	$(CDIRCL)$(CLIMOBJ)/$(CCONV).o \
	$(CDIRCL)$(CLIMOBJ)/$(CSATH).o \
	$(CDIRCL)$(CLIMOBJ)/$(CRELH).o \
	$(CDIRCL)$(CLIMOBJ)/$(CSATC).o \
	$(CDIRCL)$(CLIMOBJ)/$(CGASC).o \
	$(CDIRCL)$(CLIMOBJ)/$(CSOL).o  \
	$(CDIRCL)$(CLIMOBJ)/$(CINFR).o \
	$(CDIRCL)$(CLIMOBJ)/$(CRAYL).o \
	$(CDIRCL)$(CLIMOBJ)/$(CDSTI).o \
	$(CDIRCL)$(CLIMOBJ)/$(CDSTS).o \
	$(CDIRCL)$(CLIMOBJ)/$(CMTRX).o \
	$(CDIRCL)$(CLIMOBJ)/$(CPLAN).o \
	$(CDIRCL)$(CLIMOBJ)/$(CGRDA).o \
	$(CDIRCL)$(CLIMOBJ)/$(CADAT).o \
	$(CDIRCL)$(CLIMOBJ)/$(CINPA).o \
	$(CDIRCL)$(CLIMOBJ)/$(CRTM).o  \
	$(CDIRCL)$(CLIMOBJ)/$(CREG).o \
	$(CDIRCL)$(CLIMOBJ)/$(CRTR).o \
	$(CDIRCL)$(CLIMOBJ)/$(CATM).o \
	$(CDIRCL)$(CLIMOBJ)/$(CSET).o \
	$(CDIRCL)$(CLIMOBJ)/$(CTAU).o \
	$(CDIRCL)$(CLIMOBJ)/$(CRGC).o \
	$(CDIRCL)$(CLIMOBJ)/$(CRTC).o \
	$(CDIRCL)$(CLIMOBJ)/$(CUTL).o \
	$(CDIRCL)$(CLIMOBJ)/$(CEXT).o \
	$(CDIRCL)$(CLIMOBJ)/$(CRTX).o \
	$(CDIRCL)$(CLIMOBJ)/$(CRGX).o \
	$(CDIRCL)$(CLIMOBJ)/$(CKGS).o \
	$(CDIRCL)$(CLIMOBJ)/$(CCLD).o \
	$(CDIRCL)$(CLIMOBJ)/$(INTCO2).o \
	$(CDIRCL)$(CLIMOBJ)/$(INTH2H2).o \
	$(CDIRCL)$(CLIMOBJ)/$(INTH2N2).o \
	$(CDIRCL)$(CLIMOBJ)/$(INTIR).o \
	$(CDIRCL)$(CLIMOBJ)/$(INTO2).o \
	$(CDIRCL)$(CLIMOBJ)/$(INTSOL).o \
	$(CDIRCL)$(CLIMOBJ)/$(IRM).o \
	$(CDIRCL)$(CLIMOBJ)/$(SOLARM).o \
	$(CDIRCL)$(CLIMOBJ)/$(SOLMOX).o \
	$(CDIRCL)$(CLIMOBJ)/$(SOLOX).o 

INCALL = $(BLOK)/PHOTABLOK.inc \
	$(BLOK)/BBLOK.inc \
	$(BLOK)/CBLOK.inc \
	$(BLOK)/DBLOK.inc \
	$(BLOK)/EBLOK.inc \
	$(BLOK)/FBLOK.inc \
	$(BLOK)/GBLOK.inc \
	$(BLOK)/JBLOK.inc \
	$(BLOK)/NBLOK.inc \
	$(BLOK)/RBLOK.inc \
	$(BLOK)/SBLOK.inc \
	$(BLOK)/ZBLOK.inc \
	$(BLOK)/LTBLOK.inc \
	$(BLOK)/AERBLK.inc \
	$(BLOK)/SULBLK.inc \
	$(BLOK)/PBLOK.inc \
	$(BLOK)/QBLOK.inc \
	$(BLOK)/RRATS.inc \
	$(BLOK)/SATBLK.inc \
	$(BLOK)/WBLOK.inc \
	$(INC)/parameters.inc

###############################
# Load line
###############################

$(EXECUTABLE) :	 $(OBPATH) $(INCALL)
	$(FC) $(FCFLAG) -o $(EXECUTABLE) $(OBPATH)

###############################

# Object compile lines
###############################
$(COUPOBJ)/$(COUPLED).o : $(COUPLED).f
	$(FC) $(FCFLAG) -c $(COUPLED).f
	\mv $(COUPLED).o $(COUPOBJ)

### Compiling programs for coupling

$(COUPOBJ)/$(COUTP).o : $(CPATH0)/$(COUTP).f
	$(FC) $(FCFLAG) -c $(CPATH0)/$(COUTP).f
	\mv $(COUTP).o $(COUPOBJ)

$(COUPOBJ)/$(CINPI).o : $(CPATH0)/$(CINPI).f
	$(FC) $(FCFLAG) -c $(CPATH0)/$(CINPI).f
	\mv $(CINPI).o $(COUPOBJ)

### Photochemical model ###

$(PHOTOBJ)/$(MAIN).o :      $(MAIN).f $(INCALL)
	$(FC) $(FCFLAG) -c $(MAIN).f
	mv $(MAIN).o $(PHOTOBJ)

$(CDIRPH)$(PHOTOBJ)/$(PHOTGRID).o :      $(SUBPATH)/$(PHOTGRID).f $(INC)/parameters.inc \
			$(BLOK)/PHOTABLOK.inc
	$(FC) $(FCFLAG) -c $(SUBPATH)/$(PHOTGRID).f
	mv $(PHOTGRID).o $(CDIRPH)$(PHOTOBJ)

$(CDIRPH)$(PHOTOBJ)/$(RATES).o :     $(SUBPATH)/$(RATES).f $(INC)/parameters.inc \
			$(BLOK)/PHOTABLOK.inc $(BLOK)/RBLOK.inc
	$(FC) $(FCFLAG) -c $(SUBPATH)/$(RATES).f
	mv $(RATES).o $(CDIRPH)$(PHOTOBJ)

$(CDIRPH)$(PHOTOBJ)/$(RAIN).o :      $(SUBPATH)/$(RAIN).f $(INC)/parameters.inc \
			$(BLOK)/PHOTABLOK.inc $(BLOK)/BBLOK.inc \
			$(BLOK)/DBLOK.inc $(BLOK)/GBLOK.inc \
			$(BLOK)/NBLOK.inc $(BLOK)/WBLOK.inc
	$(FC) $(FCFLAG) -c $(SUBPATH)/$(RAIN).f
	mv $(RAIN).o $(CDIRPH)$(PHOTOBJ)

$(CDIRPH)$(PHOTOBJ)/$(AQUO).o :      $(SUBPATH)/$(AQUO).f $(INC)/parameters.inc \
			$(BLOK)/WBLOK.inc
	$(FC) $(FCFLAG) -c $(SUBPATH)/$(AQUO).f
	mv $(AQUO).o $(CDIRPH)$(PHOTOBJ)

$(CDIRPH)$(PHOTOBJ)/$(PHOT).o :      $(SUBPATH)/$(PHOT).f $(INC)/parameters.inc \
			$(BLOK)/PHOTABLOK.inc $(BLOK)/CBLOK.inc \
			$(BLOK)/JBLOK.inc $(BLOK)/QBLOK.inc \
			$(BLOK)/RBLOK.inc $(BLOK)/SBLOK.inc \
			$(BLOK)/DBLOK.inc
	$(FC) $(FCFLAG) -c $(SUBPATH)/$(PHOT).f
	mv $(PHOT).o $(CDIRPH)$(PHOTOBJ)

$(CDIRPH)$(PHOTOBJ)/$(DENS).o :      $(SUBPATH)/$(DENS).f $(INC)/parameters.inc \
			$(BLOK)/PHOTABLOK.inc
	$(FC) $(FCFLAG) -c $(SUBPATH)/$(DENS).f
	mv $(DENS).o $(CDIRPH)$(PHOTOBJ)

$(CDIRPH)$(PHOTOBJ)/$(DIFCO).o :     $(SUBPATH)/$(DIFCO).f $(INC)/parameters.inc \
			$(BLOK)/PHOTABLOK.inc $(BLOK)/AERBLK.inc \
			$(BLOK)/CBLOK.inc
	$(FC) $(FCFLAG) -c $(SUBPATH)/$(DIFCO).f
	mv $(DIFCO).o $(CDIRPH)$(PHOTOBJ)

$(CDIRPH)$(PHOTOBJ)/$(SED).o :       $(SUBPATH)/$(SED).f $(INC)/parameters.inc \
			$(BLOK)/PHOTABLOK.inc $(BLOK)/GBLOK.inc \
			$(BLOK)/NBLOK.inc $(BLOK)/AERBLK.inc \
			$(BLOK)/CBLOK.inc
	$(FC) $(FCFLAG) -c $(SUBPATH)/$(SED).f
	mv $(SED).o $(CDIRPH)$(PHOTOBJ)

$(CDIRPH)$(PHOTOBJ)/$(SAT).o :       $(SUBPATH)/$(SAT).f $(INC)/parameters.inc \
			$(BLOK)/PHOTABLOK.inc $(BLOK)/SATBLK.inc
	$(FC) $(FCFLAG) -c $(SUBPATH)/$(SAT).f
	mv $(SAT).o $(CDIRPH)$(PHOTOBJ)

$(CDIRPH)$(PHOTOBJ)/$(AERT).o :      $(SUBPATH)/$(AERT).f $(INC)/parameters.inc \
			$(BLOK)/PHOTABLOK.inc $(BLOK)/SULBLK.inc \
			$(BLOK)/AERBLK.inc
	$(FC) $(FCFLAG) -c $(SUBPATH)/$(AERT).f
	mv $(AERT).o $(CDIRPH)$(PHOTOBJ)

$(CDIRPH)$(PHOTOBJ)/$(OUT).o :       $(SUBPATH)/$(OUT).f $(INC)/parameters.inc \
			$(BLOK)/PHOTABLOK.inc $(BLOK)/BBLOK.inc \
			$(BLOK)/CBLOK.inc $(BLOK)/DBLOK.inc \
			$(BLOK)/FBLOK.inc $(BLOK)/GBLOK.inc \
			$(BLOK)/JBLOK.inc $(BLOK)/NBLOK.inc \
			$(BLOK)/RBLOK.inc $(BLOK)/SBLOK.inc \
			$(BLOK)/WBLOK.inc $(BLOK)/ZBLOK.inc \
			$(BLOK)/SATBLK.inc $(BLOK)/SULBLK.inc \
			$(BLOK)/AERBLK.inc $(BLOK)/RRATS.inc
	$(FC) $(FCFLAG) -c $(SUBPATH)/$(OUT).f
	mv $(OUT).o $(CDIRPH)$(PHOTOBJ)

$(CDIRPH)$(PHOTOBJ)/$(AERC).o :      $(SUBPATH)/$(AERC).f $(INC)/parameters.inc \
			$(BLOK)/PHOTABLOK.inc $(BLOK)/SATBLK.inc \
			$(BLOK)/SULBLK.inc
	$(FC) $(FCFLAG) -c $(SUBPATH)/$(AERC).f
	mv $(AERC).o $(CDIRPH)$(PHOTOBJ)

$(CDIRPH)$(PHOTOBJ)/$(DOCHEM).o :    $(SUBPATH)/$(DOCHEM).f $(INC)/parameters.inc \
			$(BLOK)/PHOTABLOK.inc $(BLOK)/BBLOK.inc \
			$(BLOK)/CBLOK.inc $(BLOK)/DBLOK.inc \
			$(BLOK)/GBLOK.inc $(BLOK)/NBLOK.inc \
			$(BLOK)/RBLOK.inc $(BLOK)/ZBLOK.inc \
			$(BLOK)/LTBLOK.inc $(BLOK)/SATBLK.inc \
			$(BLOK)/SULBLK.inc $(BLOK)/RRATS.inc
	$(FC) $(FCFLAG) -c $(SUBPATH)/$(DOCHEM).f
	mv $(DOCHEM).o $(CDIRPH)$(PHOTOBJ)

$(CDIRPH)$(PHOTOBJ)/$(CHEMPL).o :    $(SUBPATH)/$(CHEMPL).f $(INC)/parameters.inc \
			$(BLOK)/RBLOK.inc
	$(FC) $(FCFLAG) -c $(SUBPATH)/$(CHEMPL).f
	mv $(CHEMPL).o $(CDIRPH)$(PHOTOBJ)

$(CDIRPH)$(PHOTOBJ)/$(LTNG).o :      $(SUBPATH)/$(LTNG).f $(INC)/parameters.inc \
			$(BLOK)/BBLOK.inc $(BLOK)/LTBLOK.inc \
			$(BLOK)/NBLOK.inc $(BLOK)/PHOTABLOK.inc
	$(FC) $(FCFLAG) -c $(SUBPATH)/$(LTNG).f
	mv $(LTNG).o $(CDIRPH)$(PHOTOBJ)

$(CDIRPH)$(PHOTOBJ)/$(MSCAT).o :     $(SUBPATH)/$(MSCAT).f $(INC)/parameters.inc \
			$(BLOK)/EBLOK.inc $(BLOK)/JBLOK.inc  \
			$(BLOK)/PHOTABLOK.inc $(BLOK)/AERBLK.inc \
			$(BLOK)/CBLOK.inc  $(BLOK)/NBLOK.inc
	$(FC) $(FCFLAG) -c $(SUBPATH)/$(MSCAT).f
	mv $(MSCAT).o $(CDIRPH)$(PHOTOBJ)

$(CDIRPH)$(PHOTOBJ)/$(INIT).o :      $(SUBPATH)/$(INIT).f $(INC)/parameters.inc \
			$(BLOK)/PHOTABLOK.inc $(BLOK)/QBLOK.inc \
			$(BLOK)/JBLOK.inc $(BLOK)/DBLOK.inc
	$(FC) $(FCFLAG) -c $(SUBPATH)/$(INIT).f
	mv $(INIT).o $(CDIRPH)$(PHOTOBJ)

$(CDIRPH)$(PHOTOBJ)/$(XSEC).o :      $(SUBPATH)/$(XSEC).f $(INC)/parameters.inc \
			$(BLOK)/QBLOK.inc $(BLOK)/PBLOK.inc
	$(FC) $(FCFLAG) -c $(SUBPATH)/$(XSEC).f
	mv $(XSEC).o $(CDIRPH)$(PHOTOBJ)

$(CDIRPH)$(PHOTOBJ)/$(YSN).o :       $(SUBPATH)/$(YSN).f $(INC)/parameters.inc \
			$(BLOK)/QBLOK.inc
	$(FC) $(FCFLAG) -c $(SUBPATH)/$(YSN).f
	mv $(YSN).o $(CDIRPH)$(PHOTOBJ)

$(CDIRPH)$(PHOTOBJ)/$(MIE).o :       $(SUBPATH)/$(MIE).f $(INC)/parameters.inc
	$(FC) $(FCFLAG) -c $(SUBPATH)/$(MIE).f
	mv $(MIE).o $(CDIRPH)$(PHOTOBJ)

$(CDIRPH)$(PHOTOBJ)/$(FRACMIE).o :       $(SUBPATH)/$(FRACMIE).f $(INC)/parameters.inc
	$(FC) $(FCFLAG) -c $(SUBPATH)/$(FRACMIE).f
	mv $(FRACMIE).o $(CDIRPH)$(PHOTOBJ)

$(CDIRPH)$(PHOTOBJ)/$(LNUM).o :      $(SUBPATH)/$(LNUM).f $(INC)/parameters.inc \
			$(BLOK)/NBLOK.inc
	$(FC) $(FCFLAG) -c $(SUBPATH)/$(LNUM).f
	mv $(LNUM).o $(CDIRPH)$(PHOTOBJ)

$(CDIRPH)$(PHOTOBJ)/$(RESET).o :    $(SUBPATH)/$(RESET).f
	$(FC) $(FCFLAG) -c $(SUBPATH)/$(RESET).f
	\mv $(RESET).o $(CDIRPH)$(PHOTOBJ)

$(CDIRPH)$(PHOTOBJ)/$(TWOSTR).o :    $(SUBPATH)/$(TWOSTR).f $(INC)/parameters.inc \
			$(BLOK)/PHOTABLOK.inc $(BLOK)/CBLOK.inc \
			$(BLOK)/AERBLK.inc
	$(FC) $(FCFLAG) -c $(SUBPATH)/$(TWOSTR).f
	mv $(TWOSTR).o $(CDIRPH)$(PHOTOBJ)

$(CDIRPH)$(PHOTOBJ)/$(RAYL).o : $(SUBPATH)/$(RAYL).f
	$(FC) $(FCFLAG) -c $(SUBPATH)/$(RAYL).f
	\mv $(RAYL).o $(CDIRPH)$(PHOTOBJ)

$(CDIRPH)$(PHOTOBJ)/$(SPLI).o : $(SUBPATH)/$(SPLI).f
	$(FC) $(FCFLAG) -c $(SUBPATH)/$(SPLI).f
	\mv $(SPLI).o $(CDIRPH)$(PHOTOBJ)

$(LINPACK)/$(SAXPY).o : $(LINPATH)/$(SAXPY).f
	$(FC) $(FCFLAG) -c $(LINPATH)/$(SAXPY).f
	mv $(SAXPY).o $(LINPACK)

$(LINPACK)/$(DAXPY).o : $(LINPATH)/$(DAXPY).f
	$(FC) $(FCFLAG) -c $(LINPATH)/$(DAXPY).f
	mv $(DAXPY).o $(LINPACK)

$(LINPACK)/$(SGBFA).o : $(LINPATH)/$(SGBFA).f
	$(FC) $(FCFLAG) -c $(LINPATH)/$(SGBFA).f
	mv $(SGBFA).o $(LINPACK)

$(LINPACK)/$(DGBFA).o : $(LINPATH)/$(DGBFA).f
	$(FC) $(FCFLAG) -c $(LINPATH)/$(DGBFA).f
	mv $(DGBFA).o $(LINPACK)

$(LINPACK)/$(SGEFA).o : $(LINPATH)/$(SGEFA).f
	$(FC) $(FCFLAG) -c $(LINPATH)/$(SGEFA).f
	mv $(SGEFA).o $(LINPACK)

$(LINPACK)/$(SGBCO).o : $(LINPATH)/$(SGBCO).f
	$(FC) $(FCFLAG) -c $(LINPATH)/$(SGBCO).f
	mv $(SGBCO).o $(LINPACK)

$(LINPACK)/$(DGBCO).o : $(LINPATH)/$(DGBCO).f
	$(FC) $(FCFLAG) -c $(LINPATH)/$(DGBCO).f
	mv $(DGBCO).o $(LINPACK)

$(LINPACK)/$(SGBSL).o : $(LINPATH)/$(SGBSL).f
	$(FC) $(FCFLAG) -c $(LINPATH)/$(SGBSL).f
	mv $(SGBSL).o $(LINPACK)

$(LINPACK)/$(DGBSL).o : $(LINPATH)/$(DGBSL).f
	$(FC) $(FCFLAG) -c $(LINPATH)/$(DGBSL).f
	mv $(DGBSL).o $(LINPACK)

$(LINPACK)/$(SGESL).o : $(LINPATH)/$(SGESL).f
	$(FC) $(FCFLAG) -c $(LINPATH)/$(SGESL).f
	mv $(SGESL).o $(LINPACK)

$(LINPACK)/$(SSCAL).o : $(LINPATH)/$(SSCAL).f
	$(FC) $(FCFLAG) -c $(LINPATH)/$(SSCAL).f
	mv $(SSCAL).o $(LINPACK)

$(LINPACK)/$(DSCAL).o : $(LINPATH)/$(DSCAL).f
	$(FC) $(FCFLAG) -c $(LINPATH)/$(DSCAL).f
	mv $(DSCAL).o $(LINPACK)

$(LINPACK)/$(SDOT).o :  $(LINPATH)/$(SDOT).f
	$(FC) $(FCFLAG) -c $(LINPATH)/$(SDOT).f
	mv $(SDOT).o $(LINPACK)

$(LINPACK)/$(SASUM).o : $(LINPATH)/$(SASUM).f
	$(FC) $(FCFLAG) -c $(LINPATH)/$(SASUM).f
	mv $(SASUM).o $(LINPACK)

$(LINPACK)/$(DSUM).o :  $(LINPATH)/$(DSUM).f
	$(FC) $(FCFLAG) -c $(LINPATH)/$(DSUM).f
	mv $(DSUM).o $(LINPACK)

$(LINPACK)/$(DDOT).o :  $(LINPATH)/$(DDOT).f
	$(FC) $(FCFLAG) -c $(LINPATH)/$(DDOT).f
	mv $(DDOT).o $(LINPACK)

$(LINPACK)/$(ISMAX).o : $(LINPATH)/$(ISMAX).f
	$(FC) $(FCFLAG) -c $(LINPATH)/$(ISMAX).f
	mv $(ISMAX).o $(LINPACK)

$(LINPACK)/$(IDMAX).o : $(LINPATH)/$(IDMAX).f
	$(FC) $(FCFLAG) -c $(LINPATH)/$(IDMAX).f
	mv $(IDMAX).o $(LINPACK)

$(LINPACK)/$(SGTSL).o : $(LINPATH)/$(SGTSL).f
	$(FC) $(FCFLAG) -c $(LINPATH)/$(SGTSL).f
	mv $(SGTSL).o $(LINPACK)

$(CDIRPH)$(PHOTMAIN).o : $(CDIRPH)$(PHOTMAIN).f
	$(FC) $(FCFLAG) -c $(CDIRPH)$(PHOTMAIN).f
	\mv $(PHOTMAIN).o $(CDIRPH)


### Clima model ###
$(CDIRCL)$(CLIMOBJ)/$(CSUR).o : $(CDIRCL)$(CSUR).f
	$(FC) $(FCFLAG) -c $(CDIRCL)$(CSUR).f
	\mv $(CSUR).o $(CDIRCL)$(CLIMOBJ)

$(CDIRCL)$(CLIMOBJ)/$(CPROF).o : $(CDIRCL)$(CPATH1)/$(CPROF).f
	$(FC) $(FCFLAG) -c $(CDIRCL)$(CPATH1)/$(CPROF).f
	\mv $(CPROF).o $(CDIRCL)$(CLIMOBJ)

$(CDIRCL)$(CLIMOBJ)/$(CREAD).o : $(CDIRCL)$(CPATH1)/$(CREAD).f
	$(FC) $(FCFLAG) -c $(CDIRCL)$(CPATH1)/$(CREAD).f
	\mv $(CREAD).o $(CDIRCL)$(CLIMOBJ)

$(CDIRCL)$(CLIMOBJ)/$(CSTAR).o : $(CDIRCL)$(CPATH1)/$(CSTAR).f
	$(FC) $(FCFLAG) -c $(CDIRCL)$(CPATH1)/$(CSTAR).f
	\mv $(CSTAR).o $(CDIRCL)$(CLIMOBJ)

$(CDIRCL)$(CLIMOBJ)/$(CGRID).o : $(CDIRCL)$(CPATH1)/$(CGRID).f
	$(FC) $(FCFLAG) -c $(CDIRCL)$(CPATH1)/$(CGRID).f
	\mv $(CGRID).o $(CDIRCL)$(CLIMOBJ)

$(CDIRCL)$(CLIMOBJ)/$(CIRES).o : $(CDIRCL)$(CPATH1)/$(CIRES).f
	$(FC) $(FCFLAG) -c $(CDIRCL)$(CPATH1)/$(CIRES).f
	\mv $(CIRES).o $(CDIRCL)$(CLIMOBJ)

$(CDIRCL)$(CLIMOBJ)/$(CINTP).o : $(CDIRCL)$(CPATH1)/$(CINTP).f
	$(FC) $(FCFLAG) -c $(CDIRCL)$(CPATH1)/$(CINTP).f
	\mv $(CINTP).o $(CDIRCL)$(CLIMOBJ)

$(CDIRCL)$(CLIMOBJ)/$(CINOZ).o : $(CDIRCL)$(CPATH1)/$(CINOZ).f
	$(FC) $(FCFLAG) -c $(CDIRCL)$(CPATH1)/$(CINOZ).f
	\mv $(CINOZ).o $(CDIRCL)$(CLIMOBJ)

$(CDIRCL)$(CLIMOBJ)/$(COZON).o : $(CDIRCL)$(CPATH1)/$(COZON).f
	$(FC) $(FCFLAG) -c $(CDIRCL)$(CPATH1)/$(COZON).f
	\mv $(COZON).o $(CDIRCL)$(CLIMOBJ)

$(CDIRCL)$(CLIMOBJ)/$(CCONV).o : $(CDIRCL)$(CPATH2)/$(CCONV).f
	$(FC) $(FCFLAG) -c $(CDIRCL)$(CPATH2)/$(CCONV).f
	\mv $(CCONV).o $(CDIRCL)$(CLIMOBJ)

$(CDIRCL)$(CLIMOBJ)/$(CSATH).o : $(CDIRCL)$(CPATH2)/$(CSATH).f
	$(FC) $(FCFLAG) -c $(CDIRCL)$(CPATH2)/$(CSATH).f
	\mv $(CSATH).o $(CDIRCL)$(CLIMOBJ)

$(CDIRCL)$(CLIMOBJ)/$(CRELH).o : $(CDIRCL)$(CPATH2)/$(CRELH).f
	$(FC) $(FCFLAG) -c $(CDIRCL)$(CPATH2)/$(CRELH).f
	\mv $(CRELH).o $(CDIRCL)$(CLIMOBJ)

$(CDIRCL)$(CLIMOBJ)/$(CSATC).o : $(CDIRCL)$(CPATH2)/$(CSATC).f
	$(FC) $(FCFLAG) -c $(CDIRCL)$(CPATH2)/$(CSATC).f
	\mv $(CSATC).o $(CDIRCL)$(CLIMOBJ)

$(CDIRCL)$(CLIMOBJ)/$(CGASC).o : $(CDIRCL)$(CPATH3)/$(CGASC).f
	$(FC) $(FCFLAG) -c $(CDIRCL)$(CPATH3)/$(CGASC).f
	\mv $(CGASC).o $(CDIRCL)$(CLIMOBJ)

$(CDIRCL)$(CLIMOBJ)/$(CSOL).o : $(CDIRCL)$(CPATH3)/$(CSOL).f
	$(FC) $(FCFLAG) -c $(CDIRCL)$(CPATH3)/$(CSOL).f
	\mv $(CSOL).o $(CDIRCL)$(CLIMOBJ)

$(CDIRCL)$(CLIMOBJ)/$(CINFR).o : $(CDIRCL)$(CPATH3)/$(CINFR).f
	$(FC) $(FCFLAG) -c $(CDIRCL)$(CPATH3)/$(CINFR).f
	\mv $(CINFR).o $(CDIRCL)$(CLIMOBJ)

$(CDIRCL)$(CLIMOBJ)/$(CRAYL).o : $(CDIRCL)$(CPATH3)/$(CRAYL).f
	$(FC) $(FCFLAG) -c $(CDIRCL)$(CPATH3)/$(CRAYL).f
	\mv $(CRAYL).o $(CDIRCL)$(CLIMOBJ)

$(CDIRCL)$(CLIMOBJ)/$(CDSTI).o : $(CDIRCL)$(CPATH3)/$(CDSTI).f
	$(FC) $(FCFLAG) -c $(CDIRCL)$(CPATH3)/$(CDSTI).f
	\mv $(CDSTI).o $(CDIRCL)$(CLIMOBJ)

$(CDIRCL)$(CLIMOBJ)/$(CDSTS).o : $(CDIRCL)$(CPATH3)/$(CDSTS).f
	$(FC) $(FCFLAG) -c $(CDIRCL)$(CPATH3)/$(CDSTS).f
	\mv $(CDSTS).o $(CDIRCL)$(CLIMOBJ)

$(CDIRCL)$(CLIMOBJ)/$(CMTRX).o : $(CDIRCL)$(CPATH3)/$(CMTRX).f
	$(FC) $(FCFLAG) -c $(CDIRCL)$(CPATH3)/$(CMTRX).f
	\mv $(CMTRX).o $(CDIRCL)$(CLIMOBJ)

$(CDIRCL)$(CLIMOBJ)/$(CPLAN).o : $(CDIRCL)$(CPATH3)/$(CPLAN).f
	$(FC) $(FCFLAG) -c $(CDIRCL)$(CPATH3)/$(CPLAN).f
	\mv $(CPLAN).o $(CDIRCL)$(CLIMOBJ)

$(CDIRCL)$(CLIMOBJ)/$(CGRDA).o : $(CDIRCL)$(CPATH4)/$(CGRDA).f
	$(FC) $(FCFLAG) -c $(CDIRCL)$(CPATH4)/$(CGRDA).f
	\mv $(CGRDA).o $(CDIRCL)$(CLIMOBJ)

$(CDIRCL)$(CLIMOBJ)/$(CADAT).o : $(CDIRCL)$(CPATH4)/$(CADAT).f
	$(FC) $(FCFLAG) -c $(CDIRCL)$(CPATH4)/$(CADAT).f
	\mv $(CADAT).o $(CDIRCL)$(CLIMOBJ)

$(CDIRCL)$(CLIMOBJ)/$(CINPA).o : $(CDIRCL)$(CPATH4)/$(CINPA).f
	$(FC) $(FCFLAG) -c $(CDIRCL)$(CPATH4)/$(CINPA).f
	\mv $(CINPA).o $(CDIRCL)$(CLIMOBJ)

# RRTM files
$(CDIRCL)$(CLIMOBJ)/$(CRTM).o : $(CDIRCL)$(CPATHR)/$(CRTM).f 
	$(FC) $(FCFLAG) -c $(CDIRCL)$(CPATHR)/$(CRTM).f
	\mv $(CRTM).o $(CDIRCL)$(CLIMOBJ)

$(CDIRCL)$(CLIMOBJ)/$(CRTR).o : $(CDIRCL)$(CPATHR)/$(CRTR).f
	$(FC) $(FCFLAG) -c $(CDIRCL)$(CPATHR)/$(CRTR).f
	\mv $(CRTR).o $(CDIRCL)$(CLIMOBJ)

$(CDIRCL)$(CLIMOBJ)/$(CREG).o : $(CDIRCL)$(CPATHR)/$(CREG).f
	$(FC) $(FCFLAG) -c $(CDIRCL)$(CPATHR)/$(CREG).f
	\mv $(CREG).o $(CDIRCL)$(CLIMOBJ)

$(CDIRCL)$(CLIMOBJ)/$(CATM).o : $(CDIRCL)$(CPATHR)/$(CATM).f
	$(FC) $(FCFLAG) -c $(CDIRCL)$(CPATHR)/$(CATM).f
	\mv $(CATM).o $(CDIRCL)$(CLIMOBJ)

$(CDIRCL)$(CLIMOBJ)/$(CSET).o : $(CDIRCL)$(CPATHR)/$(CSET).f
	$(FC) $(FCFLAG) -c $(CDIRCL)$(CPATHR)/$(CSET).f
	\mv $(CSET).o $(CDIRCL)$(CLIMOBJ)

$(CDIRCL)$(CLIMOBJ)/$(CTAU).o : $(CDIRCL)$(CPATHR)/$(CTAU).f
	$(FC) $(FCFLAG) -c $(CDIRCL)$(CPATHR)/$(CTAU).f
	\mv $(CTAU).o $(CDIRCL)$(CLIMOBJ)

$(CDIRCL)$(CLIMOBJ)/$(CRGC).o : $(CDIRCL)$(CPATHR)/$(CRGC).f
	$(FC) $(FCFLAG) -c $(CDIRCL)$(CPATHR)/$(CRGC).f
	\mv $(CRGC).o $(CDIRCL)$(CLIMOBJ)

$(CDIRCL)$(CLIMOBJ)/$(CRTC).o : $(CDIRCL)$(CPATHR)/$(CRTC).f
	$(FC) $(FCFLAG) -c $(CDIRCL)$(CPATHR)/$(CRTC).f
	\mv $(CRTC).o $(CDIRCL)$(CLIMOBJ)

$(CDIRCL)$(CLIMOBJ)/$(CUTL).o : $(CDIRCL)$(CPATHR)/$(CUTL).f
	$(FC) $(FCFLAG) -c $(CDIRCL)$(CPATHR)/$(CUTL).f
	\mv $(CUTL).o $(CDIRCL)$(CLIMOBJ)

$(CDIRCL)$(CLIMOBJ)/$(CEXT).o : $(CDIRCL)$(CPATHR)/$(CEXT).f
	$(FC) $(FCFLAG) -c $(CDIRCL)$(CPATHR)/$(CEXT).f
	\mv $(CEXT).o $(CDIRCL)$(CLIMOBJ)

$(CDIRCL)$(CLIMOBJ)/$(CRTX).o : $(CDIRCL)$(CPATHR)/$(CRTX).f
	$(FC) $(FCFLAG) -c $(CDIRCL)$(CPATHR)/$(CRTX).f
	\mv $(CRTX).o $(CDIRCL)$(CLIMOBJ)

$(CDIRCL)$(CLIMOBJ)/$(CRGX).o : $(CDIRCL)$(CPATHR)/$(CRGX).f
	$(FC) $(FCFLAG) -c $(CDIRCL)$(CPATHR)/$(CRGX).f
	\mv $(CRGX).o $(CDIRCL)$(CLIMOBJ)

$(CDIRCL)$(CLIMOBJ)/$(CCLD).o : $(CDIRCL)$(CPATHR)/$(CCLD).f
	$(FC) $(FCFLAG) -c $(CDIRCL)$(CPATHR)/$(CCLD).f
	\mv $(CCLD).o $(CDIRCL)$(CLIMOBJ)

$(CDIRCL)$(CLIMOBJ)/$(CKGS).o : $(CDIRCL)$(CPATHR)/$(CKGS).f
	$(FC) $(FCFLAG) -c $(CDIRCL)$(CPATHR)/$(CKGS).f
	\mv $(CKGS).o $(CDIRCL)$(CLIMOBJ)

$(CDIRCL)$(CLIMOBJ)/$(INTCO2).o : $(CDIRCL)$(CPATH1)/$(INTCO2).f
	$(FC) $(FCFLAG) -c $(CDIRCL)$(CPATH1)/$(INTCO2).f
	\mv $(INTCO2).o $(CDIRCL)$(CLIMOBJ)

$(CDIRCL)$(CLIMOBJ)/$(INTH2H2).o : $(CDIRCL)$(CPATH1)/$(INTH2H2).f
	$(FC) $(FCFLAG) -c $(CDIRCL)$(CPATH1)/$(INTH2H2).f
	\mv $(INTH2H2).o $(CDIRCL)$(CLIMOBJ)

$(CDIRCL)$(CLIMOBJ)/$(INTH2N2).o : $(CDIRCL)$(CPATH1)/$(INTH2N2).f
	$(FC) $(FCFLAG) -c $(CDIRCL)$(CPATH1)/$(INTH2N2).f
	\mv $(INTH2N2).o $(CDIRCL)$(CLIMOBJ)

$(CDIRCL)$(CLIMOBJ)/$(INTIR).o : $(CDIRCL)$(CPATH1)/$(INTIR).f
	$(FC) $(FCFLAG) -c $(CDIRCL)$(CPATH1)/$(INTIR).f
	\mv $(INTIR).o $(CDIRCL)$(CLIMOBJ)

$(CDIRCL)$(CLIMOBJ)/$(INTO2).o : $(CDIRCL)$(CPATH1)/$(INTO2).f
	$(FC) $(FCFLAG) -c $(CDIRCL)$(CPATH1)/$(INTO2).f
	\mv $(INTO2).o $(CDIRCL)$(CLIMOBJ)

$(CDIRCL)$(CLIMOBJ)/$(INTSOL).o : $(CDIRCL)$(CPATH1)/$(INTSOL).f
	$(FC) $(FCFLAG) -c $(CDIRCL)$(CPATH1)/$(INTSOL).f
	\mv $(INTSOL).o $(CDIRCL)$(CLIMOBJ)

$(CDIRCL)$(CLIMOBJ)/$(IRM).o : $(CDIRCL)$(CPATH3)/$(IRM).f
	$(FC) $(FCFLAG) -c $(CDIRCL)$(CPATH3)/$(IRM).f
	\mv $(IRM).o $(CDIRCL)$(CLIMOBJ)

$(CDIRCL)$(CLIMOBJ)/$(SOLARM).o : $(CDIRCL)$(CPATH3)/$(SOLARM).f
	$(FC) $(FCFLAG) -c $(CDIRCL)$(CPATH3)/$(SOLARM).f
	\mv $(SOLARM).o $(CDIRCL)$(CLIMOBJ)

$(CDIRCL)$(CLIMOBJ)/$(SOLMOX).o : $(CDIRCL)$(CPATH3)/$(SOLMOX).f
	$(FC) $(FCFLAG) -c $(CDIRCL)$(CPATH3)/$(SOLMOX).f
	\mv $(SOLMOX).o $(CDIRCL)$(CLIMOBJ)

$(CDIRCL)$(CLIMOBJ)/$(SOLOX).o : $(CDIRCL)$(CPATHR3/$(SOLOX).f
	$(FC) $(FCFLAG) -c $(CDIRCL)$(CPATH3)/$(SOLOX).f
	\mv $(SOLOX).o $(CDIRCL)$(CLIMOBJ)

###############################
# Cleanup/Re-make Options
###############################

clean : 
	rm $(CDIRPH)$(PHOTMAIN).o
	rm $(CDIRPH)$(PHOTOBJ)/*.o
	rm $(CDIRCL)$(CLIMOBJ)/*.o
	rm $(COUPOBJ)/*.o
	rm $(EXECUTABLE)

new:
	$(MAKE) clean
	$(MAKE)

conclude :
	echo
	echo '================='
	echo '  Makefile done'
	echo '================='

