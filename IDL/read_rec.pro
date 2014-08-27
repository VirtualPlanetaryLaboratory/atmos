pro read_rec, fn, ifm, idat, skip, ntskip, mincol, maxcol, nx, mxrec, x
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;	p u r p o s e :						      ;;
;;								      ;;
;;	this IDL procedure a single record from FORTRAN formatted     ;;
;;      or binary tabular data file 	 			      ;;
;;								      ;;
;;      i n p u t :                                                   ;;
;;								      ;;
;;      fn - file name (string)                                       ;;
;;      ifm - input file format index (1) ascii (2) binary            ;;
;;      idat - input data type index  (1) string, (2) i*2  (3) i*4    ;;
;;             (4) real*4 (5) real*8                                  ;;
;;      skip - number of records to skip at top of file               ;;
;;      ntskip - index type of variable to skip over (see idat)       ;;
;;      mincol - first column to read                                 ;;
;;      maxcol - last column to read                                  ;;
;;      mxrec - maximum number of records to read from file (real)    ;;
;;								      ;;
;;      o u t p u t :                                                 ;;
;;								      ;;
;;      nx - number of values read from each record                   ;;
;;   mxrec - number of recores actually read from file                ;;
;;     x - vector with data from record                               ;; 
;;								      ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; 
ns = strarr(1)
;
START:
;
  nx = maxcol - mincol + 1
  nimax = mxrec
  nirec = mxrec
  x0 = fltarr(nx,nirec)
;
;     s k i p    u n n e c e s s a r y    r e c o r d s
;
nrec = 0
on_ioerror, FINISH
if(ifm eq 1) then begin
;
  openr, unit, fn, /get_lun;
;
  skipped = 0
  SKIP_RECF:
    skipped = skipped + 1.0
    if(skipped le skip) then readf, unit, ns
  if(skipped lt skip) then goto, SKIP_RECF
endif
;
if(ifm eq 2) then begin
;
  openr, unit, fn, /f77_unformatted, /get_lun
;
  skipped = 0.0
  SKIP_RECU:
    skipped = skipped + 1.0
    if(skipped le skip) then readu, unit, ns
  if(skipped lt skip) then goto, SKIP_RECU
endif
;
;    s e t    d i m e n s i o n    o f    a r r a y s
;
if (mincol gt 1) then begin
  if (ifm eq 1) then dummy = bytarr(mincol - 1)
  if (ifm eq 2) then begin
    if(ntskip eq 2) then dummy = intarr(mincol-1)
    if(ntskip eq 3) then dummy = lonarr(mincol-1)
    if(ntskip eq 4) then dummy = fltarr(mincol-1)
    if(ntskip eq 5) then dummy = dblarr(mincol-1)
  endif
endif
; 
if(idat eq 1) then begin
  arr = bytarr(nx)
endif
if(idat eq 2) then begin
  arr = intarr(nx)
endif
if(idat eq 3) then begin
  arr = lonarr(nx)
endif
if(idat eq 4) then begin
  arr = fltarr(nx)
endif
if(idat eq 5) then begin
  arr = dblarr(nx)
endif
;
;    r e a d    x     a n d    y    d a t a 
;
on_ioerror, FINISH
;
cmin = 0
cmax = nx - 1
nx0 = 0
ny0 = 0
if(ifm eq 1) then begin
  if(mincol gt 1) then begin
    nrec = 0.0
    READ_FMT_1:
      if (not EOF(unit)) then begin
        readf, unit, dummy,arr
;
;		calculate the indices of output value
;
        i = nrec
        x0(0:nx-1,i) = arr(cmin:cmax)
;	print, i,cmin,cmax,dummy,x0(0,i),x0(nx-1,i)
        nrec = nrec + 1.0
      if(nrec lt mxrec) then goto, READ_FMT_1
    endif
  endif
  if(mincol le 1) then begin
    nrec = 0.0
    READ_FMT_2:
      if (not EOF(unit)) then begin
        readf, unit, arr
;
;		calculate the indices of output value
;
        i = nrec
        x0(0:nx-1,i) = arr(cmin:cmax)
        nrec = nrec + 1.0
;	print, i,cmin,cmax,nrec,x0(0,i),x0(nx-1,i)
        if(nrec lt mxrec) then goto, READ_FMT_2
      endif
  endif
endif
;
if(ifm eq 2) then begin
  if(mincol gt 1) then begin
    nrec = 0.0
    READ_UNF_1:
      if (not EOF(unit)) then begin
        readu, unit, dummy,arr
;
;		calculate the indices of output value
;
        i = nrec
        x0(0:nx-1,i) = arr(cmin:cmax)
;	print, i,cmin,cmax,dummy,x0(0,i),x0(nx-1,i)
        nrec = nrec + 1.0
      if(nrec lt mxrec) then goto, READ_UNF_1
    endif
  endif

  if(mincol le 1) then begin
    nrec = 0.0
    READ_UNF_2:
      if (not EOF(unit)) then begin
        readu, unit, arr
;
;		calculate the indices of output value
;
        i = nrec
        x0(0:nx-1,i) = arr(cmin:cmax)
;	print, i,cmin,cmax,x0(0,i),x0(nx-1,i)
        nrec = nrec + 1.0
        if(nrec lt mxrec) then goto, READ_UNF_2
      endif
  endif
endif
;
;    c l o s e    u n i t    a n d    e x i t
;
FINISH:
mxrec = nrec
if(n_elements(unit) ne 0) then begin
  close, unit
  free_lun, unit
endif
;
i = nrec - 1
if (i lt 0) then i = 0
print, ' '
print, 'number of records   = ',mxrec
if(mxrec-1 ge 0) then begin
  print, 'x(0,0)          =',x0(0,0)
  if(nx gt 0) then print, 'x(nx-1,0)       =',x0(nx-1,0)
  if(nx gt 0 and i gt 0) then begin
    print, 'x(0,i,mxrec-1)        =',x0(0,i)
    print, 'x(nx-1,i,mxrec-1)     =',x0(nx-1,i)
    print, ' '
  endif
;
x = x0(0:nx-1,0:nirec-1)
   
endif 
;
return
end
