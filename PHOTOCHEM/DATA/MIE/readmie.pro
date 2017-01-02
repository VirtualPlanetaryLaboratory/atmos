psname='test.ps'
    th=5 & xyth = 8 & ct=4 & cs=1.5 
    ctl=4 & csl=1.25
    set_plot, 'ps'
    !p.font=0
    device, encapsul=0, /helv,/isolatin1, filename=psname
    device, /color

openw, 1, 'soltotalSUL.dat'
openw, 2, 'irtotalSUL.dat'

readcol, 'kastingIRgrid.dat', IRwn,IRwl, /silent   ;microns
readcol, 'kastingGridSolar.dat', SOLwl, /silent    ;microns

rads=['0.01','0.02','0.03','0.04','0.05','0.06','0.07','0.08','0.09','0.10','0.20','0.30','0.40','0.50','0.60','0.70','0.80','0.90','1.00']

nrad=n_elements(rads)
nsol=n_elements(SOLwl)
nir=n_elements(Irwl)

print, nsol,nrad, nir

loadct, 12
for i=0,nrad-1 do begin
;i=0
fn='h2so4_mono_'+rads[i]+'um.mie'

readcol, fn, wl,rind,cind,sig_ext,sig_ca, sig_abs, qext,qsca,qabs,pizer0,cosbar, /silent

;add points at the end so the interpolations don't go haywire
wl=[wl,500.]
;add a linear extrapolation of the 
nmie=n_elements(qext)
npts=2
X=wl[nmie-npts:nmie-1]
Y=alog10(qext[nmie-npts:nmie-1])
res=linfit(X,Y)


qext=[qext,10^(res[0]+res[1]*500.)]

Y=alog10(pizer0[nmie-npts:nmie-1])
res=linfit(X,Y)

pizer0=[pizer0,10^(res[0]+res[1]*500.)]

Y=alog10(cosbar[nmie-npts:nmie-1])
res=linfit(X,Y)
cosbar=[cosbar,10^(res[0]+res[1]*500.)]


;plot, wl, qext,  xthick=xyth,ythick=xyth,charthick=ct,charsize=cs,thick=th,title='qext: '+fn, /xlog, /ylog 
qextSOL=abs(interpol(qext, wl,SOLwl, /spline))
qextIR=abs(interpol(qext, wl,IRwl))
;oplot, SOLwl, qextSOL, psym=2, color=80,thick=th
;oplot, IRwl, qextIR, psym=1, color=200,thick=th

;print, pizer0[0:10]
;plot, wl, pizer0, xthick=xyth,ythick=xyth,charthick=ct,charsize=cs,thick=th,title='ssa: '+fn, /xlog,/ylog 
ssaSOL=abs(interpol(pizer0, wl,SOLwl, /spline))
ssaIR=abs(interpol(pizer0, wl,IRwl))
;oplot, SOLwl, ssaSOL, psym=2, color=80,thick=th
;oplot, IRwl, ssaIR, psym=1, color=200,thick=th


plot, wl, cosbar, xthick=xyth,ythick=xyth,charthick=ct,charsize=cs,thick=th,title='asym: '+fn, /xlog,/ylog 
asymSOL=abs(interpol(cosbar, wl,SOLwl, /spline))
asymIR=abs(interpol(cosbar, wl,IRwl))  ;linear interpolation works best here.  slightly overpredicts at the longest wavelengths, but does it in a better fashion than quadratic or spline

oplot, SOLwl, asymSOL, psym=2, color=80,thick=th
oplot, IRwl, asymIR, psym=1, color=200,thick=th


;create a data file for Jim's climate model, for the "solar" grid
printf, 1, '     SIZE INTERVAL   ' + strtrim(i+1,2)
printf, 1, '     WAVl            QEXT       W0        gf'
for j=0,nsol-1 do begin
printf, 1, SOLwl[j],qextSOL[j],ssaSOL[j],asymSOL[j]
endfor  

;create a data file for Jim's climate model, for the "IR" grid
printf, 2, '     SIZE INTERVAL   ' + strtrim(i+1,2)
printf, 2, '     WAVl            QEXT       W0        gf'
for j=0,nir-1 do begin
printf, 2, IRwl[j],qextIR[j],ssaIR[j],asymIR[j]
endfor  ;end loop over wavelength to create the file for Jim's climate model




endfor  ;end loop over each paricle size

close, 1
close, 2

device, /close
set_plot, 'x'
!p.font=0
!p.multi=0
littlesalmonspawn='gv '+psname+'&'
spawn, littlesalmonspawn


stop
end
