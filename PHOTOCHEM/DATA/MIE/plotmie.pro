;taken from readmie.pro, which was used to create datafiles for the
;photochemical model, as well as makes plots as a function of all
;particle wavelengths.  Here we are going hack this down to make plots
;just for 0.2 microns for presentation in the paper.


psname='test.ps'
    th=5 & xyth = 8 & ct=4 & cs=1.5 
    ctl=4 & csl=1.25
    set_plot, 'ps'
    !p.font=0
    device, encapsul=0, /helv,/isolatin1, filename=psname
    device, /color

readcol, 'kastingIRgrid.dat', IRwn,IRwl, /silent   ;microns
readcol, 'kastingGridSolar.dat', SOLwl, /silent    ;microns

loadct, 12


;i=0
fn='h2so4_mono_0.20um.mie'

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

xrange=[0.1,250]
xt='wavelength (microns)'



plot, wl, qext,  xthick=xyth,ythick=xyth,charthick=ct,charsize=cs,thick=th, /xlog, /ylog,yrange=[0.01,5],xrange=xrange, /xstyle, xtitle=xt,/ystyle 

qextSOL=abs(interpol(qext, wl,SOLwl, /spline))
qextIR=abs(interpol(qext, wl,IRwl))
oplot, SOLwl, qextSOL, psym=2, color=80,thick=th
oplot, IRwl, qextIR, psym=1, color=200,thick=th

;print, pizer0[0:10]
plot, wl, pizer0, xthick=xyth,ythick=xyth,charthick=ct,charsize=cs,thick=th, /xlog,/ylog,xtitle=xt,xrange=xrange,yrange=[1e-4,2],/xstyle, /ystyle 
ssaSOL=abs(interpol(pizer0, wl,SOLwl, /spline))
ssaIR=abs(interpol(pizer0, wl,IRwl))
oplot, SOLwl, ssaSOL, psym=2, color=80,thick=th
oplot, IRwl, ssaIR, psym=1, color=200,thick=th


plot, wl, cosbar, xthick=xyth,ythick=xyth,charthick=ct,charsize=cs,thick=th, /xlog,/ylog,xtitle=xt,xrange=xrange,yrange=[1e-4,2],/xstyle, /ystyle  
asymSOL=abs(interpol(cosbar, wl,SOLwl, /spline))
asymIR=abs(interpol(cosbar, wl,IRwl))  ;linear interpolation works best here.  slightly overpredicts at the longest wavelengths, but does it in a better fashion than quadratic or spline

oplot, SOLwl, asymSOL, psym=2, color=80,thick=th
oplot, IRwl, asymIR, psym=1, color=200,thick=th



device, /close
set_plot, 'x'
!p.font=0
!p.multi=0
littlesalmonspawn='gv '+psname+'&'
spawn, littlesalmonspawn


stop
end
