epsname='so2.ps'

    th=5 & xyth = 8 & ct=4 & cs=1.5 
    ctl=4 & csl=1.25
    set_plot, 'ps'
    !p.font=0
    device, encapsul=0, /helv,/isolatin1, filename=epsname
    device, /color

loadct, 13

readcol, 'SO2_zahnle.abs', wave,so2old,sdflkj,asdfasd, /silent
wave=wave/10.  ;convert A to nm


readcol, '2007jd009695-ds01.txt', wn,wl,s32,d32,s33,d33,s34,d34, /silent
;wl in nm, note also some negative vlues for cross sections (errors?)

;plot=1 is cross section comparison
;plot=2 is test of new grid
;plot=3 is Danielache isotopologue data
;plot=4 is test of new grid over isotopolgues 1A and 0.5A
;plot=5 is plot of isotope fractionation factors
;plot=6 shows errors
plot=5
int=1  ;for plot 5, look at integrated area under Farquhar's laser's
broadband=0 ;if broadband=1, use D2 lamp from Masterson 2011, else use laser Farquhar 2001


if plot eq 1 then begin
colors=[0,80,160,240]
psym=[0,0,0,0]
lines=[0,0,0,0]

xri=[190,210]
yri=[1e-18,2e-17]


;xri=[180,360]
;yri=[1e-22,2e-17]


plot, wave,so2old,xrange=xri,yrange=yri,/ylog,xtitle='Nanometers',charthick=ct,xthick=xyth,ythick=xyth,thick=th,ytitle='cross section (cm!U2!N)',/xstyle,/ystyle,lines=lines[0],psym=psym[0]

oplot, wave, so2old, thick=th+3, psym=4

oplot, wl, s32, thick=th, psym=psym[1], color=colors[1],lines=lines[1]

oplot, wl, s33, thick=th, psym=psym[2], color=colors[2],lines=lines[2]

oplot, wl, s34, thick=th, psym=psym[3], color=colors[3],lines=lines[3]

label=['old SO2 XS','!U32!NSO!D2!N','!U33!NSO!D!N','!U34!NSO!D!N']

   legend,label, thick=th+1,charthick=ctl,charsize=csl, /bottom,/left,lines=lines,psym=psym,color=colors,textcolor=colors,spacing=1.5

endif ;end plot 1


if plot eq 2 then begin
colors=[0,80,160,240]
psym=[0,4,1,5]
lines=[0,0,0,0]

xri=[190,210]
yri=[1e-18,2e-17]



plot, wl,s32,xrange=xri,yrange=yri,/ylog,xtitle='Nanometers',charthick=ct,xthick=xyth,ythick=xyth,thick=th,ytitle='cross section (cm!U2!N)',/xstyle,/ystyle,lines=lines[0],psym=psym[0]

oplot, wave, so2old, thick=th, psym=psym[1], color=colors[1],lines=lines[1]

readcol, 'SO2_HR2oldGrid',wave,HR2old, /silent
wave=wave/10.  ;convert to nm

oplot, wave, HR2old, thick=th, psym=psym[2], color=colors[2],lines=lines[2]

readcol, 'SO2_HRnewGrid',wave,HR2new, /silent
wave=wave/10.  ;convert to nm


oplot, wave, HR2new, thick=th, psym=psym[3], color=colors[3],lines=lines[3]

label=['high res !U32!NSO!D2!N xs','old xs','HighRes2OldGrid','HighRes2NewGrid']

   legend,label, thick=th+1,charthick=ctl,charsize=csl, /bottom,/left,lines=lines,psym=psym,color=colors,textcolor=colors,spacing=1.5

endif ;end plot 2


if plot eq 3 then begin
colors=[80,160,240]
psym=[4,1,5]
lines=[0,0,0]

xri=[200,205]
yri=[1e-18,2e-17]


;xri=[220,360]
;yri=[1e-22,2e-17]

title='Danielache data'


plot, wl,s32,xrange=xri,yrange=yri,/ylog,xtitle='Nanometers',charthick=ct,xthick=xyth,ythick=xyth,thick=th,ytitle='cross section (cm!U2!N)',/xstyle,/ystyle,lines=lines[0],psym=psym[0], /nodata,title=title

oplot, wl, s32, thick=th, color=colors[0],lines=0
oplot, wl, s32, thick=th, psym=psym[0], color=colors[0]

oplot, wl, s33, thick=th, color=colors[1],lines=0
oplot, wl, s33, thick=th, psym=psym[1], color=colors[1]


oplot, wl, s34, thick=th, color=colors[2],lines=0
oplot, wl, s34, thick=th, psym=psym[2], color=colors[2]



label=['!U32!NSO!D2!N','!U332!NSO!D2!N','!U34!NSO!D2!N']

   legend,label, thick=th+1,charthick=ctl,charsize=csl, /bottom,/left,lines=lines,psym=psym,color=colors,textcolor=colors,spacing=1.5

endif ;end plot 3


if plot eq 4 then begin
colors=[80,160,240]
psym=[1,1,1]
lines=[0,0,0]

xri=[200,205]
yri=[1e-18,2e-17]

;title='Danielache data to 1 A grid'
title='Danielache data 0.5 A grid'


plot, wl,s32,xrange=xri,yrange=yri,/ylog,xtitle='Nanometers',charthick=ct,xthick=xyth,ythick=xyth,thick=th,ytitle='cross section (cm!U2!N)',/xstyle,/ystyle,lines=lines[0],psym=psym[0], /nodata,title=title

oplot, wl, s32, thick=th, color=colors[0],lines=0
;readcol, 'SO2_HRnewGrid',wave,HR2new ;1A grid
readcol, 'HRnew32',wave,HR2new, /silent ;0.5 A grid
wave=wave/10.  ;convert to nm
oplot, wave, HR2new, thick=th, psym=psym[0], color=colors[0]


oplot, wl, s33, thick=th, color=colors[1],lines=0
readcol, 'SO2_HRnewGrid33',wave,HR2new, /silent ;1A grid
readcol, 'HRnew33',wave,HR2new, /silent ;0.5A grid
wave=wave/10.  ;convert to nm
oplot, wave, HR2new, thick=th, psym=psym[1], color=colors[1]




oplot, wl, s34, thick=th, color=colors[2],lines=0
;readcol, 'SO2_HRnewGrid34',wave,HR2new  ;1A
readcol, 'HRnew34',wave,HR2new, /silent   ;0.5A
wave=wave/10.  ;convert to nm
oplot, wave, HR2new, thick=th, psym=psym[2], color=colors[2]



label=['!U32!NSO!D2!N','!U332!NSO!D2!N','!U34!NSO!D2!N']

   legend,label, thick=th+1,charthick=ctl,charsize=csl, /bottom,/left,lines=lines,psym=psym,color=colors,textcolor=colors,spacing=1.5

endif ;end plot 4



if plot eq 5 then begin
colors=[80,240]
psym=[-1,-1]
lines=[0,0]

xri=[200,210]
;xri=[190,220]
xri=[190,215]


;title='Danielache data to 1 A grid'
title='0.5 Angstrom grid fractionation factors'
;title='longer wavelength fractionation factors (old grid)'


abscicon=1  ;plot of just the D33 factors
if abscicon eq 1 then title=''
;if abscicon eq 1 then title='isotopologue spectra shifted from !U32!NSO!D2!N'
if abscicon eq 1 then yri=[-500,500]
;if abscicon eq 1 then xri=[200,210]

if int then begin
xri=[192,194]
;xri=[185,200]
yri=[-700,700]
ystyle=9
endif

xri=[192,194]
;xri=[185,200]
yri=[-700,700]
ystyle=1

if broadband then xri=[184,220]




plot, wl,s32,xrange=xri,yrange=yri,xtitle='Nanometers',charthick=ct,charsize=cs,xthick=xyth,ythick=xyth,thick=th,ytitle='fractionation factor (permil)',/xstyle,ystyle=ystyle,lines=lines[0],psym=psym[0], /nodata,title=title,xmargin=[10,6]


readcol, 'HRnew32',wave,HR2new32, /silent ;0.5 A grid
readcol, 'HRnew33',wave,HR2new33, /silent ;0.5A grid
readcol, 'HRnew34',wave,HR2new34, /silent   ;0.5A
wave=wave/10.  ;convert to nm


if int then begin
loadct, 0
chop=where(wave ge 192.5 and wave le 193.5)

if broadband then chop=where(wave ge 184.0 and wave le 220.0)

n=n_elements(chop)
gauss=exp(-1*((findgen(n)-float(n)/2)/(float(n)/6))^2)
axis, yaxis=1,ythick=xyth,charthick=ct,charsize=cs,yrange=[0,1], ystyle=1, ytitle='laser transmission function',/save


if broadband then begin
;use Masterson broadband instead of Farquhar laser

;series 4 and 3 of Masterson Laser data (from spreadsheet in ~reinterpretingMIF
;series 4 - AA5-AA95 1800-1950 A
s4=[-3.761696e1,3.573691e-2,-8.218244e-6]

;series 3 - X5-X340  1940-2500 A
s3=[-1.989012e2,3.420474e-1,-2.206156e-4,6.365773e-8,-6.943965e-12]
chop1=where(wave ge 180.0 and wave le 194.5)
chop2=where(wave gt 194.5 and wave le 220.0)
g1=poly(wave[chop1]*10,s4)
g2=poly(wave[chop2]*10,s3)
gauss=[g1,g2]
endif


pxval=[wave[chop[0]],wave[chop],wave[chop[n-1]]]
pyval=[0,gauss,0]

;polyfill, wave[chop],gauss, color=180 ; plot the laser spectrum ORIG

polyfill, pxval,pyval, color=180 ; plot the laser spectrum


plot, wl,s32,xrange=xri,yrange=yri,xtitle='Nanometers',charthick=ct,charsize=cs,xthick=xyth,ythick=xyth,thick=th,ytitle='fractionation factor (permil)',/xstyle,ystyle=ystyle,lines=lines[0],psym=psym[0], /nodata,title=title,xmargin=[10,6], /noerase
loadct,13

endif  ;end int loop






;temp testing of directly shifted 32SO2
;HR2new33 = shift(HR2new32,1)
;HR2new34 = shift(HR2new32,2)



e33=1000.*(HR2new33/HR2new32 -1.)
e34=1000.*(HR2new34/HR2new32 -1.)
capE33=e33 - 1000.*((1+e34/1000.)^(0.515)-1)

if abscicon eq 0 then oplot, wave, e34, thick=th, psym=psym[0], color=colors[0]

oplot, wave, capE33, thick=th, psym=psym[1], color=colors[1]

n=n_elements(CapE33)-1

plots, xri,[0,0], thick=th+4, color=colors[1]

if abscicon eq 0 then begin
label=['!U34!N!9e!Dl!N!3','!U33!N!9E!Dl!N!3']
   legend,label, thick=th+1,charthick=ctl,charsize=csl, /top,/right,lines=lines,psym=psym,color=colors,textcolor=colors,spacing=1.5

endif else begin

if int eq 1 then begin
label=['!U33!N!9E!Dl!N!3','laser']
if broadband then label=['!U33!N!9E!Dl!N!3','Masterson lamp']

   legend,label, thick=th+1,charthick=ctl,charsize=csl,lines=[lines[1],lines[1]],psym=[psym[1],0],color=[colors[1],0],textcolor=[colors[1],0],spacing=1.5,pos=[193.29,-300]

endif else begin
label=['!U33!N!9E!Dl!N!3']
   legend,label, thick=th+1,charthick=ctl,charsize=csl, /bottom,/right,lines=lines[1],psym=psym[1],color=colors[1],textcolor=colors[1],spacing=1.5
endelse
endelse ;end abscicon loop


if int then begin
;chop=where(wave ge 192.5 and wave le 193.5)
;n=n_elements(chop)
;gauss=exp(-1*((findgen(n)-float(n)/2)/(float(n)/6))^2)
s32=HR2new32[chop]
n32=s32/max(s32)

print, int_tabulated(wave[chop],CapE33[chop]*gauss)

;print, int_tabulated(wave[chop],CapE33[chop]*gauss*n32)

axis, yaxis=1,/save,ythick=xyth,charthick=ct,charsize=cs,yrange=[0,1], ystyle=1, ytitle='laser transmission function'

oplot, wave[chop],gauss, thick=th
;oplot, wave[chop],n32, thick=th


yrange=[-400,400]
plot, wave[chop],gauss*CapE33[chop],xtitle='Nanometers',charthick=ct,xthick=xyth,ythick=xyth,thick=th,ytitle='fractionation factor*tranmission (permil)',/xstyle,/ystyle,lines=lines[0],psym=psym[0],title=title,charsize=cs,xmargin=[10,3],yrange=yrange

;plot, wave[chop],gauss*CapE33[chop]*n32,xtitle='Nanometers',charthick=ct,xthick=xyth,ythick=xyth,thick=th,ytitle='fractionation factor*tranmission (permil)',/xstyle,/ystyle,lines=lines[0],psym=psym[0],title=title,charsize=cs,xmargin=[10,3],yrange=yrange


plots, [192.5,193.5],[0,0], thick=th+4, color=colors[1]
xyouts, 193.0,300, ' !9D!3!U33!NS fractionation factors',charthick=cth,charsize=csl
xyouts, 193.0,250, 'convolved with 1 nm ',charthick=cth,charsize=csl
xyouts, 193.0,200, 'Gaussian about 193 nm',charthick=cth,charsize=csl
xyouts, 193.0,150, 'Integrated !9D!3!U33!NS:',charthick=cth,charsize=csl
xyouts, 193.0,80, '       + 18.8 permil',charthick=cth,charsize=csl+0.5

start=findgen(60)/2 + 183.5  ;starting array from 183.5 to 213.5

count=0
for i=0,n_elements(start)-1 do begin
chop=where(wave ge start[i] and wave le start[i]+1)
n=n_elements(chop)
gauss=exp(-1*((findgen(n)-float(n)/2)/(float(n)/6))^2)
print, start[i],start[i]+1,int_tabulated(wave[chop],CapE33[chop]*gauss),int_tabulated(wave[chop],CapE33[chop]*gauss*n32)

if start[i] eq 192 then begin

;repeat first plot
ystyle=9

plot, wl,s32,xrange=xri,yrange=yri,xtitle='Nanometers',charthick=ct,charsize=cs,xthick=xyth,ythick=xyth,thick=th,ytitle='fractionation factor (permil)',/xstyle,ystyle=ystyle,lines=lines[0],psym=psym[0], /nodata,title=title,xmargin=[10,6]

axis, yaxis=1,/save,ythick=xyth,charthick=ct,charsize=cs,yrange=[0,1], ystyle=1, ytitle='laser transmission function'

loadct,0
polyfill, wave[chop],gauss, color=180
oplot, wave[chop],gauss, thick=th
plot, wl,s32,xrange=xri,yrange=yri,xtitle='Nanometers',charthick=ct,charsize=cs,xthick=xyth,ythick=xyth,thick=th,ytitle='fractionation factor (permil)',/xstyle,ystyle=ystyle,lines=lines[0],psym=psym[0], /nodata,title=title,xmargin=[10,6], /noerase
loadct,13

oplot, wave, capE33, thick=th, psym=psym[1], color=colors[1]
plots, xri,[0,0], thick=th+4, color=colors[1]

label=['!U33!N!9E!Dl!N!3','laser']
   legend,label, thick=th+1,charthick=ctl,charsize=csl,lines=[lines[1],lines[1]],psym=[psym[1],0],color=[colors[1],0],textcolor=[colors[1],0],spacing=1.5,/bottom, /right


yrange=[-300,300]
plot, wave[chop],gauss*CapE33[chop],xtitle='Nanometers',charthick=ct,xthick=xyth,ythick=xyth,thick=th,ytitle='fractionation factor*tranmission (permil)',/xstyle,/ystyle,lines=lines[0],psym=psym[0],title=title,charsize=cs,xmargin=[10,3],yrange=yrange

plots, [192.0,193.0],[0,0], thick=th+4, color=colors[1]
xyouts, 192.3,250, ' !9D!3!U33!NS fractionation factors',charthick=cth,charsize=csl
xyouts, 192.3,210, 'convolved with 1 nm ',charthick=cth,charsize=csl
xyouts, 192.3,170, 'Gaussian about 192.5 nm',charthick=cth,charsize=csl
xyouts, 192.3,130, 'Integrated !9D!3!U33!NS:',charthick=cth,charsize=csl
xyouts, 192.3,70, '     - 28.8 permil',charthick=cth,charsize=csl+0.5




endif


endfor

chop=where(wave ge 183.5 and wave le 220.)
print, 183.5,220.,int_tabulated(wave[chop],CapE33[chop])


endif



endif ;end plot 5

if plot eq 6 then begin
colors=[80,160,240]
psym=[4,1,5]
lines=[0,0,0]

xri=[46600,47600]
yri=[1e-19,1e-17]

;xri=[32700,33700]
;yri=[1e-19,1.1e-18]

ylog=0

title='Danielache data'


plot, wn,s32,xrange=xri,yrange=yri,ylog=ylog,xtitle='wavenumbers',charthick=ct,xthick=xyth,ythick=xyth,thick=th,ytitle='cross section (cm!U2!N)',/xstyle,/ystyle,lines=lines[0],psym=psym[0], /nodata,title=title

oploterror, wn, s32, d32,thick=th, color=colors[0],lines=0,errthick=th-3,errcolor=colors[0]
;oploterror, wn, s32, thick=th, psym=psym[0], color=colors[0]

oploterror, wn, s33, d33,thick=th, color=colors[1],lines=0,errthick=th-3,errcolor=colors[1]
;oploterror, wn, s33, thick=th, psym=psym[1], color=colors[1]


oploterror, wn, s34, d34, thick=th, color=colors[2],lines=0,errthick=th-3,errcolor=colors[2]
;oploterror, wn, s34, thick=th, psym=psym[2], color=colors[2]



label=['!U32!NSO!D2!N','!U332!NSO!D2!N','!U34!NSO!D2!N']

   legend,label, thick=th+1,charthick=ctl,charsize=csl, /bottom,/left,lines=lines,psym=psym,color=colors,textcolor=colors,spacing=1.5

endif ;end plot 6







device, /inches, xsize=7, ysize=5
device, /close
set_plot, 'x'
!p.font=0
littlesalmonspawn='gv '+epsname+'&'
spawn, littlesalmonspawn




end
