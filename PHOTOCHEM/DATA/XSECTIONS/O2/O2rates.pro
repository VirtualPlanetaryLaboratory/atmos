epsname='o2prates.eps'

    th=5 & xyth = 8 & ct=4 & cs=1.5 
    ctl=4 & csl=1.25
    set_plot, 'ps'
    !p.font=0
    device, encapsul=1, /helv,/isolatin1, filename=epsname
    device, /color

loadct, 13

;option 1 = modern atmosphere (from October 2009)
;option 2 = ancient atmosphere (from October 2009)
;option 3 = early archean atm (Nov 2010) - Thuillier flux
;option 4 = late archean atm (Nov 2010) - photo.dat flux
;option 5 = late archean atm (Nov 2010) - Thuillier flux

option=5
;it would be worthwhile to re-run experiment 1 using the Thuillier flux corrections



if option eq 1 then begin
title='O2 -> O+O photorate - Modern atmosphere'
xri=[1e-30,1e-6]
yri=[0,80]
n1='modIO20'
n2='modIO21'
n3='modIO22'
n4='modIO22newConv'
label=['Allen Band Model','JFK exp sums','Yoshino2OldGrid','Yoshino2NewGrid']
endif

if option eq 2 then begin
title='O2 -> O+O photorate - Ancient atmosphere'
;I have no memory of what this atmosphere was, but it was probably
;late Archean
xri=[1e-12,1e-7]
yri=[0,80]
n1='ArchIO20'
n2='ArcheanIO21'
n3='ArchIO22OG_2'  ;redone with full O2 xs added to new SR
n4='ArchIO22NG_2' ;redone with full o2 xs added to new SR

label=['Allen Band Model (bad)','JFK exp sums','Yoshino2OldGrid','Yoshino2NewGrid']
endif


if option eq 3 then begin
;ok Nov 2010 versions (no O2, CH4=1e-4, CO2=3e-3, timega=0, Thuillier
;flux, 100 km grid


title='O2 -> O+O photorate - low O2 atmosphere'
xri=[1e-10,1e-7]
yri=[0,100]

n1='newArchIO21'
n2='newArchIO22_oldgrid'
n3='newArchIO22_0.5'
n4='newArchIO22_0.1'
n5='newArchIO22_0.05'

label=['JFK exp sums','Yoshino2OldGrid','Yoshino2_0.5Agrid','Yoshino2_0.1Agrid','Yoshino2_0.05Agrid']
endif


if option eq 4 then begin
;ok Nov 2010 versions (O2=1e-8, CH4=1e-4, CO2=3e-3, timega=0, photo.dat
;flux, 100 km grid


title='O2 -> O+O photorate - Late Archean atm - photo.dat flux'
xri=[1e-10,1e-7]
yri=[0,100]

;should check this against the new results from option5. these should
;be the same now... yes they are.  this is good

n1='newLArchIO21'
n2='newLArchIO22_oldgrid'
n3='newLArchIO22_0.5'
n4='newLArchIO22_0.1'
n5='newLArchIO22_0.05' ; ok - so 0.1A is fine...

label=['JFK exp sums','Yoshino2OldGrid','Yoshino2_0.5Agrid','Yoshino2_0.1Agrid','Yoshino2_0.05Agrid']
endif


if option eq 5 then begin
;ok Nov 2010 versions (O2=1e-8, CH4=1e-4, CO2=3e-3, timega=0,
;Thuillier flux, 100 km grid

;above was initial experiment

;next are fixed flux models need to produce O2=2e-9 and CH41e-4:
;fluxes: O2=5.6E12, CH4, 3.221e12
;going to better flux drops these significantly...

;how do the mixing ratio's change with different cross
;sections: (OLD METHOd/NEW METHOD)

;0.1A grid O2=2e-10, Ch4=5.9e-5
;0.05A grid O2=1.15e-10, CH4=4.13E-5
;0.02A grid O2=7e11, CH4=2.83e5
;oldgrid O2=5.e8, CH4=8.2e04
;IO2=1  O2=5.1e8, CH4=8.2e-4


title='O2 -> O+O photorate - Late Archean atm - Atlas1 flux'
xri=[1e-10,1e-7]
yri=[0,100]

;; original (incorrect) tests with faulty binning
;; n1='newLArchIO21_highFLUX'
;; ;n1='newLArchIO22_0.02_highFLUX'

;; n2='newLArchIO22_oldgrid_highFLUX'
;; n3='newLArchIO22_0.5_highFLUX'
;; n4='newLArchIO22_0.1_highFLUX'
;; n5='newLArchIO22_0.05_highFLUX'


;ok here we go again:
;using high res flux which is now properly binned
;starting with LR model (i.e. exponential sums)
;ok Nov 2010 versions (O2=2e-9, CH4=1e-4, CO2=3e-3, timega=0,
;Thuillier flux, 100 km grid
;this requires:  CH4 flux = 2.087350E+11  O2 flux = 2.580970E+11
;oldgrid
;0.5:            CH4 flux = 1.982749E+11  O2 flux = 2.542530E+11
;0.1             CH4 flux = 1.827849E+11  O2 flux = 2.354877E+11
;0.05
;0.02            CH4 flux = 1.643669E+11  O2 flux = 2.120638E+11


;ok all is looking well.  going to repeat the #3 experiment

;new tests - uncomment when completed (fixed mr runs)

n1='newLArchIO21_TFLUX'
;n1='newLArchIO22_0.02_TFLUX'

n2='newLArchIO22_oldgrid_TFLUX'
n3='newLArchIO22_0.5_TFLUX'
n4='newLArchIO22_0.1_TFLUX'
n5='newLArchIO22_0.05_TFLUX'



label=['JFK exp sums','Yoshino2OldGrid','Yoshino2_0.5Agrid','Yoshino2_0.1Agrid','Yoshino2_0.05Agrid']
;label=['Yosh_0.02Agrid','Yoshino2OldGrid','Yoshino2_0.5Agrid','Yoshino2_0.1Agrid','Yoshino2_0.05Agrid']
endif

orig =0  ; 0 is production loop for paper where I change around the color scheme...)

if orig eq 1 then begin


colors=[0,80,160,240]
psym=[1,5,2,4]

if option ge 3 then begin
colors=[colors,35]
psym=[psym,6]
endif

readcol, n1, z,po1d,po2
z=z/1e5

plot, po2,z,xrange=xri,yrange=yri,/xlog,xtitle='photolysis rate',charthick=ct,xthick=xyth-2,ythick=xyth-2,thick=th-3,ytitle='height (km)',/xstyle,psym=psym[0],title=title

readcol,n2,z,po1d1,po21
z=z/1e5

oplot, po21, z, thick=th, psym=psym[1], color=colors[1]


readcol,n3,z,po1d2,po22
z=z/1e5

oplot, po22, z, thick=th, psym=psym[2], color=colors[2]


readcol,n4,z,po1d3,po24
z=z/1e5

oplot, po24, z, thick=th, psym=psym[3], color=colors[3]


if option ge 3 then begin
readcol,n5,z,po1d5,po25
z=z/1e5
oplot, po25, z, thick=th, psym=psym[4], color=colors[4]
endif


   legend,label, thick=th+1,charthick=ctl,charsize=csl, /top,/left,lines=0,psym=psym,color=colors,textcolor=colors,spacing=1.5
endif else begin
;start production loop here


colors=[80,225,35,256,0]
psym=[5,2,4,1,0]

title=''
label=['Kasting exponential sum','Yoshino to original grid','Yoshino to 0.5A grid','Yoshino to 0.1A grid','Yoshino to 0.05A grid']


readcol, n1, z,po1d,po2
z=z/1e5

plot, po2,z,xrange=xri,yrange=yri,/xlog,xtitle='O!D2!N photolysis rate (s!U-1!N)',charthick=ct,xthick=xyth,ythick=xyth,thick=th,ytitle='height (km)',/xstyle,psym=psym[0],title=title, /nodata,charsize=cs


oplot, po2,z,thick=th, psym=psym[0],color=colors[0]


readcol,n2,z,po1d1,po21
z=z/1e5

oplot, po21, z, thick=th, psym=psym[1], color=colors[1]


readcol,n3,z,po1d2,po22
z=z/1e5

oplot, po22, z, thick=th, psym=psym[2], color=colors[2]


readcol,n4,z,po1d3,po24
z=z/1e5

oplot, po24, z, thick=th, psym=psym[3], color=colors[3]


readcol,n5,z,po1d5,po25
z=z/1e5
oplot, po25, z, thick=th, psym=psym[4], color=colors[4]




   legend,label, thick=th+1,charthick=ctl,charsize=csl, /top,/left,lines=0,psym=psym,color=colors,textcolor=colors,spacing=1.5

xyouts, 1.25e-10,5., 'B', charthick=cth, charsize=cs


endelse


device, /inches, xsize=7, ysize=5
device, /close
set_plot, 'x'
!p.font=0
littlesalmonspawn='gv '+epsname+'&'
spawn, littlesalmonspawn




end
