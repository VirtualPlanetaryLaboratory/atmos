epsname='o2.ps'

    th=5 & xyth = 8 & ct=4 & cs=1.5 
    ctl=4 & csl=1.25
    set_plot, 'ps'
    !p.font=0
    device, encapsul=0, /helv,/isolatin1, filename=epsname
    device, /color

loadct, 13

;readcol, 'Yoshino92.abs', wl, xs
;save, filename='yosh.sav',wl,xs



restore, 'yosh.sav'

;colors=[0,80,170,35,300]
colors=[0,80,220,35,256]
psym=[0,5,2,4,1]

xri=[179,202]
xri=[186,189]

xri=[190,191]

xri=[188.6,190.6]


yri=[1e-23,1e-18]

plot, wl,xs,xrange=xri,yrange=yri,/ylog,xtitle='wavelength (nm)',charthick=ct,xthick=xyth,ythick=xyth,thick=th,ytitle='cross section (cm!U2!N)',/xstyle,charsize=cs


readcol,'Yoshino2oldGrid',wav2,yoshint,/silent
wav2=wav2/10 ;convert from angstroms to nm

oplot, wav2, yoshint, thick=th, psym=psym[2], color=colors[2]


;readcol, 'Yosh2NewGrid', wav3,yoshnew,/silent  ;1A grid
readcol, 'Yosh20.5Grid', wav3,yoshnew,/silent   ;0.5A grid
wav3=wav3/10.

oplot, wav3, yoshnew, thick=th, psym=psym[3], color=colors[3]


readcol, 'O2XStopofmodernATMOld',inc,wav,o1d_o,o_o,/silent
wav=wav/10.  ;convert from angstroms to nm
oplot, wav, o_o, thick=th, psym=psym[1], color=colors[1]


readcol, 'Yosh2_0.1Agrid', wav,poo,/silent
wav=wav/10.
oplot, wav, poo, thick=th, psym=psym[4], color=colors[4]





label=['Yoshino 1992 data','Allen 1982 band approx','Yoshino to original grid','Yoshino to 0.5A grid','Yoshino to 0.1A grid']

;label=['Yoshino 1992 data','Yosh2 0.1Agrid','Yosh2OldGrid','Yosh2NewGrid']


   legend,label, thick=th+1,charthick=ctl,charsize=csl, /top,/right,lines=0,psym=psym,color=colors,textcolor=colors,spacing=1.5

xyouts, 188.67,1.8e-23, 'A', charthick=cth, charsize=cs


device, /inches, xsize=7, ysize=5
device, /close
set_plot, 'x'
!p.font=0
littlesalmonspawn='gv '+epsname+'&'
spawn, littlesalmonspawn




end
