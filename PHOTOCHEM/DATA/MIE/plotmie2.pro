pro plotmie2, i

;taken from readmie.pro, which was used to create datafiles for the
;photochemical model, as well as makes plots as a function of all
;particle wavelengths.  Here we are going hack this down to make plots
;just for 0.2 microns for presentation in the paper.


;plotmie2 is turning plotmie into a multiplot


;Sept 2010 - using this to compare the old and new MIE files for
;            fractal haze



psname='test.ps'
    th=5 & xyth = 8 & ct=4 & cs=1.5 
    ctl=4 & csl=1.25
    set_plot, 'ps'
    !p.font=0
    device, encapsul=0, /helv,/isolatin1, filename=psname
    device, /color

loadct, 12

;spherical
      root='fitmythol'
      sfilenames=['0001.DAT','0002.DAT','0003.DAT','0004.DAT','0005.DAT','0006.DAT','0007.DAT','0008.DAT','0009.DAT','001.DAT ','003.DAT ','005.DAT ','007.DAT ','01.DAT  ','013.DAT ','015.DAT ','017.DAT ','02.DAT  ','023.DAT ','027.DAT ','03.DAT  ','033.DAT ','037.DAT ','04.DAT  ','043.DAT ','047.DAT ','05.DAT  ','055.DAT ','06.DAT  ','07.DAT  ','08.DAT  ','09.DAT  ','1.DAT   ','2.DAT   ']


;fractal
      root2='fractal/fractopts'    
      ffilenames=['0.001um.txt','0.002um.txt','0.003um.txt','0.004um.txt','0.005um.txt','0.006um.txt','0.007um.txt','0.008um.txt','0.009um.txt','0.010um.txt','0.030um.txt','0.050um.txt','0.070um.txt','0.100um.txt','0.130um.txt','0.150um.txt','0.170um.txt','0.200um.txt','0.230um.txt','0.270um.txt','0.300um.txt','0.330um.txt','0.370um.txt','0.400um.txt','0.430um.txt','0.470um.txt','0.500um.txt','0.550um.txt','0.600um.txt','0.700um.txt','0.800um.txt','0.900um.txt','1.000um.txt','2.000um.txt']





omega=fltarr(34,108,2)
qext=omega
asym=omega


for j=0,33 do begin
readcol, root+strtrim(sfilenames[j],2), wav1,wav2,a,b,c
omega[j,*,0]=a
qext[j,*,0]=b
asym[j,*,0]=c

readcol, root2+strtrim(ffilenames[j],2), wav1,wav2,a,b,c
omega[j,*,1]=a
qext[j,*,1]=b
asym[j,*,1]=c
endfor

wav1=wav1/10 ;convert angstroms to nm

xrange=[150,800]
xt='wavelength (nm)'
xt='wavelength (!9m!3m)'
xout1=100.


device, xsize=14.06, ysize=23.7
!y.margin=[5,3]
!x.margin=[11,3] 



multiplot, [1,1],/init



multiplot, [1,3]

;need to add in y titles and y ticknames.  also mess with margins. and
;inset a,b,c

;ytn=['10!U-2!D','10!U-1!D','10!U0!D']

plot, wav1, qext[i,*,0],  xthick=xyth,ythick=xyth,charthick=ct,charsize=cs,thick=th, /xlog, /ylog,yrange=[0.01,50],xrange=xrange, /xstyle, /ystyle,ytitle='Q!DEXT!N',title='Particle radii = 0.2 !9m!3m  ACK ACK',ytickname=ytn 

;oplot, wav1, qext[i,*,1], psym=2, color=80,thick=th
oplot, wav1, qext[i,*,1], color=80,thick=th


yout1=2
xyouts,xout1,yout1,'a',charsize=cs+.5,charthick=ct+4


multiplot  ;second plot

ytn=['10!U-4!D','10!U-3!D','10!U-2!D','10!U-1!D','10!U0!D']
;print, pizer0[0:10]
plot, wav1, omega[i,*,0], xthick=xyth,ythick=xyth,charthick=ct,charsize=cs,thick=th, /xlog,/ylog,xrange=xrange,yrange=[1e-4,2],/xstyle, /ystyle , ytickname=ytn,ytitle='Single scattering albedo'
;oplot, wav1, omega[i,*,1], psym=2, color=80,thick=th
oplot, wav1, omega[i,*,1], color=80,thick=th


yout1=0.4
xyouts,xout1,yout1,'b',charsize=cs+.5,charthick=ct+4


multiplot ;3rd - this one should have xtitle

plot, wav1, asym[i,*,0], xthick=xyth,ythick=xyth,charthick=ct,charsize=cs,thick=th, /xlog,/ylog,xtitle=xt,xrange=xrange,yrange=[1e-4,2],/xstyle, /ystyle,ytickname=ytn,ytitle='Asymmetry factor'  

;oplot, wav1, asym[i,*,1], psym=2, color=80,thick=th
oplot, wav1, asym[i,*,1], color=80,thick=th


yout1=0.4
xyouts,xout1,yout1,'c',charsize=cs+.5,charthick=ct+4

multiplot, /reset





device, /inches, xsize=7, ysize=5
device, /close
set_plot, 'x'
!p.font=0
!p.multi=0
littlesalmonspawn='gv '+psname+'&'
spawn, littlesalmonspawn


stop
end
