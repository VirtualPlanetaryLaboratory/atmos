; IDL profile for plotting the spectra of an AD Leo Flare
; Modified by Antigona Segura. March, 2008
;
nfile=5
mydata = '../SAVEOUT/'

openw,10,'temp.txt'
;
;       define input file namesgv
;
file_co = strarr(nfile)
file_co(0) = mydata+'GJ581dCO292SolCon1.0/output_couple.dat'
file_co(1) = mydata+'GJ581dCO295SolCon1.0/output_couple.dat'
file_co(2) = mydata+'GJ581dCO298SolCon1.0/output_couple.dat'
file_co(3) = mydata+'Gl581-CO279-O21e-2-Solcon10/output_couple.dat'
file_co(4) = mydata+'Gl581-CO292-O21e-3-Solcon10/output_couple.dat'


; Reading ascii files
; Reading values DURING the FLARE

npt0 =60
alt = fltarr(npt0,nfile)
temp = fltarr(npt0,nfile)


mincol=1
maxcol=4
nskip=354
ntskip = 1
ifm=1
idat=4

for n=0,nfile-1 do begin
    npt0=60
    print, 'reading file =',file_co(n)
   printf,10,'reading file =',file_co(n)
   printf,10,''

read_rec,file_co(n),ifm,idat,nskip,ntskip,mincol,maxcol,nx,npt0,x

    for i=8,59 do begin
        alt(i,n) = x(1,i-8)
        temp(i,n) = x(2,i-8)
       printf,10,i,alt(i,n),temp(i,n)
   endfor
printf,10,''
        for i=7,0,-1 do begin
       alt(i,n) = alt(8,n)+(i+0.5)
        temp(i,n) = temp(8,n)
       printf,10,i,alt(i,n),temp(i,n)
    endfor


printf,10,''
endfor
;
;       create the plot
;
 loadct,6
  set_plot,'ps'
;  device,/encapsulated
  device, filename='temperature.ps'
 !p.font=1
  DEVICE, SET_FONT='Helvetica', /TT_FONT
 device,/color
  device, /landscape, font_size=28
  device,/inches,ysize=7.5
  device,/inches,xsize=10.0

;  !ymargin=[24,0]
;  !p.multi=[0,2,4,0,0]
;
;       create a linear plot - 2 plots in a page
;
linethick=6
!x.thick=linethick
!y.thick=linethick
!p.thick=linethick

;multiplot,[2,4]
  plot,[130, 300],[0, 30.],/nodata,$ 
            xstyle=1,ystyle=1,$
            title='GJ581d, 0.20 AU, P!ds!n=7.5 bar',$
            xtitle='Temperature (K)',$
            ytitle='Altitude (km)'
     oplot,temp(*,0),alt(*,0),thick=4,linestyle=0, color= 0
     oplot,temp(*,1),alt(*,1),thick=4,linestyle=0, color= 40
      oplot,temp(*,2),alt(*,2),thick=4,linestyle=0, color=80
      oplot,temp(*,3),alt(*,3),thick=4,linestyle=0, color=120 
      oplot,temp(*,4),alt(*,4),thick=4,linestyle=0, color=160
   xyouts,240,5,'0.92 CO!d2', charsize=1.,color=0
   xyouts,240,10,'0.95 CO!d2', charsize=1.,color=40
   xyouts,240,15,'0.98 CO!d2', charsize=1.,color=80
   xyouts,240,20,'0.79 CO!d2!n, 1x10!u-2!nO!d2', charsize=1.,color=120
   xyouts,240,25,'0.92 CO!d2!n, 1x10!u-3!nO!d2', charsize=1.,color=160

     

device, /close_file
;   multiplot,/reset
;!P.Multi = 0
;
end
