
readcol, '../../GRIDS/zahnle.grid',wl,wu
openw,xxx,'claire.grid', /get_lun

nworig=n_elements(wl)

for i=0,nworig-1 do begin

if wl(i) ge 1786. and wl(i) le 2222. then begin
inc=wu(i)-wl(i)  ;inc is how many angstromgs between points in the grid
print, i,wl(i), inc

start=wl(i)
finish=wu(i)
gridsize=0.5

for j=start,finish-gridsize,gridsize do begin
printf,xxx, j,j+gridsize
endfor



endif else begin
printf, xxx, wl(i),wu(i)
print, i,wl(i), wu(i)-wl(i)
endelse



endfor


close, xxx
free_lun, xxx
end
