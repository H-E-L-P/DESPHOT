pro mkspiremap,naxis1,naxis2,sx,sy,flux,prf,smap


; build residual map

paxis=size(prf)
n_src=n_elements(flux)
smap=fltarr(naxis1,naxis2)
make_2d,indgen(paxis[1])-paxis[1]/2,indgen(paxis[2])-paxis[2]/2,tx,ty
for k=0,n_src-1 do begin


pindx=findgen(paxis[1])-sx[k]+round(sx[k])
pindy=findgen(paxis[2])-sy[k]+round(sy[k])
make_2d,pindx,pindy,ipx2,ipy2   
nprf=interpolate(prf,ipx2,ipy2,cubic=-0.5)

smap[round(sx[k])+tx, round(sy[k])+ty]=smap[round(sx[k])+tx, round(sy[k])+ty]+nprf*flux[k]
endfor

end
