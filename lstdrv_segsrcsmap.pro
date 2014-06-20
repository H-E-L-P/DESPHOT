PRO lstdrv_segsrcSMAP, sx,sy,pfwhm250,pfwhm350,pfwhm500,hdr250, map250,sig250 $
                   ,hdr350, map350,sig350 $
                   ,hdr500, map500,sig500,segsrc,realin=realin,dd=dd
                        
     



prf250=psf_gaussian(fwhm=pfwhm250,npix=[11,11])

if(keyword_set(realin)) then prf250=prf250/max(prf250)
; build blend map
naxis1=sxpar(hdr250,'naxis1')
naxis2=sxpar(hdr250,'naxis2')

bim=fltarr(naxis1,naxis2)
print,'build blend map'

; 250 micron SN map
psize250=sxpar(hdr250,'cd2_2')
bim[sx,sy]=1.
flat=fltarr(49,49)+1.
cb250=convol(bim,flat)
bad=where(cb250 gt 1)
cb250[bad]=1.
;sm250=cb250*convol(map250/sqrt(sig250^2+dd^2),prf250);cb250*map250/sig250/total(prf250^2)
sm250=cb250*map250/(sqrt(sig250^2+dd^2))
blank250=where(sm250 le 0.5, nblank250,comp=source250)
sm250[blank250]=-100.
smap250=intarr(naxis1,naxis2)
smap250[blank250]=-1
smap250[source250]=0


;350 micron SN map
; convert to 500 micron parameters

psize350=sxpar(hdr350,'cd2_2')

prf350=psf_gaussian(fwhm=pfwhm350,npix=[11,11],/norm)

;cb350=convol(bim,prf350)
sm350=map350/sqrt(sig350^2+dd^2);cb350*cmap350/csig350/total(prf350^2)

hastrom,sm350,hdr350,smap350,chdr350,hdr250,interp=1
smap350=smap350*cb250
blank350=where(smap350 le 0.5, nblank350,comp=source350)
smap350[blank350]=-1.
smap350[source350]=0

;500 micron SN map
; convert to 500 micron parameters

psize500=sxpar(hdr500,'cd2_2')

prf500=psf_gaussian(fwhm=pfwhm500,npix=[11,11],/norm)
;cb500=convol(bim,prf500)
sm500=map500/sqrt(sig500^2+dd^2)
;sm500=cb*cmap500/csig500/total(prf500^2)
hastrom,sm500,hdr500,smap500,chdr500,hdr250,interp=1
smap500=smap500*cb250
blank500=where(smap500 le 0.5, nblank500,comp=source500)
smap500[blank500]=-1
smap500[source500]=0




smap=smap250*smap350*smap500
lim=where(smap ne 0)
smap[lim]=-1
;stop
tacc=1
snum=0L

; flood fill each segment
while tacc ne 0 do begin
   ; find first unallocated pixel
   tp=where(smap eq 0,tacc)
   ; fill that pixel with current segment
   ; number
   snum=snum+1L
;   smap[x_pix[tp[0]],y_pix[tp[0]]]=snum
  
; try and fill surrounding pixels
   sp=[tp[0]]
   fill=where(smap[sp] eq 0.,fm)
   k=0
   stnum=snum
   while fm ne 0 do begin
      ; check for mergers
      merge=where(smap[sp] ne stnum and smap[sp] ne 0 and smap[sp] ne -1)
      if(merge[0] ne -1) then begin
      remark=where(smap eq stnum)
      stnum=smap[sp[merge[0]]]
      smap[remark]=stnum
      endif
      smap[sp[fill]]=stnum
      ;sp=[sp[0]-1,sp,sp[n_elements(sp)-1]+1]
      ;sp[1]=sp[1]-1
      if(smap[sp[0]] ne -1)then sp=[sp[0]-1,sp]
      if(smap[sp[n_elements(sp)-1]] ne -1) then sp=[sp,sp[n_elements(sp)-1]+1]
      fill=where(smap[sp] eq 0 ,fm)      
      
     ;  icplot,smap[200:250,200:250]
    ;   wait,0.2
       if fm eq 0 then begin
          
          good=where(smap[sp] eq stnum)
          sp=sp[good]+naxis1
          
          fill=where(smap[sp] ge 0,fm) 
      
       endif
      
    endwhile
  ; stop 
endwhile
; find segment for each source
n_src=n_elements(sx)
segsrc=lonarr(n_src)
for i=0L,n_src-1 do begin
    segsrc[i]=smap[sx[i],sy[i]]
    rad=1
    while(segsrc[i] le 0) do begin
        xrad=sx[i]+lindgen(2*rad+1)-rad
        yrad=sy[i]+lindgen(2*rad+1)-rad
        make_2d,xrad,yrad,dx,dy
        tmap=smap[dx,dy]
       tst=where(tmap gt 0)
       if(tst[0] ne -1) then segsrc[i]=tmap[tst[0]]
       rad=rad+1
   endwhile

endfor

end
