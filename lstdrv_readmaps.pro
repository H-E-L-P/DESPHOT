

pro lstdrv_readmaps,pswfits,pmwfits,plwfits,sra,sdec,prfsize,fcat,f_src=f_src,name=name,astromcorr=astromcorr,fcatmod=fcatmod,noseg=noseg,mralims=mralims,mdeclims=mdeclims,segthresh=segthresh,smap=smap,sarea=sarea
;
;----------------------------------------------------------------------
; parameter checking & defaults
;----------------------------------------------------------------------
; source parameter

;image file

;output file names
if not keyword_set(name) then name="testName"
; PRF parameters
if not keyword_set(segthresh) then segthresh=0.5

; OK read in fits image and header and list of srcs from fits file
; read in image and naxis 
;blank=mrdfits(pswfits,0,bhdr250)
im250=mrdfits(pswfits,1,hdr250)*1e3
im350=mrdfits(pmwfits,1,hdr350)*1e3
im500=mrdfits(plwfits,1,hdr500)*1e3


naxis1_250=sxpar(hdr250,"naxis1")
naxis2_250=sxpar(hdr250,"naxis2")
naxis1_350=sxpar(hdr350,"naxis1")
naxis2_350=sxpar(hdr350,"naxis2")
naxis1_500=sxpar(hdr500,"naxis1")
naxis2_500=sxpar(hdr500,"naxis2")
if(keyword_set(smap)) then begin
pixsize250=3600.*sxpar(hdr250,"cd2_2")
pixsize350=3600.*sxpar(hdr350,"cd2_2")
pixsize500=3600.*sxpar(hdr500,"cd2_2")
endif else begin
pixsize250=3600.*sxpar(hdr250,"cdelt2")
pixsize350=3600.*sxpar(hdr350,"cdelt2")
pixsize500=3600.*sxpar(hdr500,"cdelt2")

endelse


if(keyword_set(astromcorr)) then begin
oldcrval1=sxpar(hdr250,'crval1')
fxaddpar,hdr250,'crval1',(oldcrval1+astromcorr[0]/3600.)
oldcrval2=sxpar(hdr250,'crval2')
fxaddpar,hdr250,'crval2',(oldcrval2+astromcorr[1]/3600.)
oldcrval1=sxpar(hdr350,'crval1')
fxaddpar,hdr350,'crval1',(oldcrval1+astromcorr[0]/3600.)
oldcrval2=sxpar(hdr350,'crval2')
fxaddpar,hdr350,'crval2',(oldcrval2+astromcorr[1]/3600.)
oldcrval1=sxpar(hdr500,'crval1')
fxaddpar,hdr500,'crval1',(oldcrval1+astromcorr[0]/3600.)
oldcrval2=sxpar(hdr500,'crval2')
fxaddpar,hdr500,'crval2',(oldcrval2+astromcorr[1]/3600.)

endif

make_2d,findgen(naxis1_250),findgen(naxis2_250),x_pix250,y_pix250

pixsize=[pixsize250,pixsize350,pixsize500]

if(keyword_set(smap)) then begin
nim250=mrdfits(pswfits,2,nhdr250)*1e3
nim350=mrdfits(pmwfits,2,nhdr350)*1e3
nim500=mrdfits(plwfits,2,nhdr500)*1e3
endif else begin
nim250=mrdfits(pswfits,6,nhdr250)*1e3
nim350=mrdfits(pmwfits,6,nhdr350)*1e3
nim500=mrdfits(plwfits,6,nhdr500)*1e3


endelse
adxy,hdr250,sra,sdec,sx250,sy250
; filter out sources not in the map, just check 250

sgood=where((sx250 gt 0) and (sx250 lt naxis1_250) and (sy250 gt 0) and (sy250 lt naxis2_250) and finite(im250[sx250,sy250]),ngood)

if(sgood(0) eq -1) then  stop
sx250=sx250[sgood]
sy250=sy250[sgood]
sra=sra[sgood]
sdec=sdec[sgood]
f_src=f_src[sgood]

n_src=ngood

adxy,hdr350,sra,sdec,sx350,sy350
adxy,hdr500,sra,sdec,sx500,sy500





bad=where(finite(im250) eq 0 or finite(nim250) eq 0 or nim250 eq 0)
if(bad[0] ne -1) then begin
im250[bad]=0.
nim250[bad]=1.
endif
bad=where(finite(im350) eq 0 or finite(nim350) eq 0 or nim250 eq 0)
if(bad[0] ne -1) then begin
im350[bad]=0.
nim350[bad]=1.
endif
bad=where(finite(im500) eq 0 or finite(nim500) eq 0 or nim250 eq 0)
if(bad[0] ne -1) then begin
im500[bad]=0.
nim500[bad]=1.
endif













; work out pseudo p-stat to weight models by
; calculate prob of source with f>fs in one beam
; work out input list area crudely by just assuming contiguous area
; from min to max ra & dec of sources
p_src=fltarr(n_src)

for k=0L,n_src-1 do begin
   tmp=where(f_src ge f_src[k],ntmp)
   ns=ntmp/sarea
   p_src[k]=!pi*(prfsize[0]/3600.)^2/4./alog(2.)*ns
endfor


;---------------
; load PRF from file
;---------------


pfwhm=prfsize/pixsize
sim250=im250;+3.2;+0.4
sim350=im350;+4.3;+0.7
sim500=im500;+5.7;+1.
segsrc=lonarr(n_src)+2

if not keyword_set(noseg) then lstdrv_segsrc,sx250,sy250,pfwhm[0],pfwhm[1],pfwhm[2],hdr250,im250,nim250,hdr350,im350,nim350,hdr500,im500,nim500,segsrc,bsn=segthresh
;stop

msegsize=max(histogram(segsrc))
if (msegsize gt 2000) then stop 
;segsrc=fltarr(n_src)+2
fstruct={xid:0L,inra:0.,indec:0.,f250:0.,e250:0.,et250:0.,chi250:0.,f350:0.,e350:0.,et350:0.,chi350:0.,f500:0.,e500:0.,et500:0.,chi500:0.,gid:0L,gsize:0,bkg250:0.,bkg350:0.,bkg500:0.}
fcat=replicate(fstruct,n_src)
fcat.inra=sra
fcat.indec=sdec
fcat.xid=long(sgood)
test250=1e10


tdiff250=1.
tdiff350=1.
tdiff500=1.
test250=1e4
test350=1e4
test500=1e4
;stop
;while(mean(abs([tdiff250,tdiff350,tdiff500])) gt 0.05) do begin
lstdrv_solvefluxes2,segsrc,sx250,sy250,sx350,sy350,sx500,sy500,pfwhm,sim250,sim350,sim500,nim250,nim350,nim500,naxis1_250,naxis2_250,naxis1_350,naxis2_350,naxis1_500,naxis2_500,p_src,fcat=fcat,reconmap=reconmap,f_src=f_src
oldtest250=test250
oldtest350=test350
oldtest500=test500


back250=sim250-reconmap.rmap250
back350=sim350-reconmap.rmap350
back500=sim500-reconmap.rmap500


;linfit,reconmap.rmap250,sim250,a,b
;kern=psf_gaussian(fwhm=15.,npix=[50,50],/norm);
;back250=convol(back250,kern)
;back350=convol(back350,kern)
;back500=convol(back500,kern)

;back250=smooth(back250,[21,21])
;back350=smooth(back350,[21,21])
;back500=smooth(back500,[21,21])


paxis=[11,11]
prf=psf_gaussian(fwhm=pfwhm[0],npix=paxis)
bmask250=fltarr(naxis1_250,naxis2_250)
make_2d,indgen(11)-5,indgen(11)-5,tx,ty
for k=0,n_src-1 do bmask250[sx250[k]+tx, sy250[k]+ty]=bmask250[sx250[k]+tx, sy250[k]+ty]+prf
;region250=fltarr(naxis1_250,naxis2_250)
;for k=0,n_src-1 do region250[sx250[k]+tx, sy250[k]+ty]=1.
region250=where(bmask250 ne 0.)

test250=linfit(bmask250[region250],(im250-reconmap.rmap250)[region250])

tdiff250=test250-oldtest250

prf=psf_gaussian(fwhm=pfwhm[1],npix=paxis)
bmask350=fltarr(naxis1_350,naxis2_350)
make_2d,indgen(11)-5,indgen(11)-5,tx,ty
for k=0,n_src-1 do bmask350[sx350[k]+tx, sy350[k]+ty]=bmask350[sx350[k]+tx, sy350[k]+ty]+prf
;testx=min(sx350)+findgen(max(sx350)-min(sx350)-15)+7
;testy=min(sy350)+findgen(max(sy350)-min(sy350)-15)+7
;make_2d,testx,testy,rx,ry
region350=where(bmask350 ne 0.)
test350=linfit(bmask350[region350],(im350-reconmap.rmap350)[region350])
tdiff350=test350-oldtest350


prf=psf_gaussian(fwhm=pfwhm[2],npix=paxis)
bmask500=fltarr(naxis1_500,naxis2_500)
make_2d,indgen(11)-5,indgen(11)-5,tx,ty
for k=0,n_src-1 do bmask500[sx500[k]+tx, sy500[k]+ty]=bmask500[sx500[k]+tx, sy500[k]+ty]+prf
;testx=min(sx500)+findgen(max(sx500)-min(sx500)-15)+7
;testy=min(sy500)+findgen(max(sy500)-min(sy500)-15)+7
;make_2d,testx,testy,rx,ry

region500=where(bmask500 ne 0.)
test500=linfit(bmask500[region500],(im500-reconmap.rmap500)[region500])
tdiff500=test500-oldtest500

sim250=im250-test250[0]
sim350=im350-test350[0]
sim500=im500-test500[0]
print,test250,tdiff250,test350,tdiff350,test500,tdiff500


print,mean(im250-sim250),mean(im350-sim350),mean(im500-sim500)
;stop;endwhile
mwrfits,blank,name+'_back.fits',bhdr250,/create
mwrfits,back250,name+'_back.fits',hdr250
mwrfits,back350,name+'_back.fits',hdr350
mwrfits,back500,name+'_back.fits',hdr500

; find missing sources
res250=(im250-reconmap.rmap250-test250[0])/nim250
res350=(im350-reconmap.rmap350-test350[0])/nim350
res500=(im500-reconmap.rmap500-test500[0])/nim500


find,res250,nx250,ny250,nf250,sharp,round,3.,pfwhm[0],[0.,3.],[-3.,3.],/silent
find,res350,nx350,ny350,nf350,sharp,round,3.,pfwhm[1],[0.,3.],[-3.,3.],/silent
find,res500,nx500,ny500,nf500,sharp,round,3.,pfwhm[2],[0.,3.],[-3.,3.],/silent
; keep 250 list and filter out sources in 350 and 500
xyad,hdr250,nx250,ny250,nra250,ndec250
xyad,hdr350,nx350,ny350,nra350,ndec350
xyad,hdr500,nx500,ny500,nra500,ndec500

;f1=where(sqrt(((nra250-nra350)*cos(ndec250*!pi/180.))^2+(ndec250-ndec350)^2) lt 20.,comp=k1)
;f2=where(sqrt(((nra250-nra500)*cos(ndec250*!pi/180.))^2+(ndec250-ndec500)^2) lt 20.,comp=k1)

nra=nra250
ndec=ndec250

for i=0,n_elements(nx350)-1 do begin
sep=sqrt(((nra-nra350[i])*cos(ndec*!pi/180.))^2+(ndec-ndec350[i])^2)
if(min(sep) gt 20./3600.) then begin
    nra=[nra, nra350[i]]
    ndec=[ndec,ndec350[i]]
endif
endfor
for i=0,n_elements(nx500)-1 do begin
sep=sqrt(((nra-nra500[i])*cos(ndec*!pi/180.))^2+(ndec-ndec500[i])^2)
if(min(sep) gt 20./3600.) then begin
    nra=[nra, nra500[i]]
    ndec=[ndec,ndec500[i]]
endif
endfor

; check if any are within 0.5 FHWM of known sources, just check 250?
nm=-1
for i=0,n_elements(nra)-1 do begin
sep=sqrt(((nra[i]-sra)*cos(ndec[i]*!pi/180.))^2+(sdec-ndec[i])^2)
if(min(sep) gt 18./3600.) then nm=[nm,i]
endfor
nm=nm[1:n_elements(nm)-1]
nra=nra[nm]
ndec=ndec[nm]

; check within mralims and mdeclims
if(keyword_set(mralims)) then begin
gnew=where(nra gt mralims[0] and nra lt mralims[1] and ndec gt mdeclims[0] and ndec lt mdeclims[1])
nra=nra[gnew]
ndec=ndec[gnew]
endif
adxy,hdr250,nra,ndec,nx250,ny250


; check if sources are actually in the map!
newgood=where((nx250 gt 0) and (nx250 lt naxis1_250) and (ny250 gt 0) and (ny250 lt naxis2_250) and finite(im250[nx250,ny250]),nnew)

if(nnew gt 0.) then begin

nx250=nx250[newgood]
ny250=ny250[newgood]
nra=nra[newgood]
ndec=ndec[newgood]
endif
adxy,hdr350,nra,ndec,nx350,ny350
adxy,hdr500,nra,ndec,nx500,ny500


nsources={nra:nra,ndec:ndec,nx250:nx250,ny250:ny250,nx350:nx350,ny350:ny350,nx500:nx500,ny500:ny500,obs250:fltarr(nnew), e250:fltarr(nnew), obs350:fltarr(nnew), e350:fltarr(nnew), obs500:fltarr(nnew), e500:fltarr(nnew)}

sx250=[sx250,nx250]
sy250=[sy250,ny250]
sx350=[sx350,nx350]
sy350=[sy350,ny350]
sx500=[sx500,nx500]
sy500=[sy500,ny500]
p_src=[p_src,fltarr(nnew)+1.]
f_src=[f_src,fltarr(nnew)+min(f_src)]
segsrc=lonarr(n_src+nnew)+2
if not keyword_set(noseg) then lstdrv_segsrc,sx250,sy250,pfwhm[0],pfwhm[1],pfwhm[2],hdr250,im250,nim250,hdr350,im350,nim350,hdr500,im500,nim500,segsrc,bsn=segthresh
;segsrc=intarr(n_src+nnew)+2.
fcatmod=replicate(fstruct,n_src+nnew)
fcatmod.inra=[sra, nra]
fcatmod.indec=[sdec,ndec]
fcatmod.bkg250=fltarr(n_src+nnew)+test250[0]
fcatmod.bkg350=fltarr(n_src+nnew)+test350[0]
fcatmod.bkg500=fltarr(n_src+nnew)+test500[0]
nids=fltarr(nnew)-99.
fcatmod.xid=[fcat.xid,nids]
lstdrv_solvefluxes2,segsrc,sx250,sy250,sx350,sy350,sx500,sy500,pfwhm,sim250,sim350,sim500,nim250,nim350,nim500,naxis1_250,naxis2_250,naxis1_350,naxis2_350,naxis1_500,naxis2_500,p_src,fcat=fcatmod,reconmap=reconmap,f_src=f_src

; work out true error
;250
paxis=[11,11]

bmask250=fltarr(naxis1_250,naxis2_250)
make_2d,indgen(11)-5,indgen(11)-5,tx,ty
for k=0,n_src-1 do bmask250[sx250[k]+tx, sy250[k]+ty]=1.
region=where(bmask250 eq 1.)
tmap250=im250-reconmap.rmap250
prf250=psf_gaussian(fwhm=pfwhm[0],npix=paxis)
smap250=convol(tmap250,prf250)/total(prf250^2)
nm250=convol(1./nim250^2,prf250^2)
err250=stddev(smap250[region])^2-1./nm250[sx250,sy250]
good=where(err250 gt 0.)
fcatmod[good].et250=sqrt(fcatmod[good].e250^2+err250[good])
;350
bmask350=fltarr(naxis1_350,naxis2_350)
make_2d,indgen(11)-5,indgen(11)-5,tx,ty
for k=0,n_src-1 do bmask350[sx350[k]+tx, sy350[k]+ty]=1.
region=where(bmask350 eq 1.)
tmap350=im350-reconmap.rmap350
prf350=psf_gaussian(fwhm=pfwhm[1],npix=paxis)
smap350=convol(tmap350,prf350)/total(prf350^2)
nm350=convol(1./nim350^2,prf350^2)
err350=stddev(smap350[region])^2-1./nm350[sx350,sy350]
good=where(err350 gt 0)
fcatmod[good].et350=sqrt(fcatmod[good].e350^2+err350[good])

;500

bmask500=fltarr(naxis1_500,naxis2_500)
make_2d,indgen(11)-5,indgen(11)-5,tx,ty
for k=0,n_src-1 do bmask500[sx500[k]+tx, sy500[k]+ty]=1.
region=where(bmask500 eq 1.)
tmap500=im500-reconmap.rmap500
prf500=psf_gaussian(fwhm=pfwhm[2],npix=paxis)
smap500=convol(tmap500,prf500)/total(prf500^2)
nm500=convol(1./nim500^2,prf500^2)
err500=stddev(smap500[region])^2-1./nm500[sx500,sy500]
good=where(err500 gt 0) 
fcatmod[good].et500=sqrt(fcatmod[good].e500^2+err500[good])


end



































