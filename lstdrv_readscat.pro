;uses a two-pass strategy to do list-driven sources extraction, also
;finds new sources in the residual maps.

;input:
; pswfits, pmwfits, plwfits:
; sra, sdec:
; prfsize:

;optional input
; name: name of the residual map
; fcatmod:
; astromcorr: astrometry correction
; mralims, mdeclims:
; segthresh:

;output

;author: Issac Roseboom
pro lstdrv_readSCAT,pswfits,pmwfits,plwfits,sra,sdec,prfsize,fcat,f_src=f_src,name=name,astromcorr=astromcorr,fcatmod=fcatmod,noseg=noseg,mralims=mralims,mdeclims=mdeclims,segthresh=segthresh
;----------------------------------------------------------------------
; parameter checking & defaults, PRF parameters
;----------------------------------------------------------------------
if not keyword_set(name) then name="testName"   ;output file names
;if not keyword_set(segthresh) then segthresh=0.5   
if not keyword_set(segthresh) then segthresh=1.0

; OK read in fits image and header and list of srcs from fits file
; read in image and naxis 
im250 = mrdfits(pswfits, 1, hdr250) * 1e3
im350 = mrdfits(pmwfits, 1, hdr350) * 1e3
im500 = mrdfits(plwfits, 1, hdr500) * 1e3

naxis1_250 = sxpar(hdr250, "naxis1")
naxis2_250 = sxpar(hdr250, "naxis2")
naxis1_350 = sxpar(hdr350, "naxis1")
naxis2_350 = sxpar(hdr350, "naxis2")
naxis1_500 = sxpar(hdr500, "naxis1")
naxis2_500 = sxpar(hdr500, "naxis2")

;pixsize250 = 3600. * sxpar(hdr250,"cdelt2") ;pixel size in arcsec
;pixsize350 = 3600. * sxpar(hdr350,"cdelt2")
;pixsize500 = 3600. * sxpar(hdr500,"cdelt2")
pixsize250 = 3600. * abs( sxpar(hdr250, "CD1_1") );pixel size in arcsec
pixsize350 = 3600. * abs( sxpar(hdr350, "CD1_1") )
pixsize500 = 3600. * abs( sxpar(hdr500, "CD1_1") )

if(keyword_set(astromcorr)) then begin ;apply astrometry correction
   oldcrval1 = sxpar(hdr250, 'crval1')
   fxaddpar, hdr250, 'crval1', (oldcrval1+astromcorr[0]/3600.)
   oldcrval2 = sxpar(hdr250, 'crval2')
   fxaddpar, hdr250, 'crval2', (oldcrval2+astromcorr[1]/3600.)
   oldcrval1 = sxpar(hdr350, 'crval1')
   fxaddpar, hdr350, 'crval1', (oldcrval1+astromcorr[0]/3600.)
   oldcrval2 = sxpar(hdr350, 'crval2')
   fxaddpar, hdr350, 'crval2', (oldcrval2+astromcorr[1]/3600.)
   oldcrval1 = sxpar(hdr500, 'crval1')
   fxaddpar, hdr500, 'crval1', (oldcrval1+astromcorr[0]/3600.)
   oldcrval2 = sxpar(hdr500, 'crval2')
   fxaddpar, hdr500, 'crval2', (oldcrval2+astromcorr[1]/3600.)
endif

make_2d, findgen(naxis1_250), findgen(naxis2_250), x_pix250, y_pix250

pixsize = [pixsize250, pixsize350, pixsize500]

;nim250 = mrdfits(pswfits, 6, nhdr250)*1e3 ;read in the uncertainty map
;nim350 = mrdfits(pmwfits, 6, nhdr350)*1e3
;nim500 = mrdfits(plwfits, 6, nhdr500)*1e3
nim250 = mrdfits(pswfits, 2, nhdr250) * 1e3 ;read in the uncertainty map
nim350 = mrdfits(pmwfits, 2, nhdr350) * 1e3
nim500 = mrdfits(plwfits, 2, nhdr500) * 1e3

adxy, hdr250, sra,sdec, sx250, sy250

; filter out sources not in the map, just check 250
sgood=where((sx250 gt 0) and (sx250 lt naxis1_250) and (sy250 gt 0) and (sy250 lt naxis2_250) and finite(im250[sx250,sy250]),n_src)
if(n_src eq 0) then stop

sx250 = sx250[sgood]
sy250 = sy250[sgood]
sra = sra[sgood]
sdec = sdec[sgood]
f_src = f_src[sgood]

adxy,hdr350,sra,sdec,sx350,sy350
adxy,hdr500,sra,sdec,sx500,sy500

bad = where(finite(im250) eq 0 or finite(nim250) eq 0 or nim250 eq 0, ncomp=n_goodpix250)
if(bad[0] ne -1) then begin
   im250[bad]=0.
   nim250[bad]=1.
endif

bad = where(finite(im350) eq 0 or finite(nim350) eq 0 or nim250 eq 0)
if(bad[0] ne -1) then begin
   im350[bad]=0.
   nim350[bad]=1.
endif

bad = where(finite(im500) eq 0 or finite(nim500) eq 0 or nim250 eq 0)
if(bad[0] ne -1) then begin
   im500[bad]=0.
   nim500[bad]=1.
endif

; work out pseudo p-stat to weight models by calculate prob of source with f>fs in one beam
; work out input list area crudely by just assuming contiguous area from min to max ra & dec of sources
p_src=fltarr(n_src)
;sarea=((max(sra)-min(sra))*cos(mean(sdec)*!pi/180.))*(max(sdec)-min(sdec)) ;change this to the number of pixels times the area of each pixel
sarea = n_goodpix250 * (pixsize250 / 3600.)^2. ;number of pixels times the area of each pixel

for k = 0L, n_src-1 do begin
   tmp = where(f_src ge f_src[k], ntmp)
   ns = ntmp/sarea
   p_src[k] = !pi*(prfsize[0]/3600.)^2/4./alog(2.)*ns
endfor
;---------------
; load PRF from file
;---------------
pfwhm = prfsize / pixsize
segsrc = lonarr(n_src) + 2     ;why?

;segment the map and decide which source is in which segment, change
;segthresh upwards to make segment smaller, default is 0.5
if not keyword_set(noseg) then lstdrv_segsrc,sx250,sy250,pfwhm[0],pfwhm[1],pfwhm[2],hdr250,im250,nim250,hdr350,im350,nim350,hdr500,im500,nim500,segsrc,bsn=segthresh

msegsize = max(histogram(segsrc))
if (msegsize gt 3000) then begin
   print, 'Segment size is too large. Stopping ...'
   stop 
endif

fstruct={xid:0L,inra:0.,indec:0.,f250:0.,e250:0.,et250:0.,chi250:0.,f350:0.,e350:0.,et350:0.,chi350:0.,f500:0.,e500:0.,et500:0.,chi500:0.,gid:0L,gsize:0,bkg250:0.,bkg350:0.,bkg500:0.}
fcat=replicate(fstruct,n_src)

fcat.inra = sra
fcat.indec = sdec
fcat.xid = long(sgood)

;;
lstdrv_solvefluxes2,segsrc,sx250,sy250,sx350,sy350,sx500,sy500,pfwhm,im250,im350,im500,nim250,nim350,nim500,naxis1_250,naxis2_250,naxis1_350,naxis2_350,naxis1_500,naxis2_500,p_src,fcat=fcat,reconmap=reconmap,f_src=f_src

paxis = [11,11]
prf250 = psf_gaussian(fwhm=pfwhm[0], npix=paxis)  ;max(prf) = 1
prf350 = psf_gaussian(fwhm=pfwhm[1], npix=paxis)
prf500 = psf_gaussian(fwhm=pfwhm[2], npix=paxis)

bmask250 = fltarr(naxis1_250, naxis2_250)
bmask350 = fltarr(naxis1_350, naxis2_350)
bmask500 = fltarr(naxis1_500,naxis2_500)

make_2d, indgen(11)-5,indgen(11)-5, tx, ty

for k=0L, n_src-1 do begin   ;put down a prf at every source location
   bmask250[sx250[k]+tx, sy250[k]+ty] += prf250 
   bmask350[sx350[k]+tx, sy350[k]+ty] += prf350
   bmask500[sx500[k]+tx, sy500[k]+ty] += prf500
endfor

region250 = where(bmask250 ne 0.)
region350 = where(bmask350 ne 0.)
region500 = where(bmask500 ne 0.)

;Result = LINFIT(X, Y). The result is a two-element vector containing the linear model parameters [A, B].
test250 = linfit(bmask250[region250], (im250 - reconmap.rmap250)[region250])
test350 = linfit(bmask350[region350], (im350 - reconmap.rmap350)[region350])
test500 = linfit(bmask500[region500], (im500 - reconmap.rmap500)[region500])

sim250 = im250 - test250[0] ;original image - background at each band
sim350 = im350 - test350[0]
sim500 = im500 - test500[0]

res250 = sim250 - reconmap.rmap250
res350 = sim350 - reconmap.rmap350
res500 = sim500 - reconmap.rmap500

mwrfits, blank, name + '.fits', bhdr250, /create ;output residual at 250, 350 and 500 microns
mwrfits, res250, name + '.fits', hdr250
mwrfits, res350, name + '.fits', hdr350
mwrfits, res500, name + '.fits', hdr500
;;=========================================================================
;;beginning of the second pass
;;=========================================================================
; find missing sources from the residual maps (=image - reconstructed map of known sources and the background)
;Find positive brightness perturbations (i.e stars) in an image 
find, res250 / nim250, nx250, ny250, nf250, sharp, round,3., pfwhm[0], [0.,3.], [-3.,3.],/silent
find, res350 / nim350, nx350, ny350, nf350, sharp, round,3., pfwhm[1], [0.,3.], [-3.,3.],/silent
find, res500 / nim500, nx500, ny500, nf500, sharp, round,3., pfwhm[2], [0.,3.], [-3.,3.],/silent

xyad, hdr250, nx250, ny250, nra250, ndec250; keep 250 list and filter out sources in 350 and 500
xyad, hdr350, nx350, ny350, nra350, ndec350
xyad, hdr500, nx500, ny500, nra500, ndec500

nra = nra250
ndec = ndec250

for i=0L,n_elements(nx350)-1 do begin
   sep=sqrt(((nra-nra350[i])*cos(ndec*!pi/180.))^2+(ndec-ndec350[i])^2)
   if(min(sep) gt 20./3600.) then begin
      nra=[nra, nra350[i]]
      ndec=[ndec,ndec350[i]]
   endif
endfor

for i=0L,n_elements(nx500)-1 do begin
   sep=sqrt(((nra-nra500[i])*cos(ndec*!pi/180.))^2+(ndec-ndec500[i])^2)
   if(min(sep) gt 20./3600.) then begin
      nra=[nra, nra500[i]]
      ndec=[ndec,ndec500[i]]
   endif
endfor

; check if any are within 0.5 FHWM of known sources, just check 250?
nm=-1
for i=0L, n_elements(nra)-1 do begin
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

;nsources={nra:nra,ndec:ndec,nx250:nx250,ny250:ny250,nx350:nx350,ny350:ny350,nx500:nx500,ny500:ny500,obs250:fltarr(nnew),e250:fltarr(nnew),obs350:fltarr(nnew),e350:fltarr(nnew),obs500:fltarr(nnew),e500:fltarr(nnew)}

;combine sources from the input list and new sources found in the residual map
sx250 = [sx250,nx250]
sy250 = [sy250,ny250]
sx350 = [sx350,nx350]
sy350 = [sy350,ny350]
sx500 = [sx500,nx500]
sy500 = [sy500,ny500]
p_src = [p_src,fltarr(nnew)+1.]
f_src = [f_src,fltarr(nnew)+min(f_src)]
segsrc = lonarr(n_src+nnew)+2

if not keyword_set(noseg) then lstdrv_segsrc,sx250,sy250,pfwhm[0],pfwhm[1],pfwhm[2],hdr250,im250,nim250,hdr350,im350,nim350,hdr500,im500,nim500,segsrc,bsn=segthresh

fcatmod = replicate(fstruct,n_src+nnew)
fcatmod.inra = [sra, nra]
fcatmod.indec = [sdec,ndec]
fcatmod.bkg250 = fltarr(n_src+nnew) + test250[0]  ;the background value at 250, 350 and 500
fcatmod.bkg350 = fltarr(n_src+nnew) + test350[0]
fcatmod.bkg500 = fltarr(n_src+nnew) + test500[0]
nids = fltarr(nnew) - 99. ;any new sources found in the residual map have xid = -99
fcatmod.xid = [fcat.xid,nids]

;sim250, sim350 and sim500 are images with background subtracted off
lstdrv_solvefluxes2,segsrc,sx250,sy250,sx350,sy350,sx500,sy500,pfwhm,sim250,sim350,sim500,nim250,nim350,nim500,naxis1_250,naxis2_250,naxis1_350,naxis2_350,naxis1_500,naxis2_500,p_src,fcat=fcatmod,reconmap=reconmap,f_src=f_src

; work out true error, 250
paxis = [11,11]

bmask250 = fltarr(naxis1_250,naxis2_250)
make_2d,indgen(11)-5,indgen(11)-5,tx,ty

for k=0L, n_src-1 do bmask250[sx250[k]+tx, sy250[k]+ty] = 1.

region=where(bmask250 eq 1.)
tmap250 = im250 - reconmap.rmap250
prf250 = psf_gaussian(fwhm=pfwhm[0],npix=paxis)
smap250 = convol(tmap250,prf250) / total(prf250^2)
nm250 = convol(1./nim250^2,prf250^2)
err250 = stddev(smap250[region])^2 - 1./nm250[sx250,sy250]
good = where(err250 gt 0.)
fcatmod[good].et250 = sqrt(fcatmod[good].e250^2 + err250[good])

;350
bmask350=fltarr(naxis1_350,naxis2_350)
make_2d,indgen(11)-5,indgen(11)-5,tx,ty
for k=0L,n_src-1 do bmask350[sx350[k]+tx, sy350[k]+ty]=1.

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

for k=0L, n_src-1 do bmask500[sx500[k]+tx, sy500[k]+ty]=1.
region=where(bmask500 eq 1.)
tmap500=im500-reconmap.rmap500
prf500=psf_gaussian(fwhm=pfwhm[2],npix=paxis)
smap500=convol(tmap500,prf500)/total(prf500^2)
nm500=convol(1./nim500^2,prf500^2)
err500=stddev(smap500[region])^2-1./nm500[sx500,sy500]
good=where(err500 gt 0) 
fcatmod[good].et500=sqrt(fcatmod[good].e500^2+err500[good])

end
