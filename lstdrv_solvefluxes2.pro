pro lstdrv_solvefluxes2,segsrc,sx250,sy250,sx350,sy350,sx500,sy500,pfwhm,im250,im350,im500,nim250,nim350,nim500,naxis1_250,naxis2_250,naxis1_350,naxis2_350,naxis1_500,naxis2_500,p_src,fcat=fcat,reconmap=rmap,modsel=modsel,f_src=f_src


paxis=[13,13]

snum=max(segsrc)
; set up arrays to contain reconstructed map 
rmap250=fltarr(naxis1_250,naxis2_250)
rmap350=fltarr(naxis1_350,naxis2_350)
rmap500=fltarr(naxis1_500,naxis2_500)
for seg=1L,snum do begin
tseg=where(segsrc eq seg,snsrc)

if(snsrc gt 0) then begin
fcat[tseg].gid=seg
fcat[tseg].gsize=snsrc
; array to contain 'detection' flag, i.e. did we select this source in
; the best fit model or not
aflag250=intarr(snsrc)
aflag350=intarr(snsrc)
aflag500=intarr(snsrc)

print,'starting segment ',seg,' size: ',snsrc


bu=dblarr(snsrc)
bl=dblarr(snsrc)
; build matrices for this segment needed by solver
make_2d,indgen(naxis1_250),indgen(naxis2_250),x_pix,y_pix
lstdrv_initsolveSP,sx250,sy250,tseg,prf250,pfwhm[0],paxis,im250,nim250,x_pix,y_pix,seg,bl,bu,p_src[tseg],psig250,pmat250,amat250,sn_map250,bvec250,xvec250,iBIC250,sx_pix=sx_pix250,sy_pix=sy_pix250
make_2d,indgen(naxis1_350),indgen(naxis2_350),x_pix,y_pix
lstdrv_initsolveSP,sx350,sy350,tseg,prf350,pfwhm[1],paxis,im350,nim350,x_pix,y_pix,seg,bl,bu,p_src[tseg],psig350,pmat350,amat350,sn_map350,bvec350,xvec350,iBIC350,sx_pix=sx_pix350,sy_pix=sy_pix350
make_2d,indgen(naxis1_500),indgen(naxis2_500),x_pix,y_pix
lstdrv_initsolveSP,sx500,sy500,tseg,prf500,pfwhm[2],paxis,im500,nim500,x_pix,y_pix,seg,bl,bu,p_src[tseg],psig500,pmat500,amat500,sn_map500,bvec500,xvec500,iBIC500,sx_pix=sx_pix500,sy_pix=sy_pix500

active=lonarr(snsrc)
if(snsrc gt 1) then begin

    pwts=(p_src[tseg])         
; maximum threshold for sum of weighted fluxes (total(pwts*sxvec)), if
; in doubt set this to something fairly big
;pwts=1./f_src[tseg]
thresh250=total(f_src[tseg]*pwts)*100.
; solve the linear equation for the source fluxes sxvec

lstdrv_modselLASSO,psig250,bvec250,pwts,active250,sxvec250,snsrc,f_src=f_src[tseg],thresh=thresh250

pwts=p_src[tseg]*exp(-sxvec250/10.)

thresh350=total(sxvec250*pwts)*100.

lstdrv_modselLASSO,psig350,bvec350,pwts,active350,sxvec350,snsrc,f_src=f_src[tseg],thresh=thresh350

pwts=p_src[tseg]*exp((-sxvec250/2.-sxvec350/2.)/10.)
thresh500=total(sxvec250*pwts)*100.
lstdrv_modselLASSO,psig500,bvec500,pwts,active500,sxvec500,snsrc,f_src=f_src[tseg],thresh=thresh500

;test=sxvec500[active500]-(bvec500/psig500[findgen(snsrc),findgen(snsrc)])[active500]

;bad=where(test gt 20)
;if(bad[0]) ne -1 then stop

endif else begin ; if there is only one source
    sxvec250=bvec250/psig250
    sxvec350=bvec350/psig350
    sxvec500=bvec500/psig500
    if(sxvec250 lt 0) then sxvec250=0.
    if(sxvec350 lt 0) then sxvec350=0.
    if(sxvec500 lt 0) then sxvec500=0.
    
endelse
; figure out which sources are active
active250=where(sxvec250 ne 0)

    active350=where(sxvec350 ne 0)
    active500=where(sxvec500 ne 0)
    
    if(active250[0] ne -1) then aflag250[active250]=1
    if(active350[0] ne -1) then aflag350[active350]=1
    if(active500[0] ne -1) then aflag500[active500]=1

nact=n_elements(active250)

; get the errors
e250=sqrt(1./psig250[indgen(snsrc),indgen(snsrc)])
e350=sqrt(1./psig350[indgen(snsrc),indgen(snsrc)])
e500=sqrt(1./psig500[indgen(snsrc),indgen(snsrc)])


make_2d,active250,active250,inax,inay
spsig250=(psig250[inax,inay])


if(n_elements(nactive250) gt 1) then begin

    spsig250=la_invert(spsig250,/double)
    e250[active250]=sqrt(spsig250[indgen(n_elements(active250)),indgen(n_elements(active250))])
; error can't be greater than the actual flux in the map
    t250=bvec250[indgen(n_elements(active250)),indgen(n_elements(active250))]/psig250[indgen(n_elements(active250)),indgen(n_elements(active250))]
    h250=where(t250 lt e250 and t250 gt sqrt(1./psig250[indgen(snsrc),indgen(snsrc)]))
    if(h250[0] ne -1) then e250[h250]=t250[h250]
endif else begin
    e250[active250]=sqrt(1./spsig250[active250,active250])
endelse
nact=n_elements(active350)
; get the errors
make_2d,active350,active350,inax,inay
spsig350=(psig350[inax,inay])
if(n_elements(active350) gt 1) then begin
spsig350=la_invert(spsig350,/double)
    e350[active350]=sqrt(spsig350[indgen(n_elements(active350)),indgen(n_elements(active350))])
    t350=bvec350[indgen(n_elements(active350)),indgen(n_elements(active350))]/psig350[indgen(n_elements(active350)),indgen(n_elements(active350))]
    h350=where(t350 lt e350 and t350 gt sqrt(1./psig350[indgen(snsrc),indgen(snsrc)]))
     if(h350[0] ne -1) then e350[h350]=t350[h350]
endif else begin
    e350[active350]=sqrt(1./spsig350[active350,active350])
endelse

nact=n_elements(active500)
; get the errors
make_2d,active500,active500,inax,inay
spsig500=(psig500[inax,inay])
if(n_elements(active500) gt 1) then begin
    spsig500=la_invert(spsig500,/double)
    e500[active500]=sqrt(spsig500[indgen(n_elements(active500)),indgen(n_elements(active500))])
    ; error can't be greater than the actual flux in the map
    t500=bvec500[indgen(n_elements(active500)),indgen(n_elements(active500))]/psig500[indgen(n_elements(active500)),indgen(n_elements(active500))]
    h500=where(t500 lt e500 and t500 gt sqrt(1./psig500[indgen(snsrc),indgen(snsrc)]))
     if(h500[0] ne -1) then e500[h500]=t500[h500]
endif else begin
    e500[active500]=sqrt(1./spsig500[active500,active500])
endelse




; reconstruct map from best fit solution

rmap250[sx_pix250,sy_pix250]=rmap250[sx_pix250,sy_pix250]+matrix_multiply(amat250,sxvec250)

rmap350[sx_pix350,sy_pix350]=rmap350[sx_pix350,sy_pix350]+matrix_multiply(amat350,sxvec350)

rmap500[sx_pix500,sy_pix500]=rmap500[sx_pix500,sy_pix500]+matrix_multiply(amat500,sxvec500)


; store flux values and errors for this segment

fcat[tseg].f250=sxvec250
fcat[tseg].e250=e250
fcat[tseg].f350=sxvec350
fcat[tseg].e350=e350
fcat[tseg].f500=sxvec500
fcat[tseg].e500=e500


; calculate chi2
chi250=fltarr(snsrc)
chi350=fltarr(snsrc)
chi500=fltarr(snsrc)

for s=0,snsrc-1 do begin
    dx=-sx250[tseg[s]]+sx_pix250
    dy=-sy250[tseg[s]]+sy_pix250
    good=where(abs(dx) le 5 and abs(dy) le 5,ngood)
    chi250[s]=total(((rmap250[sx_pix250[good],sy_pix250[good]]-sn_map250[good]))^2)/ngood
    dx=-sx350[tseg[s]]+sx_pix350
    dy=-sy350[tseg[s]]+sy_pix350
    good=where(abs(dx) le 5 and abs(dy) le 5,ngood)
    chi350[s]=total(((rmap350[sx_pix350[good],sy_pix350[good]]-sn_map350[good]))^2)/ngood
    dx=-sx500[tseg[s]]+sx_pix500
    dy=-sy500[tseg[s]]+sy_pix500
    good=where(abs(dx) le 5 and abs(dy) le 5,ngood)
    chi500[s]=total(((rmap500[sx_pix500[good],sy_pix500[good]]-sn_map500[good]))^2)/ngood
    
endfor
fcat[tseg].chi250=chi250
fcat[tseg].chi350=chi350
fcat[tseg].chi500=chi500

; put in upper limits for those that i killed
test=where(aflag250 eq 0)
if(test[0] ne -1) then begin
catin2=tseg[test]
tpmat250=(pmat250[test,test])

txvec250=[[bvec250[test]/tpmat250],[e250[catin2]]]
txvec250=max(txvec250,dim=2)

fcat[catin2].e250=txvec250

endif

test=where(aflag350 eq 0)
if(test[0] ne -1) then begin
catin2=tseg[test]
tpmat350=(pmat350[test,test])
txvec350=[[bvec350[test]/tpmat350],[e350[catin2]]]
txvec350=max(txvec350,dim=2)
fcat[catin2].e350=txvec350
endif


test=where(aflag500 eq 0)
if(test[0] ne -1) then begin
catin2=tseg[test]
tpmat500=(pmat500[test,test])
txvec500=[[bvec500[test]/tpmat500],[e500[catin2]]]
txvec500=max(txvec500,dim=2)
fcat[catin2].e500=txvec500

endif


bad=where(fcat[tseg].e500 lt 0.)
if(bad[0] ne -1) then stop
endif

DefSysv,'!xdev',exists=xd


endfor
rmap={rmap250:rmap250*(nim250),rmap350:rmap350*(nim350),rmap500:rmap500*(nim500)}

end
