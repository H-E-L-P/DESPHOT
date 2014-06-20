pro lstdrv_modselLASSO,apmat,bvec,psrc,active,xvec,snsrc,f_src=f_src,thresh=thresh

; adjust pmat for weights psrc

pmat=dblarr(snsrc,snsrc)
for i=0,snsrc-1 do begin
  pmat[i,*]=apmat[i,*]/psrc[i]
endfor



spmat=sprsin(transpose(pmat))




; determine initial negative gradient
w=dblarr(snsrc)
xvec=dblarr(snsrc)
res=bvec-matrix_multiply(pmat,xvec,/atranspose)
pres=res
w=-matrix_multiply(pmat,res)

dead=indgen(snsrc)
nact=0
lambda=-min(w,dind)
JB=total(abs(xvec))
if(not keyword_set(thresh)) then thresh=total(f_src*psrc)*1.1
status=intarr(snsrc)
status[dind]=1
active=where(status eq 1,nact,comp=inactive)

xxt=matrix_multiply(pmat,pmat,/btranspose)
t1=systime(/sec)
niter=0L
id2=-1
il=0
mid=-1
while(JB lt thresh and lambda gt 1e-30 and niter lt 10.*snsrc  and il lt 10) do begin
ttmp2=systime(/sec)
gamma=fltarr(snsrc)
tmp=active
make_2d,tmp,tmp,ina,inb
xxta=double(xxt[ina,inb])
if nact gt 1 then begin 
    ones=(fltarr(nact)+1.)
    ttmp=systime(/sec)


DEFSYSV, '!XLIB', EXISTS = usexlib
if(usexlib eq 1) then begin
; IF USING MP MACHINE AND HAVE libXID.so USE THESE LINES  
;--------------  
    ntmp=long(nact)
    blah=call_external('libXID.so','idlinvert_',$
             ntmp,xxta,value=[1,0])
    gammaA=xxta#ones

;-------------  
endif else begin
; IDL ONLY LINE
;-------------
   gammaA=la_invert(xxta)#ones
;------------
endelse






  ix=(systime(/sec)-ttmp)
endif else begin
    gammaA=1./xxta
    ix=0.
endelse
gamma[active]=gammaA


d1=-1
if(nact ne snsrc) then d1=(lambda+w[inactive])/(1.-(matrix_multiply(xxt[inactive,*],gamma)))
d2=-xvec[active]/gamma[active]

d=lambda
index=0
id1=-1

td1=where(d1 gt 1e-20 and d1 le d)
if(td1[0] ne -1) then begin
 d=min(d1[td1],tmp)
 id1=inactive[td1[tmp]]
 index=1
 
 ;endif else begin 
 ;    if(n_elements(td1) gt 1) then begin
     ; find next best 
 ;    sd=sort(d1[td1])
 ;    id1=inactive[td1[sd[1]]]
 ;    index=1
 ;endif
;endelse
 
endif
id2=-1
td2=where(d2 gt 1e-20 and d2 lt d)
if(td2[0] ne -1) then begin
 d=min(d2[td2],tmp)
 id2=active[td2[tmp]]
 index=2
endif

; need a catch for cases where the gradient is already negative at 0
vbad=where(d2 eq 0 and gamma[active] lt 0) 
if(vbad[0] ne -1)then begin
d=0
id2=active[vbad[0]]
status[active[vbad[0]]]=0
index=2
endif 


xvec=xvec+d*gamma
wtf=where(xvec lt 0)
omid=mid
if(index eq 1) then begin
    status[id1]=1
    mid=id1
endif
if(index eq 2) then begin
    status[id2]=0
    mid=id2
endif
;check for loops
 if(mid eq omid) then begin
     il=il+1
 endif else begin
     il=1
 endelse

active=where(status eq 1,nact,comp=inactive)
if(nact ne snsrc) then xvec[inactive]=0.
bx=fltarr(snsrc)
ttmp=systime(/sec)
tpmat=pmat[active,*]
bx=matrix_multiply(tpmat,xvec[active],/atranspose)

res=bvec-bx

w=-sprsax(spmat,res,/double)
rw=systime(/sec)-ttmp
JB=total(abs(xvec))


lambda=lambda-d

    
   
  niter=niter+1
etmp=systime(/sec)-ttmp2
if (niter mod 500 eq 0) then print,'iter ',niter,il,ix,rw,(etmp-ix-rw),index,id1,id2,lambda,d
;stop
endwhile


xvec=xvec/psrc

print,'finished ',(systime(/sec)-t1)

end
