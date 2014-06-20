pro testXIDcode
; simple test script to extract SPIRE fluxes at positions of 250
; micron detected sources from SCAT in version 15 ECDFS
; assumes you are somewhere on pact network
; also assumes Gaussian PRF, and map in Jy/beam

imfolder='/export/margaux/hermes/scat/maps/version16/'
catfolder='/export/margaux/hermes/scat/cats/version16/'

scat0=mrdfits(catfolder+'ECDFS-2010-11-16_Naive_PSW_SXT.fits',1, hd)
good=where(finite(scat0.flux)) ; for some reason some SCAT sources have NaN fluxes
scat=scat0[good]
sra=scat.ra
sdec=scat.dec
f_src=scat.flux

pswfits=imfolder+'ECDFS-2010-11-16_Naive_PSW.fits.gz'
pmwfits=imfolder+'ECDFS-2010-11-16_Naive_PMW.fits.gz'
plwfits=imfolder+'ECDFS-2010-11-16_Naive_PLW.fits.gz'

prfsize=[18.2,25.2,36.3]


lstdrv_readSCAT,pswfits,pmwfits,plwfits,sra,sdec,prfsize,fcat,f_src=f_src,astromcorr=astromcorr,fcat2=fcat2,name='ECDFS-2010-11-16_Naive_PSWXID_RESIDMAP'

fcat_final = replicate(fcat[0], n_elements(scat0))
fcat_final.xid = -1
fcat_final[round(fcat.xid)] = fcat

mwrfits, fcat_final, 'ECDFS-2010-11-16_Naive_PSWXID.fits', /create

stop
end
@lstdrv_readSCAT
@lstdrv_solvefluxes2
@lstdrv_segsrc
@lstdrv_modselLASSO
@lstdrv_initsolveSP
