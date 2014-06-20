pro testXIDcode
; simple test script to extract SPIRE fluxes at positions of 250
; micron detected sources from SCAT in version 15 ECDFS
; assumes you are somewhere on pact network
; also assumes Gaussian PRF, and map in Jy/beam

imfolder='/export/margaux/hermes/scat/maps/version15/'
catfolder='/export/margaux/hermes/scat/cats/version15/'

scat=mrdfits(catfolder+'ECDFS-2010-06-15_Naive_PSW_SXT.fits',1)
good=where(finite(scat.flux)) ; for some reason some SCAT sources have NaN fluxes
scat=scat[good]
sra=scat.ra
sdec=scat.dec
f_src=scat.flux

pswfits=imfolder+'ECDFS-2010-06-15_Naive_PSW.fits'
pmwfits=imfolder+'ECDFS-2010-06-15_Naive_PMW.fits'
plwfits=imfolder+'ECDFS-2010-06-15_Naive_PLW.fits'

prfsize=[18.2,25.2,36.3]


lstdrv_readSCAT,pswfits,pmwfits,plwfits,sra,sdec,prfsize,fcat,f_src=f_src,astromcorr=astromcorr,fcat2=fcat2,name='ecdfs_scat'

stop
end
@lstdrv_readSCAT
@lstdrv_solvefluxes2
@lstdrv_segsrc
@lstdrv_modselLASSO
@lstdrv_initsolveSP
