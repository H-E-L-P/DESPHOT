; pro SCAT_XID_batch
; simple test script to extract SPIRE fluxes at positions of 250
; micron detected sources from SCAT in version 15 ECDFS
; assumes you are somewhere on pact network
; also assumes Gaussian PRF, and map in Jy/beam

;.compile lstdrv_readscat
;.compile lstdrv_solvefluxes2
;.compile lstdrv_segsrc
;.compile lstdrv_modselLASSO
;.compile lstdrv_initsolveSP
;-----------------------------------

imfolder = '/archi/scratch/lw94/smap/'
catfolder = '/archi/scratch/lw94/smap/StarFinder/'
outfolder = '/archi/scratch/lw94/smap/SF_xid/'

all_fields = [ 'Abell-1689', $     ;0
               'Abell-1835', $     ;1
               'Abell-2218', $     ;2 
               'Abell-2219', $     ;3 
               'Abell-2390', $     ;4 
               'Abell-370', $      ;5
               'Cl0024+16', $      ;6
               'GOODS-N', $        ;7
               'GOODS-S', $        ;8
               'MS0451.6-0305', $  ;9 
               'MS1054.4-0321', $  ;10 
               'MS1358+62', $      ;11
               'RXJ0152.7-1357', $ ;12
               'RXJ13475-1145', $  ;13
               'ECDFS', $          ;14
               'COSMOS', $         ;15
               'ADFS', $           ;16
               'BOOTES-HerMES', $  ;17
               'CDFS-SWIRE', $     ;18
               'EGS-HerMES', $     ;19
               'ELAIS-N1-HerMES', $;20
               'ELAIS-S1-SWIRE', $ ;21
               'FLS', $            ;22
               'Groth-Strip', $    ;23
               'Lockman-East-ROSAT', $;24 
               'Lockman-North', $  ;25
               'Lockman-SWIRE', $  ;26
               'UDS', $            ;27
               'VVDS', $           ;28
               'XMM-LSS-SWIRE']    ;29

                                ;stem =  '-2010-11-16_Naive'
                                ;prefix = 'backbox150_SF_Cat_'

prefix = '_image_'
suffix = '_SMAP_v3.0.fits'
cat_prefix = 'mincorr_0.700000_backbox150_SF_Cat_'

;for ff = 0, 0 do begin

ff = 20

field = all_fields[ff]

argument=[imfolder+field+prefix+'250'+suffix, imfolder+field+prefix+'350'+suffix, imfolder+field+prefix+'500'+suffix, $
          catfolder + cat_prefix + field + prefix + '250' + suffix, outfolder + field + '_SMAP_XID.fits', $
          outfolder + field + '_RESID_MAPS', $
          outfolder + field + '_SMAP_XID_mod.fits']

;; argument=[imfolder+field+prefix+'250'+suffix, imfolder+field+prefix+'350'+suffix, imfolder+field+'500'+suffix, $
;;           catfolder + prefix + field + stem + '_PSW.fits.gz', outfolder+field+stem+'_SCAT_XID_new.fits', $
;;           outfolder+field+stem+'_RESID_MAPS_new', $
;;           outfolder + field + stem + '_SCAT_XID_mod.fits']
;-----------------------------------
;argument = COMMAND_LINE_ARGS( )

message, 'Command line arg 0'+argument[0], /inf
message, 'Command line arg 1'+argument[1], /inf
message, 'Command line arg 2'+argument[2], /inf

pswfits = argument[0]
pmwfits = argument[1]
plwfits = argument[2]
scat_cat_file = argument[3]
out_cat_file = argument[4]
out_resid_maps_stem = argument[5]
out_cat_file2 = argument[6]

scat0 = mrdfits(scat_cat_file,1, hd)
good = where(finite(scat0.flux)) ; for some reason some SCAT sources have NaN fluxes
scat = scat0[good]
sra = scat.ra
sdec = scat.dec
f_src = scat.flux

                                ;prfsize = [18.2, 25.2, 36.3]
prfsize = [18.15, 25.15, 36.3]

lstdrv_readSCAT, pswfits, pmwfits, plwfits, sra, sdec, prfsize, fcat, f_src=f_src,astromcorr=astromcorr,name=out_resid_maps_stem,fcatmod=fcatmod

fcat_final = replicate(fcatmod[0], n_elements(scat0))
fcat_final.xid = -1
fcat_final[round(fcatmod.xid)] = fcatmod

;; writing out the file with fluxes for original source list
mwrfits, fcat_final, out_cat_file , /create

;; writing out the file with fluxes for deeper source list
mwrfits, fcatmod, out_cat_file2 , /create
;endfor

end

