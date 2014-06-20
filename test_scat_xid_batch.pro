# setting up the environment variables

source HermesXIDSetup.csh

# running the script

set imfolder = '/export/margaux/hermes/scat/maps/version16/'
set catfolder = '/export/margaux/hermes/scat/cats/version16/'
set outfolder = '/export/margaux/hermes/scat/scat_xid/version16/'

set stem =  'ECDFS-2010-11-16_Naive'


echo  $imfolder$stem'_PSW.fits.gz'  $imfolder$stem'_PMW.fits.gz'  $imfolder$stem'_PLW.fits.gz' \
    $catfolder$stem'_PSW_SXT.fits'   'test.fits' 'test'


idl SCAT_XID_batch -args \
    $imfolder$stem'_PSW.fits.gz'  $imfolder$stem'_PMW.fits.gz'  $imfolder$stem'_PLW.fits.gz' \
    $catfolder$stem'_PSW_SXT.fits'   $outfolder$stem'_PSW_SCAT_XID.fits'   $outfolder$stem'_RESID_MAPS'

