# setting up the environment variables
source HermesXIDSetup.csh

# Running the script
set imfolder = "/Local/Users/lw94/hermes/maps/version16/"
set catfolder = '/Local/Users/lw94/hermes/cats/version16_SF/'
set outfolder = '/Local/Users/lw94/hermes/cats/SF_xid/version16/'
echo $imfolder

# foreach field ("Abell_1689" "Abell_1835" "Abell_2390")

set field = "Abell_1689"
echo "lingyuxxxxxx$field"
    
set stem =  '-2010-11-16_Naive'
set prefix = 'StarFinder_Cat_'
    
echo  $imfolder$field$stem'_PSW.fits.gz'  $imfolder$field$stem'_PMW.fits.gz'  $imfolder$field$stem'_PLW.fits.gz' \
      $catfolder$prefix$field$stem'_PSW_SXT.fits'   'test.fits' 'test'
	
/Applications/itt/idl/idl80/bin/idl SCAT_XID_batch -args \
	$imfolder$field$stem'_PSW.fits.gz'  $imfolder$field$stem'_PMW.fits.gz'  $imfolder$field$stem'_PLW.fits.gz' \
	$catfolder$prefix$field$stem'_PSW_SXT.fits'   $outfolder$field$stem'_PSW_SCAT_XID.fits'   $outfolder$field$stem'_RESID_MAPS' \
	$outfolder$field$stem'_PSW_SCAT_XID_EXTRAS.fits' 
    
# end
