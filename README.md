
DESPHOT, described in the papers above, has been used to find the flux for objects in a prior catalogue. In the first ([Roseboom et al. 2010](http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:1009.1658)) and second iteration of the algorithm ([Roseboom et al. 2012](http://adsabs.harvard.edu/abs/2012MNRAS.419.2758R))(XID), 24 micron and radio sources were used as prior positions. For the final iteration (DESPHOT) [Wang et al. 2014](http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:1312.0552) used the 250 micron sources detected with Starfinder as prior positions.


## Installation

DESPHOT requires IDL, including the IDLlib and astron library. Make sure your file cotaining the DESPHOT code is also in your IDLPATH

## Usage
DESPHOT requires the three SPIRE maps (ideally from SMAP) and a prior catalogue. Other parameters include the prf and segmentation threshold.

An example test script for running DESPHOT on one field can be found in testXIDcode.pro


---
***Shell scripts***
----
-   HermesAstroSetup.csh: Old shell script to setup IDL environment
-   HermesXIDSetup.csh: Lingyu's specific shel script for her version of XIDcode
-   scat_xid_batch.csh
-   test_scat_xid_batch.csh

---
***IDL code***
----

-   lstdrv_readmaps.pro: 
-   lstdrv_readscat.pro: uses a two-pass strategy to do list-driven sources extraction, also finds new sources in the residual maps.

    -   lstdrv_segsrc.pro: Segment the map
    -   lstdrv_solvefluxes2.pro: Solve the fluxes
        -   lstdrv_initsolvesp.pro: Creates pointing matrices 
        -   lstdrv_modsellasso.pro: Solves for flux using LASSO

-   lstdrv_segsrcsmap.pro: ???
-   mkspiremap.pro: ???
-   lstdrv_readsmapct.pro: ???
-   lstdrv_readsmap.pro: ???

-   test_arg.pro: simple program to test args??
-   test_scat_xid_batch.pro: Test batch script
    -   SCAT_XID_batch.pro: simple test script to extract SPIRE fluxes at positions of 250 micron detected sources from SCAT for a batch
-   testXIDcode.pro: Simple test script to extract SPIRE fluxes at positions of 250 micron detected sources from SCAT catalogue from one field
-   testXIDcodev16.pro: Simple test script to extract SPIRE fluxes at positions of 250 micron detected sources from SCAT catalogue from one field (v16)


## Contributing

1. Fork it!
2. Create your feature branch: `git checkout -b my-new-feature`
3. Commit your changes: `git commit -am 'Add some feature'`
4. Push to the branch: `git push origin my-new-feature`
5. Submit a pull request :D

## History

TODO: Write history

## Credits
[Roseboom et al. 2010](http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:1009.1658)

[Roseboom et al. 2012](http://adsabs.harvard.edu/abs/2012MNRAS.419.2758R)

[Wang et al. 2014](http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:1312.0552)


## License

TODO: Write license