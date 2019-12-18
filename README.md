# A *Swift* Survey of X-ray Clusters

This is all of the source code for the papers associated with the *Swift* Survey of X-ray Clusters.

## Software Requirements

### Python Requirements

This work relies on several third party Python modules which are *not* included with this repo.

* Astropy
* Astropy-regions
* Numpy
* Scipy
* Matplotlib
* tqdm
* emcee
* AplPy

In addition it requires the `load_catalogs` module that I wrote. It is available [here](https://github.com/boada/planckClusters/tree/master/catalogs).

### Other Software

To perform some of the specific x-ray analyses we need sets of standard x-ray astronomy tools.

* Ftools -- full install
* CAIO tools
* CALDB

A non-exhaustive list of x-ray specific tools required are:

* xrtpipeline
* xselect
* ximage
* xspec
* dmstat
* addascaspec
* grppha

## Installing

To get the third part Python modules

```
pip install -U astropy regions numpy scipy matplotlib tqdm emcee aplpy
```
should do it for you.

To setup non-python x-ray specific tooling, see the individual websites to get those tools up and running. The pipelines will often check to make sure the executables are available in your path, but they won't coach you on how to install that specific tool.

After installing all of the requirements, the notebooks should simply run.

```
git clone https://github.com/boada/swiftXRT.git
cd swiftXRT
```
And run the notebooks in order.


## Data Requirements

This project requires a bunch of extra data to work. Initially it requires the target list. See the source code in `load_catalogs` to understand what catalogs it initially loads.

Additionally, it requires a bunch of imaging data. It's not specifically required for the analysis, but if you would like to run notebooks 07 and 08 you will need that (or some) imaging. See the code in notebook `07` to get an idea of what imaging the pipeline expects.

It should be noted, getting all of the imaging is a non-trivial task.   

## Notebook overview

1. `01. Download XRT data` Downloads all of the XRT imaging for the target list loaded by `load_catalogs`
2. `02. Rereduce all SWIFT data` Reprocesses the downloaded data to make it more useful for analyses.
3. `03. Combine Reduced Products` Combines all of the individual observations into master events and exposure maps.
4. `03a. Remove Deep Low-z Fields` Remove fields which can be problematic in future processing
5. `04. Detect and Clean Sources` Source detection and false detection cleaning.
6. `05. Create and Fit Radial Profiles` Creates radial profiles and fits beta-models to each source
7. `06. Classify Sources` Classifies sources as point sources or extended sources
8. `07. Link PS1 and SDSS Imaging` Links the imaging into the data directories as a setup for the next notebook
9. `08. Make Images` Makes images, both XRT and Optical
10. `09. Extract Spectra` Extras spectra (and process through xspec) for each source

## Authors

* **Billie Thompson** - Rutgers University

## License

*Any use of the content of this project or repository for academic publication requires citation and acknowledgment.* If you are unsure how to properly cite and acknowledge this work, please contact me.

All of the code in this project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details. You are free to use and adapt it to your project, but please provide credit.

## Citation

I'll add the citation information here when it's available. Please contact me if you plan to use a portion of this work before the citation is available.

<!-- If you use any part of this work please reference our [paper](http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:1809.06378),
by using the following citation, produced by
[NASA ADS](http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:1810.12913):
```
@ARTICLE{2018arXiv180906378B,
   author = {{Boada}, S. and {Hughes}, J.~P. and {Menanteau}, F. and {Doze}, P. and
	{Barrientos}, L.~F. and {Infante}, L.},
    title = "{High Confidence Optical Confirmation of High Signal-to-Noise Planck Cluster Candidates}",
  journal = {ArXiv e-prints},
archivePrefix = "arXiv",
   eprint = {1809.06378},
 keywords = {Astrophysics - Cosmology and Nongalactic Astrophysics, Astrophysics - Astrophysics of Galaxies},
     year = 2018,
    month = sep,
   adsurl = {http://adsabs.harvard.edu/abs/2018arXiv180906378B},
  adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
``` -->
