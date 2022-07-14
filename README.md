# RAVEN_parameters
development space for hyspex spectral parameters

Repo containing code and example outputs for "CRISM-like" spectral parameters calculated on HySpex data and point spectra collected with an oreXpress spectrometer. Comparison for point spectra can be made to USGS spectral library spectra contained in /Output/Library/. Other /Output/ folders are labled by data (YYMMDD). /Output/220713/ contains example spectra collected with the oreXpress from the UTK rock garden collection (mostly metamorphic rocks with alteration minerals, including chlorite and some Al-clays.)

Calculation of spectral parameters was modelled as close as possible to the calcuations done in the CRISM Analysis Toolkit add-on to ENVI. In most cases, the calculation is verbatim, but in some cases the HySpex instrument does not go out to far enough wavelengths (past 2500 nm), so spectral parameters could not be copied exactly. This mostly affects the carbonate browse product, for which BD3000 was replaced with BDCARB. 
