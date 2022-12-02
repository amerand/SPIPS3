# Spectro Photo Interferometry of Pulsating stars

## What is this?

This is a python3.10 (replacing the old 2.7 version) implementation of a parallax of pulsation method for Cepheids stars, described in [Mérand et al. (Astronomy & Astrophysics 584-80, 2015)](http://adsabs.harvard.edu/abs/2015A%26A...584A..80M). As of November 2022, [17 refereed articles](https://ui.adsabs.harvard.edu/search/filter_property_fq_property=AND&filter_property_fq_property=property%3A%22refereed%22&fq=%7B!type%3Daqp%20v%3D%24fq_property%7D&fq_property=(property%3A%22refereed%22)&q=%20full%3A%22SPIPS%22%20year%3A2015-2050%20author%3A%22M%C3%A9rand%22&sort=date%20desc%2C%20bibcode%20desc&p_=0) have been published with SPIPS.


## Quick Start / Example

Quick Start:
 - download all files
 - in python3: `import delta_cep.py`
   - this will load spips.py
   - first time you load spips.py, lots of models will be downloaded: [ATLAS9 models](http://wwwuser.oats.inaf.it/castelli/grids.html) and [SATLAS models](http://cdsarc.u-strasbg.fr/viz-bin/Cat?J/A%2bA/554/A98).
 - run `delta_cep.show(delta_cep.p_fourier)` to show the model with Fourier parameters
 - run `delta_cep.show(delta_cep.p_splines)` to show the model with Splines parameters
 - run `delta_cep.fit(delta_cep.p_splines)` to run a fit. Check the inside of the function to see how it works.

The result of the model is shown below:
![Fig1](delta_cep.png)
- the upper left panel (a) shows the phased radial velocity data (points) and model (line)
- the middle left panel (b) shows the phased effective temperature data (points) and model (line)
- the lower left panel (c) shows the phased interferometric angular diameter data (points) and model (lines). Note that the different colors show the impact of the effects of the interferometric baseline on the diameter measurements, due to the presence of an circum-stellar envelop
- the panels on the right hand side (d to v) show the photometric data.

`delta_cep.fitsDemo(mode='export')` will export the model and the data to a FITS file (`delta_cep.fits`), which can be read using `delta_cep.fitsDemo(mode='import')`. Note that `delta_cep.fitsDemo(mode='import')` recomputes the model, rather than plotting the one in the FITS files. One can easily write their own routines to read and display the FITS file created by SPIPS, as they are self-explanatory. 

## Dependencies and known issues

- `numpy`, `scipy`, `matplotlib` and `astropy.io.fits` 

## License (BSD)

Copyright (c) 2017-2022, Antoine Mérand
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
