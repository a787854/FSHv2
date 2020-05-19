[![Language](https://img.shields.io/badge/python-2.7/3.6%2B-blue.svg)](https://www.python.org/)
![Latest version](https://img.shields.io/badge/latest%20version-v1.2-yellowgreen.svg)
[![License](https://img.shields.io/badge/License-GNU%203.0-blue.svg)](https://github.com/leiyangleon/FSH/blob/master/LICENSE)
[![Citation](https://img.shields.io/badge/DOI-10.1109/TGRS.2018.2860590-blue)](https://doi.org/10.1109/TGRS.2018.2860590)

<img src="preview.jpg" width="100%">

# Forest Stand Height (FSH) Python Scripts

- Windows (Anaconda Prompt)
- Linux
- OSX

This software performs the automated forest height inversion and mosaicking from spaceborne repeat-pass L-band HV-pol InSAR correlation magnitude data (e.g. JAXA’s ALOS-1/2, and the future NASA-ISRO’s NISAR) that have been pre-processed by JPL’s ROI_PAC and/or ISCE programs.

Produced by the University of Massachusetts Microwave Remote Sensing Laboratory. 

Yang Lei (ylei@caltech.edu, leiyangfrancis@gmail.com), Paul Siqueira (siqueira@umass.edu).

***Besides contacting the authors via email for technical questions, it is also recommended to [refer to the closed issues or open a new issue on GitHub](https://github.com/leiyangleon/FSH/issues).***


===================================================================================

## Documentation

### 1. [Installation](./docs/installation.md)

### 2. [How to use](./docs/how_to_run.md)

***Please note that users will have to be able to generate their own interferograms using JPL's ROIPAC or ISCE (insarApp.py workflow up to v2.2; stripmapApp.py workflow for v2.2+) software in order to use the FSH software. Information and tutorials regarding ISCE and interferogram processing can be found at [here](http://earthdef.caltech.edu/projects/isce_forum/boards). To those interested in pursuing their own interferogram processing, we used modified ISCE [scripts](./ISCE_processing_scripts) which are provided with the FSH software, and described [here](./docs/isce_preprocessing.md).*** 


***Here, we provide three test examples for FSH using JAXA's ALOS PALSAR interferograms that have been generated by ROIPAC and ISCE (with both insarApp.py and stripmapApp.py workflows in v2.2), respectively.***

- [test_example_ISCE_insarApp](./docs/test_example_ISCE_insarApp.md)
- [test_example_ISCE_stripmapApp](./docs/test_example_ISCE_stripmapApp.md)
- [test_example_ROIPAC](./docs/test_example_ROIPAC.md)

### 3. [List of Python files and corresponding descriptions](./docs/reference.md)

### 4. [Update notes](./docs/update.md)




===================================================================================



  
## [License](./docs/license_note.md)


===================================================================================

## Citation
Lei, Y., Siqueira, P., Torbick, N., Ducey, M., Chowdhury, D. and Salas, W., 2019. [Generation of Large-Scale Moderate-Resolution Forest Height Mosaic With Spaceborne Repeat-Pass SAR Interferometry and Lidar](https://ieeexplore.ieee.org/document/8439086). IEEE Transactions on Geoscience and Remote Sensing, 57(2), pp.770-787.

===================================================================================

## [Related publications and book chapters](./docs/reference_list.md)

===================================================================================

## Contributors
- Yang Lei
- Tracy Whelen
- Paul Siqueira
- Diya Chowdhury
- Gerard Ruiz Carregal
- Thannarot Kunlamai
- Tim Mayer
- Helen Baldwin
- Simon Kraatz
