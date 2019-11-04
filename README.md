# 2019 Tambien Group

Code and data associated with the Park et al. (2019) study on the Tambien Group of Ethiopia.

Park, Y., Swanson-Hysell, N.L., MacLennan, S.A., Maloof, A.C., Gebreslassie, M., Tremblay, M.M., Schoene, B., Alene, M., Anttila, E.S.C., Tesema, T., Haileab, B., 2019, The lead-up to the Sturtian Snowball Earth: Neoproterozoic chemostratigraphy time-calibrated by the Tambien Group of Ethiopia. *https://doi.org/10.1130/B35178.1*

## Environment

Some notebooks require the pyGPlates module (https://www.gplates.org/docs/pygplates/) that enables the functionality of the GPlates software package to be programmatically accessed using Python. With the exception of the pyGPlates module, which needs to be installed locally and added to the Python path, the computational environment for notebooks that require Python 3 is specified within `2019_Tambien_Group_python3.yml`. The computational environment for notebooks that require Python 2 is specified within `2019_Tambien_Group_python2.yml`.

## Repository Structure

* `Code/`
    * This folder contains all code used to perform the analyses described in this study.
    * `Composite_Chemostratigraphy/`
        * This folder contains code and data that was used to generate the Tonian-Cryogenian composite chemostratigraphy.
    * `Geochronology/`
        * This folder contains code and data that was used to generate figures associated with geochronologic data.
    * `Global_Weathering_Model/`
        * This folder contains code and data that was used to paleogeographically reconstruct large igneous provinces and ophiolites, as well as code that was used to simulate global seawater chemistry.
    * `Tambien_Stratigraphy/`
        * This folder contains code and data that was used to analyze and plot all lithostratigraphic and chemostratigraphic data developed from the Tambien Group.
* `Manuscript/`
    * This folder contains a copy of the manuscript and data repository, alongside raw figure files.
* `pystrat/`
    * This folder contains a snapshot of the Python module pystrat (https://github.com/yuempark/pystrat), which is required to run some of the notebooks within this repository.

## Contact

All code included within this repository should work 'out of the box.' However, if any difficulties are encountered while attempting to execute the code, or if you are seeking some guidance on how to adapt the code for your own study, please do not hesitate to contact any of the authors. Contact information is listed on the manuscript file.
