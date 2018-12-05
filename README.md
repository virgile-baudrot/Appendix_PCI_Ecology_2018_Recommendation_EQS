[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1972932.svg)](https://doi.org/10.5281/zenodo.1972932)

# Appendix_PCI_Ecology_2018_Recommendation_EQS

Data and Script used in the article recommended in [PCI Ecology 2018](---)

Authors: Virgile Baudrot and Sandrine Charles

**Recommendations to address uncertainties in environmental risk assessment using toxicokinetics-toxicodynamics models**. 

[DOI: 10.24072/pci.ecology.100007](10.24072/pci.ecology.100007)

Also on [BioRxiv](https://www.biorxiv.org/content/early/2018/11/22/356469)

# Reproduction of the numerical analysis

1. Data and their collection are desribed in the manuscript section `Material and methods` -> `Data from experimental toxicity tests`.
All data set are in the files: `data_Car_cst.rda`, `data_Car_var.rda`, `data_Cyp_cst.rda`, `data_Cyp_var.rda`, `data_Dim_cst.rda`, `data_Dim_var.rda`, `data_Mal_cst.rda`, `data_Mal_var.rda`, `data_PRZ_cst.rda`, `data_PRZ_var.rda`.

2. Every computing and figures can be done using scripts in file (other R script are called from this file): [Supplementary_Material_Recommendation_EQS.Rmd](Supplementary_Material_Recommendation_EQS.Rmd)

Note that we `save` and `load` all along the script since computing can be long, and object heavy to keep in RAM.

Feel free to ask (by using [issues](https://github.com/virgile-baudrot/Appendix_PCI_Ecology_2018_Recommendation_EQS/issues) for instance) if any part of the code is not well working.


