Description
---
Breast cancer data are downloaded from the cBioPortal for Cancer Genomics (http://www.cbioportal.org) [1], including somatic mutation (n = 2,369), EXP (n = 1,904), and CNA (n = 2,173). It contains the METABRIC breast cancer cohort assembled from 2509 primary breast cancer patients with 548 matched normals in the United Kingdom and Canada [2]. 

The files here contain:
1. Link_to_download_full_data.txt: A text file has a shared google drive link to download all the raw data mentioned.

2. raw_data.zip: a subset of the raw data nescessary for running examples from the full one above.

3. process_data_DrGA.R: R code examples to process the raw data from scratch before inputting them to DrGA.

4. code_DrGA.RData:  Processed data are available to run DrGA immediately.

Reference
---
```sh
1. Cerami, E. et al. The cBio cancer genomics portal: an open platform for exploring multidimensional cancer genomics data. Cancer Discov. 2, 401–404. https://doi.org/10.1158/2159-8290.cd-12-0095 (2012).

2. Pereira, B. et al. The somatic mutation profiles of 2,433 breast cancers refine their genomic and transcriptomic landscapes. Nat. Commun. 7, 11479. https://doi.org/10.1038/ncomms11479 (2016).
```
