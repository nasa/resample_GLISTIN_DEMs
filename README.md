# glistin
This set of codes is a toolbox to download and resample GLISTIN-A radar interferometer digital elevation models from NASA's Oceans Melting Greenland campaign. To use these codes,
1. Determine a directory on your local machine where the Level 1 data and resampled output products will be stored. This directory is named "dataFolder" in the scripts provided.
2. Determine which GLISTIN-A DEMs will be downloaded. Refer to Fenty_et_al_2019 for a map indicating the location of swath indices.
3. Download the Level 1 data from NASA's UAVSAR website. Note that the resampling script relies on a particular directory structure to locate the Level 1 data on the local machine within "dataFolder". Please use the download_level_1_data.py code to ensure this structure is properly constructed.
4. Resample (via bin-averaging) the DEM data to a resolution of your choosing using resample_pipeline.py.
