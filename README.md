# glistin
This set of codes is a toolbox to download and resample GLISTIN-A radar interferometer digital elevation models from NASA's Oceans Melting Greenland campaign. To use these codes,
1. Determine a directory on your local machine where the Level 1 data and resampled output products will be stored. This directory is named "dataFolder" in the scripts provided.
2. Determine which GLISTIN-A DEMs will be downloaded. Refer to Fenty_et_al_2019 for a map indicating the location of swath indices.
3. Download the Level 1 data from NASA's UAVSAR website. Note that the resampling script relies on a particular directory structure to locate the Level 1 data on the local machine within "dataFolder". Please use the download_level_1_data.py code to ensure this structure is properly constructed.
4. Resample (via bin-averaging) the DEM data to a resolution of your choosing using resample_pipeline.py.



To download data: use download_level_1_data.py to download level 1 data from https://uavsar.jpl.nasa.gov/

Keywords:
- dataFolder (-d): Directory where data will be downloaded. (Required)
- indexList (-i): List of swath indices to download. (Optional. Default is to download all 81 swaths)
- yearList (-y): List of years to download each swath. (Optional. Default is to download swaths in all available years (2016-2019))
- downloadData (-a): Choose 1 if you would like to download data, otherwise choose 0. (Optional. Default is download data.)
- downloadMetadata (-m): Choose 1 if you would like to download metadata, otherwise choose 0. (Optional. Default is download metadata.)



To resample data: use resample_pipeline.py

Keywords:
- dataFolder (-d): Directory where data will be downloaded. (Required)
- resolution (-r): Resolution for the resampling. Default is 50 m.
- indexList (-i): List of swath indices to resample. (Optional. Default is to resample all 81 swaths)
- yearList (-y): List of years to resample each swath. (Optional. Default is to download swaths in all available years (2016-2019))
