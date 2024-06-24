# Suppression Difficulty Index (SDI)

Suppression Difficulty Index (SDI) is a spatial decision support product developed by a group of wildland firefighting experts in Spain (Rodríguez y Silva et al. 2014, 2020) that has been adapted for use in the United States by the USDA Forest Service Rocky Mountain Research Station. The formulation for SDI is premised on the concept that fire behavior acts as a driving force in the numerator and conditions that facilitate suppression act as resisting forces in the denominator (Equation 1). The numerator, also called the energy behavior index (_Ice_), grows in response to conditions that increase fire behavior, so SDI will increase where and when fuel, topography, and weather conditions promote extreme fire behaviors. The denominator sums the values for sub-indices representing how accessibility (_Ia_), penetrability (_Ip_), and fireline creation (_Ic_) influence suppression opportunities across the landscape. Higher values for sub-indices in the denominator indicate better conditions for suppression, which lowers SDI for a given energy behavior index value in the numerator.   

**Equation 1**  _SDI = Ice/(Ia + Ip + Ic + 1)_

Feedback from early users in the United States suggested that SDI underestimated the difficulty of suppressing fire on steep slopes, especially ones with sparse fuels that were mapped as non-burnable by LANDFIRE. While steep slopes with little fuel (e.g., scree slopes, cliff bands, etc.) may act as natural barriers to fire, these areas are not ideal locations to deploy hand crews or machines to construct containment line. Starting in the 2021 season, an additional slope factor was applied to increase SDI on steep slopes regardless of fire behavior or other operational factors.

**Equation 2**  _SDI = Ice/(Ia + Ip + Ic + 1) + Slope Factor_

SDI provides a snapshot of how potential fire behavior and operational factors combine to influence suppression difficulty across broad landscapes. SDI tends to align well with firefighter expectations for low suppression difficulty near major roads, significant waterbodies, and flatter terrain, and for high suppression difficulty in steep, forested areas without recent history of wildfire or fuel treatment. SDI is moderately sensitive to the fire behavior scenario used to model the energy behavior index (numerator); SDI will rise with increasingly extreme fire weather inputs, but the spatial patterns and their interpretations for fire management do not drastically change across the range of typical large fire weather conditions. Hence, it is usually not necessary to model SDI for incident specific fire weather conditions.

Full details on the standard SDI calculation methods can be found in Risk Management Assistance (RMA) SharePoint site.

**Instructions for Use**

Two methods are provided to calculate SDI: 1) a Python script that uses ESRI's arcpy module for geoprocessing and 2) an R script that uses the terra package for geoprocessing. The SDI Python script was developed and tested using the arcpy module associated with ArcGIS Pro 3.1.4 (ESRI 2023) and Python version 3.9.16 (Python Software Foundation 2023). The SDI R script was developed and tested using the terra package version 1.7.39 (Hijmans 2023) in the R language for statistics and computing version 4.3.1 (R Core Team 2023). The scripts do not include advanced error handling. The end user is responsible for verifying that their inputs are of the correct type, format, and units. The only lines in the scripts that a user should need to modify are the working directory paths. Change the working directory path to match the location of the script and inputs table (SDI inputs csv).

Model inputs are specified in a table that will be read in by the script including:
1) Path to a digital elevation model (DEM) in meters
2) Path to a categorical fire behavior fuel model 40 (FBFM40) raster using the Scott and Burgan (2005) classes
3) Path to a flame length (FL) raster in meters
4) Path to heat per unit area (HUA) raster in kilojoules per square meter
5) Path to study area (SA) shapefile to analyze - should include at least a 1-km buffer to minimize edge effects
6) Path to streets/roads shapefile - assumed to be clipped to the study area
7) Path to the trails shapefile - assumed to be clipped to the study area
8) Output directory path - must exist
9) Output file name for final output

For Python, all file paths should use backward slashes, e.g., "D:\SDI_INPUTS\LA20_Elev_220.tif". For R, all file paths should use forward slashes, e.g., "D:/SDI_INPUTS/LA20_Elev_220.tif". It is assumed that all raster inputs have the same extent, resolution, cell alignment, and spatial projection. Raster reprojection and resampling is not automatically completed in the script - this is considered a data preparation task for the user.

The output includes rasters for each sub-index, the numerator, and Suppression Difficulty Index. 
