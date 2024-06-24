# Created by: Ben Gannon
# Created on: 06/22/2022
# Last Updated: 06/24/2024

'''
This script will model Suppression Difficulty Index (SDI) as defined in
Rodriguez y Silva et al. 2020 plus a universal slope modification.

This script was adapted from the earlier work of Matt Panunto, Quresh
Latif, and Pyrologix.

This script assumes that the streets and trails shapefiles and the flame
length (FL) and heat per unit area (HPA) rasters have already been
clipped to the study area.
'''

#-> Import modules
import arcpy, os, sys, traceback, shutil, numpy as np
from arcpy.sa import *

#-> Check out spatial analyst extension
arcpy.CheckOutExtension('Spatial')

#-> Define workspace
workspace = r'E:\SDI_v2\SDI' # Change this to match your directory structure
arcpy.env.scratchWorkspace = workspace + r'\SCRATCH.gdb' # Change this to match your directory structure

#-> Create variables

# Input
settings = np.genfromtxt(r'E:\SDI_v2\SDI\SDI_inputs_ESRI.csv',
                         delimiter=',',skip_header=1,dtype=None,encoding=None)
inDEM = settings[settings[:, 0] == 'DEM'][0,1] # LANDFIRE DEM raster (meters)
inASP = settings[settings[:, 0] == 'ASP'][0,1] # LANDFIRE aspect raster (degrees)
inFBFM40 = settings[settings[:, 0] == 'FBFM40'][0,1] # LANDFIRE fire behavior fuel model raster
inFL = settings[settings[:, 0] == 'FL'][0,1] # FlamMap flame length raster (meters)
inHUA = settings[settings[:, 0] == 'HUA'][0,1] # FlamMap heat per unit area raster (kJ/m2)
SA = settings[settings[:, 0] == 'SA'][0,1] # Study area shapefile for clipping
streets = settings[settings[:, 0] == 'streets'][0,1] # Streets/roads shapefile
trails = settings[settings[:, 0] == 'trails'][0,1] # Trails shapefile
opath = settings[settings[:, 0] == 'opath'][0,1] # SDI raster output path
oname = settings[settings[:, 0] == 'oname'][0,1] # SDI raster output name
outProj = 102039 # USA Contiguous Albers Equal Area Conic USGS Version

#####-----> Main body of analysis
print('####----> Suppression Difficulty Index')
print
arcpy.env.overwriteOutput = True
arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(outProj)

#-> Clip CONUS rasters to study area and create raster objects
# Assumes that FL and HUA are clipped to study area
DEM = arcpy.env.scratchWorkspace + r'\DEM'
arcpy.Clip_management(inDEM,'',DEM,SA)
DEM = Raster(DEM)
ASP = arcpy.env.scratchWorkspace + r'\ASP'
arcpy.Clip_management(inASP,'',ASP,SA)
ASP = Raster(ASP)
FBFM40 = arcpy.env.scratchWorkspace + r'\FBFM40'
arcpy.Clip_management(inFBFM40,'',FBFM40,SA)
FBFM40 = Raster(FBFM40)
FL = Raster(inFL)
HUA = Raster(inHUA)

#-> Set snap and extent environments to DEM
arcpy.env.snapRaster = DEM
arcpy.env.extent = DEM

#-> ENERGY BEHAVIOR SUB-INDEX
print('#-> Calculating Energy Behavior Sub Index')

# Reclassify Flammap flamelength raster values based on Table 2 of Rodriguez y Silva et al
# Convert flame length units to meters before reclassifying, if necessary
FL_recl = Con((FL <= 0.5), 1,
              Con((FL > 0.5) & (FL <= 1), 2,
              Con((FL > 1)   & (FL <= 1.5), 3,
              Con((FL > 1.5) & (FL <= 2), 4,
              Con((FL > 2)   & (FL <= 2.5), 5,
              Con((FL > 2.5) & (FL <= 3), 6,
              Con((FL > 3)   & (FL <= 3.5), 7,
              Con((FL > 3.5) & (FL <= 4), 8,
              Con((FL > 4)   & (FL <= 4.5), 9,
              Con((FL > 4.5), 10))))))))))
print('Reclassified flame lengths to assigned index values')

# Reclassify Flammap heat per unit area raster values based on Table 2 of Rodriguez y Silva et al
# Convert heat per unit area raster from kJ/m2 to kcal/m2
HUA_kcal = HUA * 0.2388459
HUA_recl = Con((HUA_kcal <= 380),1,
               Con((HUA_kcal > 380)  & (HUA_kcal <= 1265), 2,
               Con((HUA_kcal > 1265) & (HUA_kcal <= 1415), 3,
               Con((HUA_kcal > 1415) & (HUA_kcal <= 1610), 4,
               Con((HUA_kcal > 1610) & (HUA_kcal <= 1905), 5,
               Con((HUA_kcal > 1905) & (HUA_kcal <= 2190), 6,
               Con((HUA_kcal > 2190) & (HUA_kcal <= 4500), 7,
               Con((HUA_kcal > 4500) & (HUA_kcal <= 6630), 8,
               Con((HUA_kcal > 6630) & (HUA_kcal <= 8000), 9,
               Con((HUA_kcal > 8000), 10))))))))))
print('Reclassified heat per unit area to assigned index values')

# Calculate Energy Behavior raster
# The focal statistics smoothing operation replaces the previous propFuel multplication factor
eb_rast = (2 * Float(FL_recl) * Float(HUA_recl) / (Float(FL_recl) + Float(HUA_recl)))
print('Created energy behavior raster')
# Smooth, but set non-burnable pixels to value of 1
neighborhood = NbrCircle(56.41896,'MAP')
eb_smooth = Con(FBFM40 > 100,FocalStatistics(eb_rast,neighborhood,'MEAN','DATA'),1)
#eb_smooth.save(arcpy.env.scratchWorkspace + r'\energy_behavior')
print('Smoothed energy behavior raster')
print

#-> ACCESSIBILITY SUB INDEX
# Run Euclidean Distance on the streets layer
# Classify areas within 100 m of a road as a 10, decreasing to 1 where the nearest road is >900 m away
print('#-> Calculating accessibility sub index')
streets_lyr = arcpy.MakeFeatureLayer_management(streets,'streets_lyr')
arcpy.SelectLayerByAttribute_management(streets_lyr,'NEW_SELECTION',where_clause='"FuncSpdCat" < 8')
print('Selected roads with speedclass < 8')
roaddist = EucDistance(streets_lyr,'',30)
print('Calculated distance to roads')
accessibility = Con((roaddist <= 100),10,
                Con((roaddist > 100) & (roaddist <= 200), 9,
                Con((roaddist > 200) & (roaddist <= 300), 8,
                Con((roaddist > 300) & (roaddist <= 400), 7,
                Con((roaddist > 400) & (roaddist <= 500), 6,
                Con((roaddist > 500) & (roaddist <= 600), 5,
                Con((roaddist > 600) & (roaddist <= 700), 4,
                Con((roaddist > 700) & (roaddist <= 800), 3,
                Con((roaddist > 800) & (roaddist <= 900), 2,
                1)))))))))
#accessibility.save(arcpy.env.scratchWorkspace + r'\accessibility')
print('Reclassified distance to roads raster to accessibility index value')
print
   
#-> MOBILITY SUB-INDEX
# Skipping, setting mobility = 1 in calculation

#-> PENETRABILITY SUB-INDEX
print('#-> Calculating penetrability sub index')

# Convert Aspect Raster (degrees) into assigned values
# Changed from the paper after 2018 Cordoba visit
aspect_class = Con((ASP >= 337.5), 10, # North Facing
                   Con((ASP >= -1)    & (ASP < 22.5), 10, # North Facing
                   Con((ASP >= 22.5)  & (ASP < 67.5), 8,  # Northeast Facing
                   Con((ASP >= 292.5) & (ASP < 337.5), 7, # Northwest Facing
                   Con((ASP >= 67.5)  & (ASP < 112.5), 6, # East Facing
                   Con((ASP >= 247.5) & (ASP < 292.5), 5, # West Facing
                   Con((ASP >= 112.5) & (ASP < 157.5), 4, # Southeast Facing
                   Con((ASP >= 202.5) & (ASP < 247.5), 3, # Southwest Facing
                   Con((ASP >= 157.5) & (ASP < 202.5), 1))))))))) # South Facing 
print('Converted aspect raster into assigned values')

# Calculate Percent Slope
slope = Slope(DEM,'PERCENT_RISE')
print('Calculated percent slope')

# Convert Percent Slope Raster into assigned values
slope_class = Con((slope >= 0) & (slope < 6), 10,
                  Con((slope >= 6) & (slope < 11), 9,
                  Con((slope >= 11) & (slope < 16), 8,
                  Con((slope >= 16) & (slope < 21), 7,
                  Con((slope >= 21) & (slope < 26), 6,
                  Con((slope >= 26) & (slope < 31), 5,
                  Con((slope >= 31) & (slope < 36), 4,
                  Con((slope >= 36) & (slope < 41), 3,
                  Con((slope >= 41) & (slope < 46), 2,
                  Con((slope >= 46), 1))))))))))
print('Converted percent slope into assigned values')

# Reclassify fuels into assigned RTC values using reclass table
# This will convert the fuel types into a fuel control difficulty weight based on estimated
# fireline production rates (meters per mile) for a 20-person crew. Note that this is an
# approximation of the original table in Rodriguez y Silva et al 2014 for Scott and Burgan
# (2005) fuel models.
rtcReclass = ('101 10;102 10;103 10;104 10;105 10;106 10;107 3;108 3;109 3;'
             '121 10;122 10;123 10;124 3;'
             '141 4;142 4;143 3;144 3;145 2;146 4;147 2;148 2;149 2;'
             '161 8;162 8;163 4;164 4;165 4;'
             '181 8;182 10;183 8;184 4;185 2;186 10;187 2;188 8;189 8;'
             '201 4;202 4;203 1;204 1;205 1')
fuel_cntrl = Reclassify(FBFM40,'VALUE',rtcReclass,'NODATA')
#fuel_cntrl = ReclassByASCIIFile(FBFM40,rtc_lookup,'NODATA')
print('Reclassified fuels into assigned fireline production values')

# Merge streets speed cat 8 with trails
streets_lyr = arcpy.MakeFeatureLayer_management(streets,'streets_lyr')
arcpy.SelectLayerByAttribute_management(streets_lyr,'NEW_SELECTION',where_clause='"FuncSpdCat" > 7')
trails_streets = arcpy.env.scratchWorkspace + r'\trails_streets'
arcpy.Merge_management([trails,streets_lyr],trails_streets)
print('Merged street speed cat 8 with trails')

# Calculate length of trails within 1 ha (56.41896m radius) moving window.
trail_len = LineStatistics(trails_streets,'NONE','30','56.41896','LENGTH')
print('Calculated length of trails within 1 ha')

# Convert pre-suppression trails raster into assigned values
trail_len_class = Con((trail_len >= 0) & (trail_len < 10), 1,
                      Con((trail_len >= 10) & (trail_len < 20), 2,
                      Con((trail_len >= 20) & (trail_len < 30), 3,
                      Con((trail_len >= 30) & (trail_len < 40), 4,
                      Con((trail_len >= 40) & (trail_len < 50), 5,
                      Con((trail_len >= 50) & (trail_len < 60), 6,
                      Con((trail_len >= 60) & (trail_len < 70), 7,
                      Con((trail_len >= 70) & (trail_len < 80), 8,
                      Con((trail_len >= 80) & (trail_len < 90), 9,
                      Con((trail_len >= 90), 10))))))))))
print('Reclassified road/trail lengths to assigned values')

# Create raster for Penetrability Sub Index
# The focal statistics smoothing operation replaces the previous propFuel multplication factor
penetrability = (Float(slope_class) + Float(fuel_cntrl) + Float(aspect_class) + (2*Float(trail_len_class)))/5
print('Calculated penetrability index')
neighborhood = NbrCircle(56.41896,'MAP')
pen_smooth = FocalStatistics(penetrability,neighborhood,'MEAN','DATA')
#pen_smooth.save(arcpy.env.scratchWorkspace + r'\penetrability')
print('Smoothed penetrability index')
print

#-> FIRELINE OPENING/CREATION SUB-INDEX
print('#-> Calculating fireline opening/creation sub index')

# Create Slope Adjustment Raster for hand work
slope_adjust_hand = Con((slope >= 0) & (slope < 16), 1,
                        Con((slope >= 16) & (slope < 31), 0.8,
                        Con((slope >= 31) & (slope < 46), 0.6,
                        Con((slope >= 46), 0.5))))
print('Calculated slope adjustment raster for hand work')

# Create Slope Adjustment Raster for machine work
slope_adjust_mach = Con((slope >= 0) & (slope < 16), 1,
                        Con((slope >= 16) & (slope < 21), 0.8,
                        Con((slope >= 21) & (slope < 26), 0.7,
                        Con((slope >= 26) & (slope < 31), 0.6,
                        Con((slope >= 31) & (slope < 36), 0.5,
                        Con((slope >= 36), 0))))))
print('Calculated slope adjustment raster for machine work')

# Calculate Fireline Opening Sub Index
firelineopening = (Float(fuel_cntrl)*slope_adjust_hand) + (Float(fuel_cntrl)*slope_adjust_mach)
#firelineopening.save(arcpy.env.scratchWorkspace + r'\firelineopening')
print('Calculated fireline opening index')
print

#-> Kit O'Connor (KO) UNIVERSAL SLOPE MOBILITY HAZARD ADJUSTMENTS
print('#-> Calculating KO universal slope mobility hazard adjustment')

# Create slope hazard raster    
slope_haz = Con((slope >= 0) & (slope < 10), 0.01,
                Con((slope >= 10) & (slope < 20), 0.1,
                Con((slope >= 20) & (slope < 30), 0.2,
                Con((slope >= 30) & (slope < 40), 0.4,
                Con((slope >= 40) & (slope < 50), 0.6,   
                Con((slope >= 50) , 0.8,))))))
print('Calculated slope mobility hazard')
print
    
#-> SDI CALCULATION
print('#-> Calculing SDI from sub-indices')

# Create SDI raster without slope modification
# mobility = 1, 2020 publication version
SDI = Float(eb_smooth) / (Float(accessibility) + Float(1) + Float(pen_smooth) + Float(firelineopening))
SDI = Con(IsNull(SDI) & (FBFM40 < 100),Float(0),SDI) # Setting non-burnable pixels to zero
print('Calculated SDI without slope modification')

# Optional: save denominator for inspection
#sdi_denom = (Float(accessibility) + Float(1) + Float(pen_smooth) + Float(firelineopening))
#sdi_denom.save(arcpy.env.scratchWorkspace + r'\denominator')

# Create SDI raster with universal slope modification
SDI_SLP = SDI + slope_haz
print('Calculated SDI with slope modification')

# Set water to zero SDI
SDI_SLP_H2O = Con(FBFM40==98,Float(0),SDI_SLP)
print('Set water to zero SDI')
    
# Combine using base model for roads & trails, and slope modification for everything else
trails_rds = arcpy.env.scratchWorkspace + r'\trails_st_all'
arcpy.Merge_management([trails,streets],trails_rds)
arcpy.AddField_management(trails_rds,'ZONE','SHORT')
arcpy.CalculateField_management(trails_rds,'ZONE','1','PYTHON_9.3','')
trails_rds_rast = arcpy.env.scratchWorkspace + r'\trails_st_all_rast'
arcpy.FeatureToRaster_conversion(trails_rds,'ZONE',trails_rds_rast,SDI_SLP)
fSDI = Con(IsNull(trails_rds_rast),SDI_SLP_H2O,SDI)
print('Combined base and slope adjusted SDI for roaded/trailed or unroaded/untrailed areas')
print

# Save as integer
fSDI100 = Int(fSDI*Float(100))
fSDI100.save(opath + '\\' + oname)
print('Saved as integer')
print

print('Script complete!')
