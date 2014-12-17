# Calculate Zonal Statistics for Overlapping Zones
# ---------------------------------------------------------------------------
# Import system modules

import sys, string, os, gc, shutil
import arcpy
from arcpy import env
from arcpy.sa import *
#enable garbage collection
gc.enable()



# buffer the points

#points_to_buffer = 'C:\\GIS\\PA_predicts_mapping\\all_sites_taxa_split_coarse_06_2014_prj_WGS84_merc.shp'
#output_buffer = 'C:\\GIS\\PA_predicts_mapping\bfr_6_14.shp'

# make sure raster is in the same projection as buffer file





# Get the system TEMP variable value

tempDIR = "C:/GIS/PA_predicts_mapping/temporary"
#create scratch workspace name
# much much quicker when the working folder is in the home folder
tempWorkspace = "tempWorkspaceForStats_acc"

# make temporary workspace if one does not exist
if not arcpy.Exists(tempDIR + os.sep + tempWorkspace):
    arcpy.CreateFolder_management(tempDIR, tempWorkspace)
    
arcpy.env.workspace = tempDIR + os.sep + tempWorkspace


# Check out any necessary licenses
arcpy.CheckOutExtension("Spatial")

# Set geoprocessor object property to overwrite existing output, by default
arcpy.env.overwriteOutput = True

# stop adding outputs to the map
arcpy.env.addOutputsToMap = False



#projected versions
#inputFeatureZone = "C:/GIS/PA_predicts_mapping/bfr_6_14.shp"
#valueRaster = "C:/GIS/PA_predicts_mapping/acc_50k_merc"

inputFeatureZone = "C:/GIS/PA_predicts_mapping/bfr_11_14_moll_1000m.shp"

valueRaster = "C:/GIS/PA_predicts_mapping/acc_50k.tif"
#valueRaster = "C:/GIS/Data/GMTED2010/mn30_grd"
#valueRaster = "C:/GIS/Data/SEDAC_GRUMP/gl_grumpv1_pdens_00_grid_30/gluds00ag"
#valueRaster = "C:/GIS/Data/GAEZ/plate46/plate46_null1"
#valueRaster = "C:\GIS\PA_predicts_mapping\slope30_WGS"

# make a separate feature class for each study

def getValueList (inputTable, field):
    valueList = []
    valueSet = set()
    rows = arcpy.SearchCursor(inputTable)
    for row in rows:
        value = row.getValue(field)
        if value not in valueSet:
            valueList.append(str(value))
        valueSet.add(value)
    valueList.sort()
    return valueList


all_SS = getValueList(inputFeatureZone, "SS")



# Make a layer from the feature class
# need to have layer to do processes below
zone = "zone"
arcpy.MakeFeatureLayer_management(inputFeatureZone, zone)




i = 1

# start loop by creating feature layer for study of interest
# crashes after 150 - 250 iterations if going through individual polygons, so split by study
for n in all_SS:
    obj = "\"SS\" = '%s'" %(n)
    arcpy.SelectLayerByAttribute_management (zone,"NEW_SELECTION",obj)
    studyName = str(n)
    arcpy.MakeFeatureLayer_management(zone, studyName)

    #then use this layer for iterative zonal stats
    #Loop through feature

    all_SSS = getValueList(studyName, "SSS")
    
    for m in all_SSS:

        #select the site to process
        obj_SSS = "\"SSS\" = '%s'" %(m)
        arcpy.SelectLayerByAttribute_management(studyName, "NEW_SELECTION", obj_SSS)
   
        #create a unique name for the zone, by appending number to input name given - use this name later below
        uniqueZone = arcpy.CreateUniqueName("single_zone.shp", arcpy.env.workspace)
        uniqueTable = arcpy.CreateUniqueName("ZSasT.dbf", arcpy.env.workspace)
        #create a temporary feature to use for zonal stats as table - function takes only selected features
        arcpy.CopyFeatures_management(studyName, uniqueZone)
        outZSaT = ZonalStatisticsAsTable(uniqueZone, "SSS", valueRaster, uniqueTable, "DATA", "ALL")

        # delete temporary feature
        arcpy.Delete_management(uniqueZone)

        #combine results for the second iteration onwards
        if m == all_SSS[0]: mergedTable = uniqueTable
        else: arcpy.Append_management(inputs = uniqueTable, target = mergedTable)

        # delete the individual record after second iteration
        if m != all_SSS[0]: arcpy.Delete_management(uniqueTable)

        print i
        i = i+1
        print m

    # delete the study layer
    arcpy.Delete_management(studyName)

    #add studytable to all studies
    if n == all_SS[0]: allStudiesTable = mergedTable
    else: arcpy.Append_management(inputs = mergedTable, target = allStudiesTable)

    # delete the individual record after second iteration
    if n != all_SS[0]: arcpy.Delete_management(mergedTable)

#free memory
del zone

