

import sys, string, os, gc, shutil
import arcpy
from arcpy import env
from arcpy.sa import *
#enable garbage collection
gc.enable()




# Get the system TEMP variable value

tempDIR = "C:/GIS/PA_predicts_mapping"
# create scratch workspace name
# much much quicker when the working folder is in the home folder
tempWorkspace = "dist_to_boundary_same_PH_11_14"

# make temporary workspace if one does not exist
if not arcpy.Exists(tempDIR + os.sep + tempWorkspace):
    arcpy.CreateFolder_management(tempDIR, tempWorkspace)
    
arcpy.env.workspace = tempDIR + os.sep + tempWorkspace


# Check out any necessary licenses
arcpy.CheckOutExtension("Spatial")


# load predicts matched landuse sites
#or ALL sites (needed for PA analysis)
in_features = "C:/GIS/PA_predicts_mapping/PA_11_2014_moll.shp"



#make a layer of the WDPA polygons so that it can be selected
in_features_layer = "in_features_layer"
arcpy.MakeFeatureLayer_management(in_features, in_features_layer)


#load PA layer
WDPA_moll_flattened = "C:/GIS/Data/WDPA/WDPAmollmergeJULY/WDPAmollmergeJuly_no_biosphere_inscribed_adopted_dissolved_CHECKED.shp"

#make a layer of the WDPA polygons so that it can be selected
WDPA_layer = "WDPA_layer"
arcpy.MakeFeatureLayer_management(WDPA_moll_flattened, WDPA_layer)


#define function to get unique values for a particular field
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


#get a list of each habitat within each study
PH = getValueList(in_features_layer, "Prdmnn_")


# select the points from a study
# make into layer



for n in PH:
    obj = "\"Prdmnn_\" = '%s'" %(n)
    arcpy.SelectLayerByAttribute_management(in_features_layer,"NEW_SELECTION",obj)
    PH_layer = str(n)
    arcpy.MakeFeatureLayer_management(in_features_layer, PH_layer)

# select from WDPA by intersection with that layer

# Replace a layer/table view name with a path to a dataset (which can be a layer file) or create the layer/table view within the script
# The following inputs are layers or table views: "WDPAmollJuly14_dissolved_merged", "matched_landuse_taxa_split_coarse_07_2014"
    arcpy.SelectLayerByLocation_management(WDPA_layer,"INTERSECT",PH_layer,"#","NEW_SELECTION")
    PA_layer = "PA_" + n
    arcpy.MakeFeatureLayer_management(WDPA_layer, PA_layer)

    # make polyline so that distance inside PA are calculated too
    PA_lines = arcpy.CreateUniqueName("PAlines_" + n + ".shp", arcpy.env.workspace)
   
    arcpy.PolygonToLine_management(PA_layer, PA_lines, "IGNORE_NEIGHBORS")

    # get dist to boundary for all points
    arcpy.Near_analysis(PH_layer, PA_lines, "500000 Meters")
    
    #near_table = "C:/GIS/PA_predicts_mapping/dist_to_boundary/table_"+n+".shp" 
    #arcpy.GenerateNearTable_analysis(SS_layer, PA_lines, near_table, "500000 Meters")

    #join to layer
    #arcpy.AddJoin_management(PH_layer, "FID", near_table, "IN_FID", "KEEP_ALL")

    #save this information
    #arcpy.CreateFeatureclass_management(arcpy.env.workspace, "points_"+n+".shp", "", SS_layer) #doesnt work
    arcpy.FeatureClassToShapefile_conversion (PH_layer, arcpy.env.workspace)

    print(PH_layer)
    
    # delete the temporary layers
    arcpy.Delete_management(PH_layer)
    arcpy.Delete_management(PA_layer)
    arcpy.Delete_management(PA_lines)

    
# combine files with dist
# list all files in working directory
all_files = arcpy.ListFeatureClasses()

all_dists = "C:/GIS/PA_predicts_mapping/dist_to_boundary_same_PH_11_14/all_dists.shp"
arcpy.Merge_management(all_files, all_dists)










