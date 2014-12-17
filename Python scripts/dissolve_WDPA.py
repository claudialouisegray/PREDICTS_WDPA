

import sys, string, os, gc, shutil
import arcpy
from arcpy import env
from arcpy.sa import *
#enable garbage collection
gc.enable()


#spatial join to country
# select by country
#dissolve
#get all files
#merge
#dissolve



# Get the system TEMP variable value

tempDIR = "C:/GIS/PA_predicts_mapping"
#create scratch workspace name
# much much quicker when the working folder is in the home folder
tempWorkspace = "WDPA_dissolve1"

# make temporary workspace if one does not exist
if not arcpy.Exists(tempDIR + os.sep + tempWorkspace):
    arcpy.CreateFolder_management(tempDIR, tempWorkspace)
    
arcpy.env.workspace = tempDIR + os.sep + tempWorkspace


# Check out any necessary licenses
arcpy.CheckOutExtension("Spatial")

#load WDPA layer
in_features = "C:/GIS/Data/WDPA/WDPAmollmergeJULY/WDPAmollmergeJULY_no_biosphere_inscribed_adopted.shp"




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
ISO3 = getValueList(in_features, "ISO3")

#go through whole WDPA in a loop, selecting each level of ISO and dissolving that


i = 1

#make a layer of the WDPA polygons so that it can be selected
in_features_layer = "in_features_layer"
arcpy.MakeFeatureLayer_management(in_features, in_features_layer)

# start from what I already have
#ISO3_incomplete = ISO3[44:len(ISO3)-1]

for n in ISO3:
    obj = "\"ISO3\" = '%s'" %(n)
    arcpy.SelectLayerByAttribute_management (in_features_layer,"NEW_SELECTION",obj)
    countryCode = str(n)
    arcpy.MakeFeatureLayer_management(in_features_layer, countryCode)

    # repair this layer if needed
    arcpy.RepairGeometry_management(countryCode)
    
    #then dissolve this layer
    dissolve_name =  arcpy.CreateUniqueName(str(n), arcpy.env.workspace)
    #dissolve_name =  arcpy.CreateUniqueName("dissolve.shp", arcpy.env.workspace)
    #arcpy.Dissolve_management (countryCode, dissolve_name, "ISO3")

    arcpy.Dissolve_management(countryCode, dissolve_name,"#","ISO3 FIRST","SINGLE_PART","DISSOLVE_LINES")

    print i
    i = i+1
    print n

    # delete the study layer
    arcpy.Delete_management(countryCode)

    # add polygon to the first one
    #if n == ISO3[0]: all_dissolves = dissolve_name
    #else: arcpy.Append_management(inputs = dissolve_name, target = all_dissolves)


# list all files in working directory
all_files = arcpy.ListFeatureClasses()

# go through and dissolve each again to get rid of novel divides not present before
#for a in all_files:
#    dissolve_name2 =  "C:/GIS/PA_predicts_mapping/WDPA_2nd_dissolve/" + a
#    file_name = "C:/GIS/PA_predicts_mapping/WDPA_dissolve/" + a
#    arcpy.RepairGeometry_management(file_name)
#    arcpy.Dissolve_management(file_name, dissolve_name2,"#","#","SINGLE_PART","DISSOLVE_LINES")

    
# merge all 256 flattened shapefiles together
arcpy.env.workspace = "C:/GIS/PA_predicts_mapping/WDPA_dissolve1"
all_files2 = arcpy.ListFeatureClasses()
WDPA_moll_flattened = "C:/GIS/Data/WDPA/WDPAmollmergeJULY/WDPA_mollJuly14_no_biosphere_inscribed_adopted_flattened.shp"
arcpy.Merge_management(all_files2, WDPA_moll_flattened)

    
# repair geometry of this layer incase of errors with dissolve
arcpy.RepairGeometry_management(WDPA_moll_flattened)

# dissolve combined layer incase any PAs still overlap... new breaks are put in, not acceptable
out_feature_class = "C:/GIS/Data/WDPA/WDPAmollmergeJULY/WDPAmollJuly14_dissolved.shp"
arcpy.Dissolve_management (WDPA_moll_flattened, out_feature_class,"","", "SINGLE_PART","DISSOLVE_LINES")

#try making all one layer (can use this for dist to nearest PA)
#doesnt work, not enough memory
#out_feature_class = "C:/GIS/Data/WDPA/WDPAmollmergeJULY/WDPAmollJuly14_dissolve_multipart.shp"
#arcpy.Dissolve_management (WDPA_moll_flattened, out_feature_class,"","", "MULTI_PART","DISSOLVE_LINES")



