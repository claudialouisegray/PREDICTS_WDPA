
# Import arcpy module
import arcpy


# Local variables:
shapefile_from_R = "C:/GIS/PA_predicts_mapping/PA_11_2014.shp"
projected_points = "C:\\GIS\\PA_predicts_mapping\\PA_11_2014_moll.shp"
output = "C:\\GIS\\PA_predicts_mapping\\bfr_11_14_moll_1000m_test.shp"

# Process: Define Projection
arcpy.DefineProjection_management(shapefile_from_R, "GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]]")

# Process: Project to mollweide
arcpy.Project_management(shapefile_from_R, projected_points, "PROJCS['WGS_1984_World_Mercator',GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Mercator'],PARAMETER['False_Easting',0.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',0.0],PARAMETER['Standard_Parallel_1',0.0],UNIT['Meter',1.0]]", "", "GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]]")

#buffer
arcpy.Buffer_analysis(projected_points,output,"1000 Meters","FULL","ROUND","NONE","#")
