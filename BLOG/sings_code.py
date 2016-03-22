#Imports
from __future__ import division, print_function

import numpy
from astropy import wcs
from astropy.io import fits
import sys
import os
import arcpy

#Environment Variables
arcpy.env.pyramid = "NONE"
arcpy.env.rasterStatistics  = "NONE"

#Paths and Filenames
output_image_path = "C:/PROJECTS/R&D/ASTROARC/SINGS/TIFF/"
input_image_path = "C:/PROJECTS/R&D/ASTROARC/SINGS/FITS/"
composite_image_path = output_image_path + "COMPOSITES"
num_bands = 4
instrument = "IRAC"
galaxy = "ngc5194"
composite_filename = galaxy + "_" + instrument + "_" + str(num_bands) + "Bands.tif"
composite_path_and_name = os.path.join(composite_image_path, composite_filename)

#Make directories if they don't exist
if not os.path.exists(output_image_path):
    os.mkdir(output_image_path)

if not os.path.exists(composite_image_path):
    os.mkdir(composite_image_path)

#For all four images per galaxy
for num in range(1,num_bands+1):
    input_filename = galaxy + "_v7.phot." + str(num)+ ".fits"
    temp_filename = input_filename[:-4] + 'rev.tif'
    output_filename = input_filename[:-4] + 'rev.WCS.tif'
    temp_file_and_path = os.path.join(output_image_path, temp_filename)
    output_file_and_path = os.path.join(output_image_path, output_filename)
    input_file_and_path = os.path.join(input_image_path, input_filename)

    print('Processing file: ' + input_filename)

    #Read the image and header data
    hdulist = fits.open(input_file_and_path)
    w = wcs.WCS(hdulist[0].header)
    image = hdulist[0].data

    #Some pixel coordinates of interest.
    pixcrd = numpy.array([[0, 0], [0, image.shape[0]],[image.shape[1], 0],
            [image.shape[1], image.shape[0]]], numpy.float_)
    #Convert to string
    pixcrd_string = '0 0;' + '0 ' + str(image.shape[0]) + ';' + str(image.shape[1]) + ' 0;' + str(image.shape[1]) + ' ' + str(image.shape[0])

    # Convert pixel coordinates to world coordinates
    world = w.wcs_pix2world(pixcrd, 1)
    #Convert to string to project/warp in arcpy
    coordinate_string = str(-1*world[0,0]) + ' ' + str(world[0,1]) + ';' + str(-1*world[1,0]) + ' ' + str(world[1,1]) + ';' + str(-1*world[2,0]) + ' ' + str(world[2,1]) + ';' + str(-1*world[3,0]) + ' ' + str(world[3,1])

    # Convert the same coordinates back to pixel coordinates.
    pixcrd2 = w.wcs_world2pix(world, 1)

    print('Converting numpy array to TIFF.')

    #Converting raster to TIFF image
    assert numpy.max(numpy.abs(pixcrd - pixcrd2)) < 1e-6
    #Reversing the array
    flipim = numpy.flipud(image)
    #Writing the array to a raster
    double_image = numpy.float64(flipim)
    myRaster = arcpy.NumPyArrayToRaster(double_image)
    myRaster.save(temp_file_and_path)
    #Giving raster the sky coordinates
    source_pnt = pixcrd_string
    target_pnt = coordinate_string
    arcpy.Warp_management(temp_file_and_path, source_pnt, target_pnt,
            output_file_and_path, "POLYORDER1","BILINEAR")
    #Deleting the unreference image
    arcpy.Delete_management(temp_file_and_path)
    print(output_file_and_path)

#Composite bands
arcpy.env.workspace = output_image_path
rasters = arcpy.ListRasters(galaxy + "*", "TIF")
arcpy.CompositeBands_management(rasters,composite_path_and_name)

#Create the mosaic and add images to it
coordinate_sys = "PROJCS['WGS_1984_Web_Mercator_Auxiliary_Sphere',GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Mercator_Auxiliary_Sphere'],PARAMETER['False_Easting',0.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',0.0],PARAMETER['Standard_Parallel_1',0.0],PARAMETER['Auxiliary_Sphere_Type',0.0],UNIT['Meter',1.0]];-20037700 -30241100 10000;-100000 10000;-100000 10000;0.001;0.001;0.001;IsHighPrecision"
imagery_spatial_ref = "GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]];-400 -400 1000000000;-100000 10000;-100000 10000;8.98315284119522E-09;0.001;0.001;IsHighPrecision"
mosaic_gdb = r"C:\PROJECTS\R&D\ASTROARC\SINGS\Spitzer.gdb"
mosaic_name = "SINGS"
mosaic_dataset = os.path.join(mosaic_gdb, mosaic_name)
arcpy.CreateMosaicDataset_management(mosaic_gdb, mosaic_name, coordinate_sys,
        num_bands="", pixel_type="", product_definition="NONE",
        product_band_definitions="")
arcpy.AddRastersToMosaicDataset_management(mosaic_dataset, "Raster Dataset",
        composite_image_path, "UPDATE_CELL_SIZES", "UPDATE_BOUNDARY",
        "NO_OVERVIEWS", "", "0", "1500", imagery_spatial_ref, "#", "SUBFOLDERS",
        "ALLOW_DUPLICATES", "NO_PYRAMIDS", "NO_STATISTICS",
        "NO_THUMBNAILS", "#", "NO_FORCE_SPATIAL_REFERENCE")

print("Done.")