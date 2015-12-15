from __future__ import division, print_function

import numpy
from astropy import wcs
from astropy.io import fits
#import h5py
import arcpy
import os
from os import listdir
from os.path import isfile, join

#Define where files are and where they are going
output_path = "C:\\PROJECTS\\R&D\\ASTROARC\\SCOSMOS\\TIFF"
input_path = "C:\\PROJECTS\\R&D\\ASTROARC\\SCOSMOS\\FITS"
composite_path = output_path + "\\" + "COMPOSITES"
numBands = 3
program = 'SCOSMOS'

if not os.path.exists(output_path):
    os.makedirs(output_path)
if not os.path.exists(composite_path):
    os.makedirs(composite_path)

#Get list of files that will be converted
fits_files = [ f for f in listdir(input_path) if isfile(join(input_path,f)) ]

print(fits_files)

counter = 1
#For every file...
for fits_file in fits_files:
    
    print("Processing file: " + str(counter) + ", " + fits_file + "...")
    counter = counter + 1

    filename = fits_file
    output_filename = filename[:-4] + 'rev.tif'
    geo_output_filename = filename[:-4] + 'rev.WCS.tif'
    outputfile_and_path = output_path + "\\" + output_filename
    geo_file_and_path = output_path + "\\" + geo_output_filename
    file_and_path = input_path + "\\" + filename
    composite_filename = composite_path + "\\" + program + "_" + str(numBands) + "Bands.tif"
    
    if not os.path.exists(geo_file_and_path):    
    
        print(outputfile_and_path)
        print(geo_file_and_path)
        hdulist = fits.open(file_and_path)
        
        w = wcs.WCS(hdulist[0].header)
        image = hdulist[0].data


        # Some pixel coordinates of interest.
        pixcrd = numpy.array([[0, 0], [0, image.shape[0]],[image.shape[1], 0], [image.shape[1], image.shape[0]]], numpy.float_)
        pixcrd_string = '0 0;' + '0 ' + str(image.shape[0]) + ';' +  str(image.shape[1]) + ' 0;' + str(image.shape[1]) + ' ' + str(image.shape[0])
        #print(pixcrd_string)

        # Convert pixel coordinates to world coordinates
        # The second argument is "origin" -- in this case we're declaring we
        # have 1-based (Fortran-like) coordinates.
        world = w.wcs_pix2world(pixcrd, 1)
        coordinate_string = str(-1*world[0,0]) + ' ' + str(world[0,1]) + ';' + str(-1*world[1,0]) + ' ' + str(world[1,1]) + ';' + str(-1*world[2,0]) + ' ' + str(world[2,1]) + ';' + str(-1*world[3,0]) + ' ' + str(world[3,1]) 

        # Convert the same coordinates back to pixel coordinates.
        pixcrd2 = w.wcs_world2pix(world, 1)

        # These should be the same as the original pixel coordinates, modulo
        # some floating-point error.
        assert numpy.max(numpy.abs(pixcrd - pixcrd2)) < 1e-6
        flipim= numpy.flipud(image)
        double_image = numpy.float64(flipim)
        myRaster = arcpy.NumPyArrayToRaster(double_image)
        myRaster.save(outputfile_and_path)
        source_pnt = pixcrd_string
        #"'0 0';'0 2237';'2408 0';'2408 2237'"
        target_pnt = coordinate_string
        #"'202.92642053 47.01343105';'202.93045379 47.47945483';'202.19065519 47.01400053';'202.18819197 47.48003368'"
        arcpy.Warp_management(outputfile_and_path, source_pnt, target_pnt, geo_file_and_path, "POLYORDER1","BILINEAR")
        arcpy.Delete_management(outputfile_and_path)
        
arcpy.env.workspace = output_path
rasters = arcpy.ListRasters("*", "TIF")
print(rasters)
arcpy.CompositeBands_management(rasters,composite_filename)
print("Done.")