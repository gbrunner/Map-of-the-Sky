## This script converts SINGS galaxies from FITS files to TIFF files and also
## projects the TIFF files into a WCS so that the data can be consumed in a 
## GIS (ArcGIS).  It also uses Composite Bands to stack the 4 IRAC bands into
## a single image
##
##Created by Gregory Brunner (gregbrunn@gmail.com) on 7/17/2015

from __future__ import division, print_function

import numpy
from astropy import wcs
from astropy.io import fits
import os, sys
#import h5py
import arcpy

output_path = 'C:/PROJECTS/R&D/ASTROARC/SINGS/TIFF/'
input_path = 'C:/PROJECTS/R&D/ASTROARC/SINGS/FITS/'
numBands = 4

#galaxy_list = ['ngc0024','ngc0628', 'ngc2841', 'ngc2976', 'ngc3184', 'ngc3351', 'ngc3521', 'ngc3627', 'ngc3938', 'ngc4321', 'ngc4569', 'ngc4579', 'ngc4725', 'ngc4736', 'ngc4826', 'ngc6946', 'ngc7331', 'ngc0925', 'ngc5055']
#galaxy_list = ['ngc0024', 'ngc0337', 'ngc0584', 'ngc0628', 'ngc0855', 'ngc0925',
#            'ngc1097', 'ngc1266', 'ngc1291', 'ngc1316', 'ngc1377', 'ngc1404', 
#            'ngc1482', 'ngc1512', 'ngc1566', 'ngc1705', 'ngc2403', 'ngc2798',
#            'ngc2841', 'ngc2915', 'ngc2976', 'ngc3049', 'ngc3031', 'ngc3034',
#            'ngc3190', 'ngc3184', 'ngc3198', 'ngc3265', 'ngc3351', 'ngc3521',
#            'ngc3621', 'ngc3627', 'ngc3773', 'ngc3938', 'ngc4125', 'ngc4236',
#            'ngc4254', 'ngc4321', 'ngc4450', 'ngc4536', 'ngc4552', 'ngc4559',
#            'ngc4569', 'ngc4579', 'ngc4594', 'ngc4625', 'ngc4631', 'ngc4725',
#            'ngc4736', 'ngc4826', 'ngc5033', 'ngc5055', 'ngc5194', 'ngc5195',
#            'ngc5408', 'ngc5474', 'ngc5713', 'ngc5866', 'ngc6822',
#            'ngc6946', 'ngc7331', 'ngc7552', 'ngc7793']
            
galaxy_list = ['ic2574', 'ic4710', 'mrk33', 'hoii', 'm81dwa', 'ddo053', 'hoi', 'hoix', 'm81dwb', 
               'ddo154', 'ddo165', 'tol89']   
            
arcpy.env.pyramid = "NONE"
arcpy.env.rasterStatistics  = "NONE"
            
for galaxy in galaxy_list:
    #galaxy = "ngc5195"
    instrument = "IRAC"
    composite_filename = output_path + "COMPOSITES" + "\\" + galaxy + "_" + instrument + "_" + str(numBands) + "Bands.tif"
    if not os.path.exists(composite_filename):    
        for num in range(1,numBands+1):
            filename = galaxy + "_v7.phot." + str(num)+ ".fits"
            output_filename = filename[:-4] + 'rev.tif'
            geo_output_filename = filename[:-4] + 'rev.WCS.tif'
            outputfile_and_path = output_path + output_filename
            geo_file_and_path = output_path + geo_output_filename
            file_and_path = input_path + filename
            print(outputfile_and_path)
            print(geo_file_and_path)
            hdulist = fits.open(file_and_path)
            #print(hdulist.info())
            w = wcs.WCS(hdulist[0].header)
            image = hdulist[0].data

            # Print out the "name" of the WCS, as defined in the FITS header
            #print(w.wcs.name)

            # Print out all of the settings that were parsed from the header
            #w.wcs.print_contents()

            # Some pixel coordinates of interest.
            pixcrd = numpy.array([[0, 0], [0, image.shape[0]],[image.shape[1], 0], [image.shape[1], image.shape[0]]], numpy.float_)
            pixcrd_string = '0 0;' + '0 ' + str(image.shape[0]) + ';' +  str(image.shape[1]) + ' 0;' + str(image.shape[1]) + ' ' + str(image.shape[0])

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
            target_pnt = coordinate_string
            arcpy.Warp_management(outputfile_and_path, source_pnt, target_pnt, geo_file_and_path, "POLYORDER1","BILINEAR")
            arcpy.Delete_management(outputfile_and_path)
        
        arcpy.env.workspace = output_path
        rasters = arcpy.ListRasters(galaxy + "*", "TIF")
        print(rasters)
        arcpy.CompositeBands_management(rasters,composite_filename)
        
print("Done.")