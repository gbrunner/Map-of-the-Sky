from __future__ import division, print_function

#imports
import numpy
from astropy import wcs
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
import arcpy
import os
from os import listdir
from os.path import isfile, join

arcpy.env.pyramid = "NONE"
arcpy.env.rasterStatistics  = "NONE"

#Define where files are and where they are going
output_path = "\\\\STLGREG\\ASTROARC\\MIPSGAL\\FITS_GALACTIC_REV"
input_path = "C:\\PROJECTS\\R&D\\ASTROARC\\MIPSGAL\\FITS"
#instrument = "MIPS"
projection = 'GALACTIC'

if not os.path.exists(output_path):
    os.makedirs(output_path)

#Get list of files that will be converted
fits_files = [ f for f in listdir(input_path) if isfile(join(input_path,f)) ]

counter = 1
#For every file...
for fits_file in fits_files:
    
    print("Processing file: " + str(counter) + ", " + fits_file + "...")
    counter = counter + 1

    filename = fits_file
    output_filename = filename[:-4] + 'rev.tif'
    geo_output_filename = filename[:-4] + 'rev.' + projection + '.tif'
    outputfile_and_path = output_path + "\\" + output_filename
    geo_file_and_path = output_path + "\\" + geo_output_filename
    file_and_path = input_path + "\\" + filename

    if not os.path.exists(geo_file_and_path):

        #Print statements just to be sure everything is right    
        #print(outputfile_and_path)
        #print(geo_file_and_path)

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
        #print(world)        
        if (projection == 'WCS'):       
            coordinate_string = str(-1*world[0,0]) + ' ' + str(world[0,1]) + ';' + str(-1*world[1,0]) + ' ' + str(world[1,1]) + ';' + str(-1*world[2,0]) + ' ' + str(world[2,1]) + ';' + str(-1*world[3,0]) + ' ' + str(world[3,1]) 
            #print(coordinate_string)
        if (projection == 'GALACTIC'):
        # Convert the same coordinates back to pixel coordinates.
            pixcrd2 = w.wcs_world2pix(world, 1)
            coordinate = SkyCoord(world,frame='icrs', unit='deg')
            coord = coordinate.galactic
            print(coord.to_string('decimal'))           
            coord_string =   coord.to_string('decimal')          

            longitude = []
            latitude = []
            #print(coord_string[0].split()[1])
            longitude.append(float(coord_string[0].split()[0]))
            longitude.append(float(coord_string[1].split()[0]))
            longitude.append(float(coord_string[2].split()[0]))
            longitude.append(float(coord_string[3].split()[0]))      
            latitude.append(float(coord_string[0].split()[1]))
            latitude.append(float(coord_string[1].split()[1]))
            latitude.append(float(coord_string[2].split()[1]))
            latitude.append(float(coord_string[3].split()[1]))
            if longitude[0] >= 1:
                coordinate_string = str(longitude[0]) + ' ' + str(latitude[0]) + ';' + str(longitude[1]) + ' ' + str(latitude[1]) + ';' + str(longitude[2]) + ' ' + str(latitude[2]) + ';' + str(longitude[3]) + ' ' + str(latitude[3])
            else:
                longitude[2] = longitude[2]-360
                longitude[3] = longitude[3]-360
                coordinate_string = str(longitude[0]) + ' ' + str(latitude[0]) + ';' + str(longitude[1]) + ' ' + str(latitude[1]) + ';' + str((longitude[2])) + ' ' + str(latitude[2]) + ';' + str((longitude[3])) + ' ' + str(latitude[3])
            #print(coordinate_string)
            
            if longitude[0] >= 180:
                longitude[0] = (longitude[0]-360)
                longitude[1] = (longitude[1]-360)              
                longitude[2] = (longitude[2]-360)
                longitude[3] = (longitude[3]-360)
                #coordinate_string = str(longitude[0]) + ' ' + str(latitude[0]) + ';' + str(longitude[1]) + ' ' + str(latitude[1]) + ';' + str(longitude[2]) + ' ' + str(latitude[2]) + ';' + str(longitude[3]) + ' ' + str(latitude[3])
            
            longitude[0] = -1*(longitude[0])
            longitude[1] = -1*(longitude[1])              
            longitude[2] = -1*(longitude[2])
            longitude[3] = -1*(longitude[3])
            coordinate_string = str(longitude[0]) + ' ' + str(latitude[0]) + ';' + str(longitude[1]) + ' ' + str(latitude[1]) + ';' + str(longitude[2]) + ' ' + str(latitude[2]) + ';' + str(longitude[3]) + ' ' + str(latitude[3])
           
            print(coordinate_string)

        # These should be the same as the original pixel coordinates, modulo
        # some floating-point error.
        assert numpy.max(numpy.abs(pixcrd - pixcrd2)) < 1e-6
        flipim = numpy.flipud(image)
        double_image = numpy.float64(flipim)
        myRaster = arcpy.NumPyArrayToRaster(double_image)
        myRaster.save(outputfile_and_path)
        source_pnt = pixcrd_string
        #"'0 0';'0 2237';'2408 0';'2408 2237'"
        target_pnt = coordinate_string

        arcpy.Warp_management(outputfile_and_path, source_pnt, target_pnt, geo_file_and_path, "POLYORDER1","BILINEAR")
        arcpy.Delete_management(outputfile_and_path)
#        
#
print("Done.")