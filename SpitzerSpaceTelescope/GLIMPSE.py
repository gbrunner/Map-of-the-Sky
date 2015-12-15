##Automate downloading the GLIMPSE archive
#from __future__ import division, print_function
#Native imports
import os
from os import listdir
from os.path import isfile, join
import urllib.request
#TODO: import multiprocessing

#Anaconda imports
import numpy
from bs4 import BeautifulSoup
import requests
from astropy import wcs
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord

#ArcGIS imports
import arcpy

class GLIMPSE(object):

    def __init__(self, write_to_folder, webpage, output_TIFF_path, projection, output_gdb, output_mosaic):
        self.write_to_folder = write_to_folder
        self.webpage = webpage
        self.output_TIFF_path = output_TIFF_path
        self.composite_path = join(self.output_TIFF_path, "COMPOSITES")
        self.projection = projection
        self.output_gdb = output_gdb
        self.output_mosaic = output_mosaic
        self.numBands = 4
        self.instrument = "IRAC"
        self.mosaic_spatial_ref = "PROJCS['WGS_1984_Web_Mercator_Auxiliary_Sphere',GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Mercator_Auxiliary_Sphere'],PARAMETER['False_Easting',0.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',0.0],PARAMETER['Standard_Parallel_1',0.0],PARAMETER['Auxiliary_Sphere_Type',0.0],UNIT['Meter',1.0]];-20037700 -30241100 10000;-100000 10000;-100000 10000;0.001;0.001;0.001;IsHighPrecision"
        self.imagery_spatial_ref = "GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]];-400 -400 1000000000;-100000 10000;-100000 10000;8.98315284119522E-09;0.001;0.001;IsHighPrecision"
        self.filelist = []
        self.linklist = []

    #Creates a list of fits files and their download locations
    def get_file_and_link_lists(self):
        print('Creating lists of links and output files...')

        with urllib.request.urlopen(self.webpage) as url:
            s = url.read()

        soup = BeautifulSoup(s, 'html.parser')

        for link in soup.find_all('a'):
            first_path = link.get('href')
            with urllib.request.urlopen(join(self.webpage, first_path)) as second_url:
                ds = second_url.read()
            deep_soup = BeautifulSoup(ds, 'html.parser')
            for file_link in deep_soup.find_all('a'):
                temp_link = file_link.get('href')
                if (temp_link[-5:] == '.fits'):
                    self.filelist.append(file_link.get('href'))
                    self.linklist.append(join(self.webpage, first_path + file_link.get('href')))

    #Download the GLIMPSE archive
    def download_glimpse(self):
        self.get_file_and_link_lists()

        print('Starting Download Process...')

        if not os.path.exists(self.write_to_folder):
            os.makedirs(self.write_to_folder)

        num_files = str(len(self.filelist))
        print(num_files)
        for file in enumerate(self.filelist):
            output_filename = join(self.write_to_folder,file[1])
            if not os.path.exists(output_filename):
                print('Downloading file ' + str(file[0]+1) + ' of ' + num_files + ', ' + output_filename + ' from ' + self.linklist[file[0]])
                r = requests.get(self.linklist[file[0]]) #(webpage+link)
                with open(output_filename, 'wb') as f:
                    for chunk in r.iter_content(chunk_size=1024):
                        if chunk: # filter out keep-alive new chunks
                            f.write(chunk)
                            f.flush()
            else:
                print(output_filename + ' already exists. Skipping...')
        print("Done Downloading GLIMPSE Data.")

    #Images need to be 'shifted' to -180 <= longitude <= 180 to work right in
    #ArcGIS and ArcGIS.com
    def normalize_x(self, val):
        if val > 180:
            return val-360
        elif val <-180:
            return val+360
        else:
            return val

    #TODO: write a function that forms the string that gets passed to
    #Warp_management()
    def create_warp_string():
        pass

    #converts fits files to tiff files
    def convert_to_tiff(self):
        arcpy.env.pyramid = "NONE"
        arcpy.env.rasterStatistics  = "NONE"

        if not os.path.exists(self.output_TIFF_path):
            os.makedirs(self.output_TIFF_path)
            #Get list of files that will be converted
        fits_files = [ f for f in listdir(self.write_to_folder) if isfile(join(self.write_to_folder,f)) ]

        #For every file...
        for fits_file in enumerate(fits_files):
            print("Processing file: " + str(fits_file[0]+1) + ", " + fits_file[1] + "...")
            filename = fits_file[1]
            output_filename = filename[:-4] + 'rev.tif'
            geo_output_filename = filename[:-4] + 'rev.' + self.projection + '.tif'
            outputfile_and_path = join(self.output_TIFF_path, output_filename)
            geo_file_and_path = join(self.output_TIFF_path, geo_output_filename)
            file_and_path = join(self.write_to_folder,filename)
            if not os.path.exists(geo_file_and_path):
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
                #GLIMPSE is Backwards! The data is originially in GALACTIC so
                #There is no need to change the projection
                longitude = []
                latitude = []
                if (self.projection == 'GALACTIC'):
                    coord_string = str(-1*world[0,0]) + ' ' + str(world[0,1]) + ';' + str(-1*world[1,0]) + ' ' + str(world[1,1]) + ';' + str(-1*world[2,0]) + ' ' + str(world[2,1]) + ';' + str(-1*world[3,0]) + ' ' + str(world[3,1])
                    #coord_string = str(world[0,0]) + ' ' + str(world[0,1]) + ';' + str(world[1,0]) + ' ' + str(world[1,1]) + ';' + str(world[2,0]) + ' ' + str(world[2,1]) + ';' + str(world[3,0]) + ' ' + str(world[3,1])
                    print(coord_string)
                    arg_space = coord_string.replace(';', ' ')
                    args = arg_space.split(' ')
                    #Fixing the Coordinates
                    longitude.append(self.normalize_x(float(args[0])))
                    longitude.append(self.normalize_x(float(args[2])))
                    longitude.append(self.normalize_x(float(args[4])))
                    longitude.append(self.normalize_x(float(args[6])))
                    latitude.append(float(args[1]))
                    latitude.append(float(args[3]))
                    latitude.append(float(args[5]))
                    latitude.append(float(args[7]))
                if (self.projection == 'WCS'):
                    # Convert the same coordinates back to pixel coordinates.
                    #pixcrd2 = w.wcs_world2pix(world, 1)
                    coordinate = SkyCoord(world,frame='icrs', unit='deg')
                    coord = coordinate.galactic
                    print(coord.to_string('decimal'))
                    coord_string =   coord.to_string('decimal')
                    longitude.append(self.normalize_x(float(coord_string[0].split()[0])))
                    longitude.append(self.normalize_x(float(coord_string[1].split()[0])))
                    longitude.append(self.normalize_x(float(coord_string[2].split()[0])))
                    longitude.append(self.normalize_x(float(coord_string[3].split()[0])))
                    latitude.append(float(coord_string[0].split()[1]))
                    latitude.append(float(coord_string[1].split()[1]))
                    latitude.append(float(coord_string[2].split()[1]))
                    latitude.append(float(coord_string[3].split()[1]))

                coordinate_string = str(longitude[0]) + ' ' + str(latitude[0]) + ';' + str(longitude[1]) + ' ' + str(latitude[1]) + ';' + str(longitude[2]) + ' ' + str(latitude[2]) + ';' + str(longitude[3]) + ' ' + str(latitude[3])
                print(coordinate_string)
                # These should be the same as the original pixel coordinates, modulo
                # some floating-point error.
                #assert numpy.max(numpy.abs(pixcrd - pixcrd2)) < 1e-6
                flipim = numpy.flipud(image)
                double_image = numpy.float64(flipim)
                myRaster = arcpy.NumPyArrayToRaster(double_image)
                myRaster.save(outputfile_and_path)
                source_pnt = pixcrd_string
                target_pnt = coordinate_string
                arcpy.Warp_management(outputfile_and_path, source_pnt, target_pnt, geo_file_and_path, "POLYORDER1","BILINEAR")
                arcpy.Delete_management(outputfile_and_path)

    #Uses composite bands to create 4 band images of the galactic plane
    def create_multiband_images(self):
        arcpy.env.pyramid = "PYRAMIDS"
        arcpy.env.rasterStatistics  = "STATISTICS"

        arcpy.env.workspace = self.output_TIFF_path
        rasList = arcpy.ListRasters("*I1*")
        if not os.path.exists(self.composite_path):
            os.mkdir(self.composite_path)
        for ras in rasList:
            composite_filename = join(self.composite_path,  ras[:22] + self.instrument + "_" + str(self.numBands) + "Bands.tif")
            if not os.path.exists(composite_filename):
                rasters = ras[:22] + 'I1.rev.' + self.projection+'.tif;' + ras[:22] + 'I2.rev.'+ self.projection+'.tif;' + ras[:22] + 'I3.rev.'+ self.projection+'.tif;' + ras[:22] + 'I4.rev.'+ self.projection+'.tif'
                print('Creating Composite Band Raster from: ' + rasters)
                arcpy.CompositeBands_management(rasters,composite_filename)
                for n in range(self.numBands):
                    print('Deleting ' + ras[:22] + 'I' + str(n+1) +'.rev.' + self.projection+'.tif')
                    arcpy.Delete_management(ras[:22] + 'I' + str(n+1) +'.rev.' + self.projection+'.tif')

    #Creates the ArcGIS mosaic dataset for GLIMPSE legacy program data
    def create_mosaic_dataset(self):
        if not arcpy.Exists(self.output_gdb):
            arcpy.CreateFileGDB_management(os.path.dirnam(output_gdb),os.path.basename(output_gdb))

        if not arcpy.Exists(join(self.output_gdb, self.output_mosaic)):
            arcpy.CreateMosaicDataset_management(self.output_gdb, self.output_mosaic,
                                    self.mosaic_spatial_ref, num_bands=self.numBands,
                                    pixel_type="", product_definition="NONE",
                                    product_band_definitions="")

        arcpy.AddRastersToMosaicDataset_management(join(self.output_gdb, self.output_mosaic),
                                    "Raster Dataset", self.composite_path,
                                    "UPDATE_CELL_SIZES", "UPDATE_BOUNDARY",
                                    "UPDATE_OVERVIEWS", "", "0",
                                    "1500", self.imagery_spatial_ref,
                                    "#", "SUBFOLDERS",
                                    "ALLOW_DUPLICATES",
                                    "NO_PYRAMIDS",
                                    "NO_STATISTICS",
                                    "NO_THUMBNAILS", "#",
                                    "NO_FORCE_SPATIAL_REFERENCE")

write_to_folder = 'C:/PROJECTS/R&D/ASTROARC/GLIMPSE/I/0.6_mosaics_v2.0'
webpage = "http://irsa.ipac.caltech.edu/data/SPITZER/GLIMPSE/images/I/0.6_mosaics_v2.0"
output_TIFF_Path = 'C:/PROJECTS/R&D/ASTROARC/GLIMPSE/I/0.6_mosaics_v2.0/TIFF/GALACTIC'
projection = 'GALACTIC'
output_gdb = "C:/PROJECTS/R&D/ASTROARC/GLIMPSE.gdb"
output_mosaic = "GLIMPSE_I"
runGLIMPSE = GLIMPSE(write_to_folder, webpage, output_TIFF_Path, projection, output_gdb, output_mosaic)
#runGLIMPSE.download_glimpse()
#runGLIMPSE.convert_to_tiff()
#runGLIMPSE.create_multiband_images()
runGLIMPSE.create_mosaic_dataset()
