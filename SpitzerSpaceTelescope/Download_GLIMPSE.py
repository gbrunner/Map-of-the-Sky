##Automate downloading the GLIMPSE archive
from __future__ import division, print_function


from bs4 import BeautifulSoup
#import http
import urllib.request
import os #, sys
import requests


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

class GLIMPSE(object):
   
   def __init__(self, write_to_folder, webpage, output_TIFF_path, projection):
       self.write_to_folder = write_to_folder #'C:/PROJECTS/R&D/ASTROARC/GLIMPSE/II/0.6_mosaics_v2.0/'
       #webpage = "http://irsa.ipac.caltech.edu/data/SPITZER/GLIMPSE/images/I/0.6_mosaics_v2.0"
       self.webpage = webpage #"http://irsa.ipac.caltech.edu/data/SPITZER/GLIMPSE/images/II/0.6_mosaics_v2.0"
       self.output_TIFF_path = output_TIFF_path
       self.projection = projection
       self.filelist = []
       self.linklist = []

   def getFileAndLinkLists(self):
       print('Creating lists of links and output files...')

       with urllib.request.urlopen(self.webpage) as url:
           s = url.read()
       
       soup = BeautifulSoup(s, 'html.parser')

    
       for link in soup.find_all('a'):
           first_path = link.get('href')
           with urllib.request.urlopen(self.webpage + '/' + first_path) as second_url:
               ds = second_url.read()
           deep_soup = BeautifulSoup(ds, 'html.parser')
           for file_link in deep_soup.find_all('a'):      
               temp_link = file_link.get('href')        
               if (temp_link[-5:] == '.fits'):
                   self.filelist.append(file_link.get('href'))
                   self.linklist.append(self.webpage + '/' + first_path + file_link.get('href'))


   def DownloadGLIMPSE(self):
     
       self.getFileAndLinkLists()       
       
       print('Starting Download Process...')

       if not os.path.exists(self.write_to_folder):
           os.makedirs(self.write_to_folder)
                   
       num_files = str(len(self.filelist))
       print(num_files)
       counter = 0
       for file in self.filelist:
           output_filename = self.write_to_folder + '/' + file
           if not os.path.exists(output_filename):
               print('Downloading file ' + str(counter+1) + ' of ' + num_files + ', ' + output_filename + ' from ' + self.linklist[counter])
               r = requests.get(self.linklist[counter]) #(webpage+link)
               with open(output_filename, 'wb') as f:
                   for chunk in r.iter_content(chunk_size=1024): 
                       if chunk: # filter out keep-alive new chunks
                           f.write(chunk)
                           f.flush()
           else:
               print(output_filename + ' already exists. Skipping...')
           counter = counter+1
       print("Done Downloading GLIMPSE Data.")
       
   def ConvertToTIFF(self):
      arcpy.env.pyramid = "NONE"
      arcpy.env.rasterStatistics  = "NONE"       
       
      if not os.path.exists(self.output_TIFF_path):
          os.makedirs(self.output_TIFF_path)
          #Get list of files that will be converted
      fits_files = [ f for f in listdir(self.write_to_folder) if isfile(join(self.write_to_folder,f)) ]

      counter = 1
      #For every file...
      for fits_file in fits_files:    
          print("Processing file: " + str(counter) + ", " + fits_file + "...")
          counter = counter + 1
          filename = fits_file
          output_filename = filename[:-4] + 'rev.tif'
          geo_output_filename = filename[:-4] + 'rev.' + self.projection + '.tif'
          outputfile_and_path = self.output_TIFF_path + "\\" + output_filename
          geo_file_and_path = self.output_TIFF_path + "\\" + geo_output_filename
          file_and_path = self.write_to_folder + "\\" + filename
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
              if (self.projection == 'WCS'):       
                  coordinate_string = str(-1*world[0,0]) + ' ' + str(world[0,1]) + ';' + str(-1*world[1,0]) + ' ' + str(world[1,1]) + ';' + str(-1*world[2,0]) + ' ' + str(world[2,1]) + ';' + str(-1*world[3,0]) + ' ' + str(world[3,1]) 
                  #print(coordinate_string)
              if (self.projection == 'GALACTIC'):
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
      
   def CreateMultibandImages(self):
       arcpy.env.workspace = self.output_TIFF_path
       rasList = arcpy.ListRasters("*I1*")
       instrument = "IRAC"
       numBands = 4
       os.mkdir(self.output_TIFF_path + "\\" + "COMPOSITES" + "\\")
       for ras in rasList:
           composite_filename = self.output_TIFF_path + "\\" + "COMPOSITES" + "\\" + ras[:22] + instrument + "_" + str(numBands) + "Bands.tif"        
           if not os.path.exists(composite_filename):    
               rasters = ras[:22] + 'I1.rev.GALACTIC.tif;' + ras[:22] + 'I2.rev.GALACTIC.tif;' + ras[:22] + 'I3.rev.GALACTIC.tif;' + ras[:22] + 'I4.rev.GALACTIC.tif'
               print(rasters)
               arcpy.CompositeBands_management(rasters,composite_filename)
               
   def CreateMosaicDataset(self):
       arcpy.CreateMosaicDataset_management(workspace, mosaic_name, coordinate_sys, 
                                     num_bands="", pixel_type="", product_definition="NONE", 
                                     product_band_definitions="")
       arcpy.AddRastersToMosaicDataset_management(mosaic_dataset, "Raster Dataset", 
                                           path_to_files, "UPDATE_CELL_SIZES", 
                                           "UPDATE_BOUNDARY", "NO_OVERVIEWS", 
                                           "", "0", 
                                           "1500", imagery_spatial_ref, 
                                           "#", "SUBFOLDERS", 
                                           "ALLOW_DUPLICATES", 
                                           "NO_PYRAMIDS", 
                                           "NO_STATISTICS", 
                                           "NO_THUMBNAILS", "#", 
                                           "NO_FORCE_SPATIAL_REFERENCE")
#        arcpy.DefineMosaicDatasetNoData_management(mosaic_dataset, num_bands="4", 
#                                           bands_for_nodata_value="ALL_BANDS 0", 
#                                           bands_for_valid_data_range="", where_clause="", 
#                                           Composite_nodata_value="NO_COMPOSITE_NODATA")
#        arcpy.BuildFootprints_management(mosaic_dataset, where_clause="", 
#                                 reset_footprint="RADIOMETRY", min_data_value="1", 
#                                 max_data_value="1.84467440737096E+19", 
#                                 approx_num_vertices="80", shrink_distance="0", 
#                                 maintain_edges="NO_MAINTAIN_EDGES", 
#                                 skip_derived_images="SKIP_DERIVED_IMAGES", 
#                                 update_boundary="UPDATE_BOUNDARY", request_size="2000", 
#                                 min_region_size="100", simplification_method="NONE", 
#                                 edge_tolerance="", max_sliver_size="20", 
#                                 min_thinness_ratio="0.05")
#        arcpy.BuildOverviews_management(mosaic_dataset, where_clause="", 
#                                define_missing_tiles="DEFINE_MISSING_TILES", 
#                                generate_overviews="GENERATE_OVERVIEWS", 
#                                generate_missing_images="GENERATE_MISSING_IMAGES", 
#                                regenerate_stale_images="REGENERATE_STALE_IMAGES")

write_to_folder = 'C:/PROJECTS/R&D/ASTROARC/GLIMPSE/I/0.6_mosaics_v2.0'
webpage = "http://irsa.ipac.caltech.edu/data/SPITZER/GLIMPSE/images/I/0.6_mosaics_v2.0"
output_TIFF_Path = 'C:/PROJECTS/R&D/ASTROARC/GLIMPSE/I/0.6_mosaics_v2.0/TIFF/GALACTIC'
projection = 'GALACTIC'
runGLIMPSE = GLIMPSE(write_to_folder, webpage, output_TIFF_Path, projection)
runGLIMPSE.DownloadGLIMPSE()
#runGLIMPSE.ConvertToTIFF()
#runGLIMPSE.CreateMultibandImages()