##Build SINGS Mosaic Dataset

from __future__ import division, print_function

import numpy
from astropy import wcs
from astropy.io import fits
import requests
import arcpy
import os, sys

class SINGS(object):
   
   def __init__(self, input_fits_path, output_tiff_path, workspace, mosaic_name, mosaic_dataset, 
                path_to_files, coordinate_sys, imagery_spatial_ref):
       arcpy.env.pyramid = "NONE"
       arcpy.env.rasterStatistics  = "NONE"
       self.input_fits_path = input_fits_path
       self.output_tiff_path = output_tiff_path

       self.workspace = "C:/PROJECTS/R&D/ASTROARC/SINGS/SINGS.gdb"
       self.mosaic_name = "SINGS_4Band_IRAC_ALL_WMAS"
       self.path_to_files = "C:/PROJECTS/R&D/ASTROARC/SINGS/TIFF"
       self.mosaic_dataset = workspace + "/" + mosaic_name
       self.coordinate_sys = "PROJCS['WGS_1984_Web_Mercator_Auxiliary_Sphere',GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Mercator_Auxiliary_Sphere'],PARAMETER['False_Easting',0.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',0.0],PARAMETER['Standard_Parallel_1',0.0],PARAMETER['Auxiliary_Sphere_Type',0.0],UNIT['Meter',1.0]];-20037700 -30241100 10000;-100000 10000;-100000 10000;0.001;0.001;0.001;IsHighPrecision"
       self.imagery_spatial_ref = "GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]];-400 -400 1000000000;-100000 10000;-100000 10000;8.98315284119522E-09;0.001;0.001;IsHighPrecision"
       self.write_to_folder = input_fits_path
       self.instrument = "IRAC"
       self.numBands = 4
       self.base_url = "http://irsa.ipac.caltech.edu/data/SPITZER/SINGS/galaxies/"
       self.galaxy_list = galaxy_list

   def DownloadSINGS(self):
       for galaxy in self.galaxy_list:
           for band in range(1, self.numBands+1):
               url = self.base_url + galaxy + "/" + self.instrument + "/" + galaxy + "_v7.phot." + str(band) + ".fits"
               print(url)
               output_filename = galaxy+"_v7.phot."+str(band)+".fits"
               local_filename = self.input_fits_path + output_filename
               print(local_filename)               
               if not os.path.exists(local_filename):
                   r = requests.get(url)            
                   with open(local_filename, 'wb') as f:
                       for chunk in r.iter_content(chunk_size=1024): 
                           if chunk: # filter out keep-alive new chunks
                               f.write(chunk)
                               f.flush()
        #return local_filename
       print("Done Downloading SINGS FITS Files.")

   def ConvertFITStoTIFF(self):
       for galaxy in self.galaxy_list:
           composite_filename = self.path_to_files + '\\' + "COMPOSITES" + "\\" + galaxy + "_" + self.instrument + "_" + str(self.numBands) + "Bands.tif"
           if not os.path.exists(composite_filename):    
               for num in range(1,self.numBands+1):
                   filename            = galaxy + "_v7.phot." + str(num)+ ".fits"
                   output_filename     = filename[:-4] + 'rev.tif'
                   geo_output_filename = filename[:-4] + 'rev.WCS.tif'
                   outputfile_and_path = self.output_tiff_path + output_filename
                   geo_file_and_path   = self.output_tiff_path + geo_output_filename
                   file_and_path       = self.input_fits_path + filename
                   print(outputfile_and_path)
                   print(geo_file_and_path)
                   hdulist = fits.open(file_and_path)
                   #print(hdulist.info())
                   w = wcs.WCS(hdulist[0].header)
                   image = hdulist[0].data
                   
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
       print("Done Converting FITS to TIFF.")
       
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
                                           
                                           
output_path = 'C:/PROJECTS/R&D/ASTROARC/SINGS/TIFF/'
input_path = 'C:/PROJECTS/R&D/ASTROARC/SINGS/FITS/'
numBands = 4
galaxy_list = ['ngc0024', 'ngc0337', 'ngc0584', 'ngc0628', 'ngc0855', 'ngc0925',
            'ngc1097', 'ngc1266', 'ngc1291', 'ngc1316', 'ngc1377', 'ngc1404', 
            'ngc1482', 'ngc1512', 'ngc1566', 'ngc1705', 'ngc2403', 'ngc2798',
            'ngc2841', 'ngc2915', 'ngc2976', 'ngc3049', 'ngc3031', 'ngc3034',
            'ngc3190', 'ngc3184', 'ngc3198', 'ngc3265', 'ngc3351', 'ngc3521',
            'ngc3621', 'ngc3627', 'ngc3773', 'ngc3938', 'ngc4125', 'ngc4236',
            'ngc4254', 'ngc4321', 'ngc4450', 'ngc4536', 'ngc4552', 'ngc4559',
            'ngc4569', 'ngc4579', 'ngc4594', 'ngc4625', 'ngc4631', 'ngc4725',
            'ngc4736', 'ngc4826', 'ngc5033', 'ngc5055', 'ngc5194', 'ngc5195',
            'ngc5408', 'ngc5474', 'ngc5713', 'ngc5866', 'ngc6822', 'ngc6946', 
            'ngc7331', 'ngc7552', 'ngc7793','ic2574', 'ic4710', 'mrk33', 'hoii', 
            'm81dwa', 'ddo053', 'hoi', 'hoix', 'm81dwb', 'ddo154', 'ddo165', 'tol89']   