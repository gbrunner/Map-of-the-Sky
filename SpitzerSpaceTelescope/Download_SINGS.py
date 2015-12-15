##This script downloads a list of FITS files of nearby galaxies collected 
##under SINGS
##
##Created by Gregory Brunner (gregbrunn@gmail.com) on 7/17/2015

import requests
import os, sys

#galaxy_list = ['ngc0628', 'ngc2841', 'ngc2976', 'ngc3184', 'ngc3351', 'ngc3521', 'ngc3627', 'ngc3938', 'ngc4321', 'ngc4569', 'ngc4579', 'ngc4725', 'ngc4736', 'ngc4826', 'ngc6946', 'ngc7331', 'ngc0925', 'ngc5055']
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
            
write_to_folder = 'C:/PROJECTS/R&D/ASTROARC/SINGS/FITS/'
for galaxy in galaxy_list:

    instrument = "IRAC"
    numBands = 4

    for band in range(1, numBands+1):
        url = "http://irsa.ipac.caltech.edu/data/SPITZER/SINGS/galaxies/"+galaxy+"/"+instrument + "/" +galaxy+"_v7.phot."+str(band)+".fits"
        print(url)
        output_filename = galaxy+"_v7.phot."+str(band)+".fits"
        local_filename = write_to_folder + output_filename
        print(local_filename)
        
        if not os.path.exists(local_filename):
            r = requests.get(url)            
            with open(local_filename, 'wb') as f:
                for chunk in r.iter_content(chunk_size=1024): 
                    if chunk: # filter out keep-alive new chunks
                        f.write(chunk)
                        f.flush()
#return local_filename
print("Done.")