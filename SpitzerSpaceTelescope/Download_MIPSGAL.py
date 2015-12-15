
from bs4 import BeautifulSoup
import http
import urllib.request
import os, sys
import requests

#print(s)
filelist = []
write_to_folder = 'C:/PROJECTS/R&D/ASTROARC/MIPSGAL'
webpage = "http://irsa.ipac.caltech.edu/data/SPITZER/MIPSGAL/images/mosaics24/"

with urllib.request.urlopen(webpage) as url:
    s = url.read()

soup = BeautifulSoup(s, 'html.parser')
    
for link in soup.find_all('a'):
    if (len(link.get('href')) == 19):
        filelist.append(link.get('href'))
        #print(webpage + link.get('href'))

print(str(len(filelist))+ " Files to Download")


for link in filelist:
    output_filename = write_to_folder + '/' + link
    print(output_filename)
    print(webpage+link)
    if not os.path.exists(output_filename):
        r = requests.get(webpage+link)
        with open(output_filename, 'wb') as f:
            for chunk in r.iter_content(chunk_size=1024): 
                if chunk: # filter out keep-alive new chunks
                    f.write(chunk)
                    f.flush()

print("Done.")        
