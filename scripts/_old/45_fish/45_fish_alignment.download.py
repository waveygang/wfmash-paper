# Prerequisites on Octopus
# guix package -A write_something
# guix install python-pandas
# guix install python-urllib3

import sys

path_45_fish_fixed = sys.argv[1]
dir_output = sys.argv[2]


# Load information
import pandas as pd

fish_df = pd.read_excel(open(path_45_fish_fixed, 'rb'))


# Construct URLs and download the genomes

# Directory and ID structure
# https://ftp.ncbi.nlm.nih.gov/genomes/all/README.txt
# URL example: ftp.ncbi.nlm.nih.gov/genomes/all/GCA/904/848/185/GCA_904848185.1_fAcaLat1.1/GCA_904848185.1_fAcaLat1.1_genomic.fna.gz

import os
import urllib.request

for index, row in fish_df.iterrows():
    #print(row['ToL_id'], row['Accession'])
    
    accession_ToLId = row['Accession'] + '_' + row['ToL_id']
    
    prefix, levels_and_version = row['Accession'].split('_')
    level1 = levels_and_version[:3]
    level2 = levels_and_version[3:6]
    level3 = levels_and_version[6:9]
    assembly_version = levels_and_version[9:]

    ftp_url = f'https://ftp.ncbi.nlm.nih.gov/genomes/all/{prefix}/{level1}/{level2}/{level3}/{accession_ToLId}/{accession_ToLId}_genomic.fna.gz'
    
    print(ftp_url)

    urllib.request.urlretrieve(ftp_url, os.path.join(dir_output, f'{accession_ToLId}_genomic.fna.gz'))
