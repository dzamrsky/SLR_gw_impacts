# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 15:25:12 2022

@author: daniel
"""

import tarfile
import os
import shutil

input_dir = r'/projects/0/qt16165/_dzamrsky/_A4_SRM_models/_SLR_models'
main_dir = r'/projects/0/qt16165/_dzamrsky/_A4_SRM_models/_conc_plots_2000AD'
temp_dir = r'/projects/0/qt16165/_dzamrsky/_A4_SRM_models/temp'

# look for the distributions at 30000 and copy those to to main directory
for path, subdirs, files in os.walk(input_dir):
    for name in files:
        if '_img_to_vid' in name and 'avg_sim' in path:
            print(os.path.join(path, name))
            
            tr_file_dir = os.path.join(path, name)
            
            # open the tarfile and extract contents
            tar = tarfile.open(tr_file_dir)
            tar.extractall(path = temp_dir)
            tar.close()
            
            #   get the SRM and COSCAT number
            srm_name = [i for i in path.split('/') if '_summary' in i][0].split('_summary_')[1]

            # look for the distributions at 30000 and copy those to to main directory
            for path, subdirs, files in os.walk(temp_dir):
                for name in files:
                    if 'RCP_45_img_30000' in name:
                        print(os.path.join(path, name))
                        shutil.copyfile(os.path.join(path, name), os.path.join(main_dir, srm_name + '.png'))

            # look for the distributions at 30000 and copy those to to main directory
            for path, subdirs, files in os.walk(temp_dir):
                for name in files:
                    if 'RCP_45_img_500' in name:
                        print(os.path.join(path, name))
                        shutil.copyfile(os.path.join(path, name), os.path.join(main_dir, srm_name + '_IC.png'))

            shutil.rmtree(temp_dir)