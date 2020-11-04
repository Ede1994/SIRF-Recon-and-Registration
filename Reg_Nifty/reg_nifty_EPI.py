# Copyright University College London Hospital, 2020
# Author: Eric Einspaenner, Institute of Nuclear Medicine
# For internal research only.

#%% import all necessary modules
import os
import shutil
import re
import numpy
import math
from art import tprint
import sirf.Reg as Reg
import sirf.Reg as Eng_ref
import sirf.Reg as Eng_flo


#%% functions
# function for sorting lists, containing strings with numbers (https://stackoverflow.com/a/48413757)
def sorted_alphanumeric(data):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(data, key=alphanum_key, reverse=False)

# registration of EPI images
def reg_epi(ref, flo_path):
    i = 0
    epi_reg = Reg.NiftyAladinSym()
    epi_reg.set_reference_image(ref)
    for image in sorted_alphanumeric(os.listdir(flo_path)):
        flo_file = flo_path + image
        print(flo_file)
        flo = Eng_flo.ImageData(flo_file)
        epi_reg.set_floating_image(flo)
        epi_reg.process()
        tm_epi = epi_reg.get_transformation_matrix_forward()
        tm_epi.write(path_EPI + 'tm_epi_' + str(i))
        i += 1


#%% data path and set filenames

py_path = os.getcwd()

# avoid cluttering of files, delete working-folder and start from scratch
working_folder = py_path + '/working'

# change the current working directory to the given path
os.chdir(working_folder)


#%% create folders for results

#%% Create folders for results

path_EPI = working_folder + '/EPI/'

if not os.path.exists(path_EPI):
    os.makedirs(path_EPI, mode=0o770)
    print('Create Folder: {}'.format(path_EPI))


#%% Registration EPI, delivers transformation matrices 
# define reference image (first image) and float-path

tprint('Start Registration of EPIs')

# refernce file
epi_data_path = '/home/eric/Dokumente/PythonProjects/SIRF/UKL_data/EPI/20191009/1/'
ref_file = epi_data_path + 'epi_0.nii'
ref = Eng_ref.ImageData(ref_file)

# float files
flo_path = epi_data_path

# Niftyreg with EPI images
reg_epi(ref, flo_path)

tprint('Finish Registration')


translation_values = []
rotation_values = []

#%% get TM from Reg and calculate mean translation amplitude for every EPI volume
for file in sorted_alphanumeric(os.listdir(path_EPI)):
    print(path_EPI + file)
    tm_fwD_arr = numpy.loadtxt(path_EPI + file)
    translation = numpy.sqrt((tm_fwD_arr[0][3]) ** 2 + (tm_fwD_arr[1][3]) ** 2 + (tm_fwD_arr[2][3]) ** 2)
    translation_values.append(translation)
    rotation = numpy.sqrt((math.atan2(tm_fwD_arr[2][1], tm_fwD_arr[2][2])) ** 2 + (math.atan2(-tm_fwD_arr[2][0], numpy.sqrt(tm_fwD_arr[2][1] ** 2 + tm_fwD_arr[2][2] ** 2))) ** 2 + (math.atan2(tm_fwD_arr[1][0], tm_fwD_arr[0][0])) ** 2)
    rotation_values.append(rotation)


#%% save list of translation amplitude values as txt-file
with open('translation.txt', 'w') as f:
    for item in translation_values:
        f.write("%s\n" % item)

with open('rotation.txt', 'w') as f:
    for item in rotation_values:
        f.write("%s\n" % item)

tprint('DONE!')