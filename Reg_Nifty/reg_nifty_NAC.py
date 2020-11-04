# Copyright University College London Hospital, 2020
# Author: Eric Einspaenner, Institute of Nuclear Medicine
# For internal research only.

#%% import all necessary modules
import os
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

# registration of NAC images
def reg_nac(ref, flo_path, list_NACs):
    i = 0
    nac_reg = Reg.NiftyAladinSym()
    nac_reg.set_reference_image(ref)
    for image in sorted_alphanumeric(list_NACs):
        nac_file = flo_path + image
        print(nac_file)
        flo = Eng_flo.ImageData(nac_file)
        nac_reg.set_floating_image(flo)
        nac_reg.process()
        tm_nac = nac_reg.get_transformation_matrix_forward()
        tm_nac.write(path_moco + 'tm_nac_' + str(i))
        i += 1


#%% data path and set filenames

py_path = os.getcwd()

# avoid cluttering of files, delete working-folder and start from scratch
working_folder = py_path + '/working'

# change the current working directory to the given path
os.chdir(working_folder)


#%% create folders for results

path_NAC = working_folder + '/recon/NAC/'
path_moco = working_folder + '/moco/'

if not os.path.exists(path_NAC):
    os.makedirs(path_NAC, mode=0o770)
    print('Create Folder: {}'.format(path_NAC))
if not os.path.exists(path_moco):
    os.makedirs(path_moco, mode=0o770)
    print('Create Folder: {}'.format(path_moco))


#%% Nifty registration NAC
# define reference image (first image) and float-path, NAC

# list of all NACs
list_NACs = [f for f in os.listdir(path_NAC) if f.endswith(".nii")]

tprint('Start Reg for NAC')

ref_file = path_NAC + 'NAC_0.nii'
ref = Eng_ref.ImageData(ref_file)

flo_path = path_NAC

# Niftyreg with NAC images
reg_nac(ref, flo_path, list_NACs)

tprint('Finish Reg for NAC')


translation_values = []
rotation_values = []

#%% get TM from Reg and calculate mean translation amplitude for every NAC volume
for file in sorted_alphanumeric(os.listdir(path_moco)):
    print(path_moco + file)
    tm_fwD_arr = numpy.loadtxt(path_moco + file)
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