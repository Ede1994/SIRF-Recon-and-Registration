#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Copyright University College London Hospital, 2020
Author: Eric Einspaenner, Institute of Nuclear Medicine
For internal research only.
"""
import os
import re
import numpy
from art import tprint
import sirf.STIR as Pet
import sirf.Reg as Reg
import sirf.Reg as Eng_ref
import sirf.Reg as Eng_flo


#%% Functions

# function for sorting lists, containing strings with numbers (https://stackoverflow.com/a/48413757)
def sorted_alphanumeric(data):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(data, key=alphanum_key, reverse=False)

# generate a list, contains the time intervalls (sum of frames)
def sum_list(lists):
    sum_liste = []
    add = int(0)
    sum_liste.append(add)
    for frame in lists[0:]:
        add += frame
        sum_liste.append(add)
    return sum_liste


#%% Definition of different path

py_path = os.getcwd()

# avoid cluttering of files, delete working-folder and start from scratch
working_folder = py_path + '/working2'

# change the current working directory to the given path
os.chdir(working_folder)

#%% Create folders for results

path_EPI = working_folder + '/EPI/'
path_mu = working_folder + '/mu/'
path_mu2 = working_folder + '/mu2/'

if not os.path.exists(path_mu):
    os.makedirs(path_mu, mode=0o770)
    print('Create Folder: {}'.format(path_mu))
if not os.path.exists(path_mu2):
    os.makedirs(path_mu2, mode=0o770)
    print('Create Folder: {}'.format(path_mu2))


#%% Define time frames (read frames.txt) and time intervals

# path to folder with frames.txt
frames_path = py_path + '/UKL_data/frames/'

# read frames.txt
with open(frames_path + 'frames.txt', 'r') as f:
    time_frames = [int(line.rstrip()) for line in f]
print('Frames: {}'.format(time_frames))

# number of motion steps
num_motion_steps = len(time_frames)

# define time intervals
time_intervals = sum_list(time_frames)
print('Time intervals: {}'.format(time_intervals))


# define the tm in the middle of the motion states
num_tm = []
for time, dur in zip(time_intervals, time_frames):
    t = int(time/2 + dur/4)
    num_tm.append(t)
print('Correct TM: {}'.format(num_tm))


#%% Files
    
attn_file = py_path + '/UKL_data/mu_Map/stir_mu_map.hv'  # .nii possible, requires ITK
print('mu-Map: {}'.format(attn_file))

# template for acq_data
template_acq_data = Pet.AcquisitionData('Siemens_mMR', span=11, max_ring_diff=16, view_mash_factor=1)
template_acq_data.write('template.hs')


#%% resample mu-Map into correct space and transform via invers tm
tprint('Start Resampling')

attn_image = Pet.ImageData(attn_file)

# template refernce
template_image = template_acq_data.create_uniform_image(1.0)

# EPI refernce file
epi_data_path = py_path + '/UKL_data/EPI/1/'
epi_file = epi_data_path + 'epi_0.nii'
epi = Eng_ref.ImageData(epi_file)

# define space matrices
tm_fwd = numpy.loadtxt(py_path + '/UKL_data/tm_epi/reg_NAC_EPI.txt')
tm_inv = numpy.loadtxt(py_path + '/UKL_data/tm_epi/reg_NAC_EPI_inv.txt')

tm_shift = numpy.loadtxt(py_path + '/UKL_data/tm_epi/tm_trans.txt')

# settings for attn resampler
resamplers_attn = Reg.NiftyResample()
resamplers_attn.set_reference_image(epi)
#resamplers_attn.set_floating_image(attn_image)
resamplers_attn.set_padding_value(0)
resamplers_attn.set_interpolation_type_to_linear()

i = 0
for num in num_tm:
    print('Begin resampling mu-Maps: {}'.format(path_EPI + 'tm_epi_' + str(num) + '.txt'))
    
    # read tm-matrix as numpy array
    matrix = numpy.loadtxt(path_EPI + 'tm_epi_inv_' + str(num) + '.txt')
    
    # create affine transformation from numpy array
    tm = Reg.AffineTransformation(matrix)
    tm2 = Reg.AffineTransformation(tm_shift)

    resamplers_attn.clear_transformations()
    resamplers_attn.set_floating_image(attn_image)
    resamplers_attn.add_transformation(tm2)
    new_attn = resamplers_attn.forward(attn_image)
    
    
    resamplers_attn.clear_transformations()
    resamplers_attn.set_floating_image(new_attn)
    resamplers_attn.add_transformation(tm)
    new_attn2 = resamplers_attn.forward(new_attn)
    
    new_attn2.write(path_mu + 'stir_mu_map_in_recon_space_' + str(i))
    Reg.ImageData(new_attn2).write(path_mu + 'mu_' + str(i) + '.nii')
    
    print('Finish resampling mu-Maps: {}'.format(i))
    
    i += 1

tprint('Finish Resampling')


#%% resampling back to STIR space

# list of all uMaps
list_uMap = [f for f in os.listdir(path_mu) if f.endswith(".hv")]

for i, mu in zip(range(len(list_uMap)), sorted_alphanumeric(list_uMap)):
    print('Begin resampling mu-Maps: {}'.format(mu))
    
    attn_image = Pet.ImageData(path_mu + mu)

    # settings for attn resampler
    resamplers_attn2 = Reg.NiftyResample()
    resamplers_attn2.set_reference_image(template_image)
    resamplers_attn2.set_floating_image(attn_image)
    resamplers_attn2.set_padding_value(0)
    resamplers_attn2.set_interpolation_type_to_linear()
    new_attn = resamplers_attn2.forward(attn_image)
    new_attn.write(path_mu2 + 'stir_mu_map_in_recon_space_' + str(i))
    Reg.ImageData(new_attn).write(path_mu2 + 'mu_' + str(i) + '.nii')
    
    print('Finish resampling mu-Maps: {}'.format(i))

    
