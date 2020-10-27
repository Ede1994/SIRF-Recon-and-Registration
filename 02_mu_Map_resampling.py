#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Copyright University College London Hospital, 2020
Author: Eric Einspaenner, Institute of Nuclear Medicine
For internal research only.
"""
import os
import numpy
from art import tprint
import sirf.STIR as Pet
import sirf.Reg as Reg


#%% Functions

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
working_folder = py_path + '/working'

# change the current working directory to the given path
os.chdir(working_folder)

#%% Create folders for results

path_EPI = working_folder + '/EPI/'
path_mu = working_folder + '/mu/'

if not os.path.exists(path_mu):
    os.makedirs(path_mu, mode=0o770)
    print('Create Folder: {}'.format(path_mu))


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
template_image = template_acq_data.create_uniform_image(1.0)

# define space matrices
tm_fwd = numpy.loadtxt(py_path + '/UKL_data/tm_epi/reg_NAC_EPI.txt')
tm_inv = numpy.loadtxt(py_path + '/UKL_data/tm_epi/reg_NAC_EPI_inv.txt')

# settings for attn resampler
resamplers_attn = Reg.NiftyResample()
resamplers_attn.set_reference_image(template_image)
resamplers_attn.set_floating_image(attn_image)
resamplers_attn.set_padding_value(0)
resamplers_attn.set_interpolation_type_to_linear()

i = 0
for num in num_tm:
    print('Begin resampling mu-Maps: {}'.format(path_EPI + 'tm_epi_inv_' + str(num) + '.txt'))
    
    # read matrix and calculate invers
    #matrix = numpy.loadtxt(path_EPI + 'tm_epi_' + str(num) + '.txt')
    matrix = numpy.loadtxt(path_EPI + 'tm_epi_inv_' + str(num) + '.txt')
    
    # tm space transformation: EPI -> NAC
    # transform tm into PET space: T_nac = R⁻1 * T_epi * R
    matrix1 = tm_fwd.dot(matrix)
    matrix2 = matrix1.dot(tm_inv)
    
    # create affine transformation from numpy array
    tm = Reg.AffineTransformation(matrix2)

    resamplers_attn.clear_transformations()
    resamplers_attn.add_transformation(tm)
    new_attn = resamplers_attn.forward(attn_image)
    new_attn.write(path_mu + 'stir_mu_map_in_recon_space_' + str(i))
    Reg.ImageData(new_attn).write(path_mu + 'mu_' + str(i) + '.nii')

    print('Finish resampling mu-Maps: {}'.format(i))
    
    i += 1

tprint('Finish Resampling')
