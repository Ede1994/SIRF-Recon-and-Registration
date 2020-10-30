#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Copyright University College London Hospital, 2020
Author: Eric Einspaenner, Institute of Nuclear Medicine
For internal research only.
"""
import os
import shutil
from art import tprint
import sirf.STIR as Pet


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
print(py_path)

# data path to list-mode and normalisation files
data_path_LM = '/home/eric/Dokumente/PythonProjects/SIRF/UKL_data/LM/20190723/20190723/'

# avoid cluttering of files, delete working-folder and start from scratch
working_folder = py_path + '/working'
if os.path.exists(working_folder):
    shutil.rmtree(working_folder)
if not os.path.exists(working_folder):
    os.makedirs(working_folder, mode=0o770)

# change the current working directory to the given path
os.chdir(working_folder)

# input files
list_file = data_path_LM + '1.3.12.2.1107.5.2.38.51008.30000019072305114947500000006.l.hdr'
print('LM data: {}'.format(list_file))

# output filename prefixes
sino_file = 'sino'


#%% Create folders for results

path_sino = working_folder + '/sino/'
path_rando = working_folder + '/rando/'

if not os.path.exists(path_sino):
    os.makedirs(path_sino, mode=0o770)
    print('Create Folder: {}'.format(path_sino))
if not os.path.exists(path_rando):
    os.makedirs(path_rando, mode=0o770)
    print('Create Folder: {}'.format(path_rando))


#%% create template and set lm2sino converter
# template for acq_data
template_acq_data = Pet.AcquisitionData('Siemens_mMR', span=11, max_ring_diff=16, view_mash_factor=1)
template_acq_data.write('template.hs')


#%% Create listmode-to-sinograms converter object

lm2sino = Pet.ListmodeToSinograms()

# set input, output and template files
lm2sino.set_input(list_file)
lm2sino.set_output_prefix(sino_file)
lm2sino.set_template('template.hs')


#%% Define time frames (read frames.txt) and time intervals
    
# path to folder with frames.txt
frames_path = '/home/eric/Dokumente/PythonProjects/SIRF/UKL_data/frames/'

# read frames.txt
with open(frames_path + 'frames.txt', 'r') as f:
    time_frames = [int(line.rstrip()) for line in f]
print('Frames: {}'.format(time_frames))

# number of motion steps
num_motion_steps = len(time_frames)

# define time intervals
time_intervals = sum_list(time_frames)
print('Time intervals: {}'.format(time_intervals))


#%% Create sinograms and randoms

tprint('Start Sinograms')

for i in range(len(time_intervals)-1):
    print('Begin: Frame {}'.format(i))
    print('Time interval: {} - {}'.format(time_intervals[i], time_intervals[i+1]))

    # listmode-to-sinogram
    lm2sino.set_time_interval(time_intervals[i], time_intervals[i+1])
    lm2sino.set_up()
    lm2sino.process()
    acq_data = lm2sino.get_output()
    acq_data.write(path_sino + 'sino'+str(i))

    # randoms estimate
    randoms = lm2sino.estimate_randoms()
    randoms.write(path_rando + 'rando' + str(i))

    print('Finish successful: Frame {}'.format(i))

tprint('Finish Sinograms')