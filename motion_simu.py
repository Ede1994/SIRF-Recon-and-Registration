# Copyright University College London Hospital, 2020
# Author: Eric Einspaenner, Institute of Nuclear Medicine
# For internal research only.

#%% import all necessary modules
import os
import shutil
import sys
import glob
import re
import numpy
import matplotlib.pyplot as plt
from math import cos, sin, radians, degrees
from art import *
from sirf import SIRF
import sirf.STIR as Pet
import sirf.Reg as Reg
import sirf.Reg as Eng_ref
import sirf.Reg as Eng_flo
from nipype.interfaces import fsl


#%% functions
# convert degree in radians and calculate cos and sin
def trig(angle):
    r = radians(angle)
    return cos(r), sin(r)


#%% define data path and working dir
# data path to recon file
recon_file = '/home/rich/Documents/ReCo/working2/floates/test.nii'

# avoid cluttering of files, delete working-folder and start from scratch
working_folder = '/home/rich/Documents/ReCo/working2'
#if os.path.exists(working_folder):
#    shutil.rmtree(working_folder)
if not os.path.exists(working_folder):
    os.makedirs(working_folder, mode=0o770)

# change the current working directory to the given path
os.chdir(working_folder)


#%% define ref and float images
# ref = float, for motion simulation
ref_file = recon_file
ref = Eng_ref.ImageData(ref_file)

flo_file = recon_file
flo = Eng_flo.ImageData(flo_file)


#%% define motion parameters
# rotation
xC, xS = trig(0)
yC, yS = trig(0)
zC, zS = trig(30)

# translation
dX = 0
dY = 0
dZ = 0


#%% define TM
Translate = numpy.array([[1, 0, 0, dX],
                        [0, 1, 0, dY],
                        [0, 0, 1, dZ],
                        [0, 0, 0, 1]])

Rotate_X = numpy.array([[1, 0, 0, 0],
                        [0, xC, -xS, 0],
                        [0, xS, xC, 0],
                        [0, 0, 0, 1]])

Rotate_Y = numpy.array([[yC, 0, -yS, 0],
                        [0, 1, 0, 0],
                        [yS, 0, yC, 0],
                        [0, 0, 0, 1]])

Rotate_Z = numpy.array([[zC, -zS, 0, 0],
                        [zS, zC, 0, 0],
                        [0, 0, 1, 0],
                        [0, 0, 0, 1]])

TM = numpy.dot(Rotate_Z, numpy.dot(Rotate_Y, numpy.dot(Rotate_X, Translate)))


#%% start simulation motion
tprint('Start TM')

# define tm matrix
tm = Reg.AffineTransformation(TM)
tm_inv = tm.get_inverse()

# apply TM to ref image and resample an image with this informations
print('Begin resampling')
resampler = Reg.NiftyResample()
resampler.set_reference_image(ref)
resampler.set_floating_image(flo)
resampler.add_transformation(tm)
resampler.set_padding_value(0)
resampler.set_interpolation_type_to_linear()
output = resampler.forward(flo)
output.write(working_folder + '/floates/test_motion.nii')

print('Resampling successful')