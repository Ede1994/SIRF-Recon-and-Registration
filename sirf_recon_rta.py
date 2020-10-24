# Copyright University College London Hospital, 2020
# Author: Eric Einspaenner, Institute of Nuclear Medicine
# For internal research only.

#%% import all necessary modules
import os
import re
import shutil
import numpy
from art import tprint
from sirf import SIRF
import sirf.STIR as Pet
import sirf.Reg as Reg
import sirf.Reg as Eng_ref
import sirf.Reg as Eng_flo


#%% functions

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


# registration of EPI images
def reg_epi(ref, flo_path):
    i = 0
    resampler2 = Reg.NiftyAladinSym()
    resampler2.set_reference_image(ref)
    for image in sorted_alphanumeric(os.listdir(flo_path)):
        flo_file = flo_path + image
        print(flo_file)
        flo = Eng_flo.ImageData(flo_file)
        resampler2.set_floating_image(flo)
        resampler2.process()
        tm_epi = resampler2.get_transformation_matrix_forward()
        tm_epi.write(path_EPI + 'tm_epi_' + str(i))
        i += 1


# Set the STIR verbosity
Pet.set_verbosity(1)


#%% data path and set filenames

py = 'sirf_recon_resample.py'
py_path = os.getcwd()
print(py_path)

# data path to list-mode and normalisation files
data_path_LM = py_path + '/UKL_data/LM/20190723/20190723/'
#data_path_MR = py_path + '/UKL_data/Processed/MedPhys_MoCo_Test_20190417_7/'

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
norm_file = data_path_LM + '1.3.12.2.1107.5.2.38.51008.30000019072305114947500000008.n.hdr'
print('Norm data: {}'.format(norm_file))
attn_file = py_path + '/UKL_data/mu_Map/stir_mu_map.hv'  # .nii possible, requires ITK
print('mu-Map: {}'.format(attn_file))

# output filename prefixes
sino_file = 'sino'


#%% settings for reconstruction
# set acq_model
acq_model = Pet.AcquisitionModelUsingRayTracingMatrix()
acq_model.set_num_tangential_LORs(5)

# set recon, OSEM
recon = Pet.OSMAPOSLReconstructor()
num_subsets = 7
num_subiterations = 4
recon.set_num_subsets(num_subsets)
recon.set_num_subiterations(num_subiterations)

# definitions for attenuation
attn_acq_model = Pet.AcquisitionModelUsingRayTracingMatrix()

# definitions for detector sensitivity modelling
asm_norm = Pet.AcquisitionSensitivityModel(norm_file)
acq_model.set_acquisition_sensitivity(asm_norm)


#%% redirect STIR messages to some files
# you can check these if things go wrong
msg_red = Pet.MessageRedirector('info.txt', 'warn.txt')


#%% create template and set lm2sino converter
# template for acq_data
template_acq_data = Pet.AcquisitionData('Siemens_mMR', span=11, max_ring_diff=16, view_mash_factor=2)
template_acq_data.write('template.hs')


#%% create listmode-to-sinograms converter object
lm2sino = Pet.ListmodeToSinograms()

# set input, output and template files
lm2sino.set_input(list_file)
lm2sino.set_output_prefix(sino_file)
lm2sino.set_template('template.hs')


#%% create folders for results
path_sino = working_folder + '/sino/'
path_rando = working_folder + '/rando/'
path_recon = working_folder + '/recon/'
path_NAC = working_folder + '/recon/NAC/'
path_EPI = working_folder + '/EPI/'
path_mu = working_folder + '/mu/'
path_tm = working_folder + '/tm/'
path_moco = working_folder + '/moco/'
if not os.path.exists(path_sino):
    os.makedirs(path_sino, mode=0o770)
    print('Create Folder: {}'.format(path_sino))
if not os.path.exists(path_rando):
    os.makedirs(path_rando, mode=0o770)
    print('Create Folder: {}'.format(path_rando))
if not os.path.exists(path_recon):
    os.makedirs(path_recon, mode=0o770)
    print('Create Folder: {}'.format(path_recon))
if not os.path.exists(path_NAC):
    os.makedirs(path_NAC, mode=0o770)
    print('Create Folder: {}'.format(path_NAC))
if not os.path.exists(path_EPI):
    os.makedirs(path_EPI, mode=0o770)
    print('Create Folder: {}'.format(path_EPI))
if not os.path.exists(path_mu):
    os.makedirs(path_mu, mode=0o770)
    print('Create Folder: {}'.format(path_mu))
if not os.path.exists(path_tm):
    os.makedirs(path_tm, mode=0o770)
    print('Create Folder: {}'.format(path_tm))
if not os.path.exists(path_moco):
    os.makedirs(path_moco, mode=0o770)
    print('Create Folder: {}'.format(path_moco))


#%% define time frames (read frames.txt) and time intervals
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


#%% Registration EPI, delivers transformation matrices 
# define reference image (first image) and float-path

tprint('Start Registration of EPIs')

# refernce file
epi_data_path = py_path + '/UKL_data/EPI/1/'
ref_file = epi_data_path + 'epi_0.nii'
ref = Eng_ref.ImageData(ref_file)

# float files
flo_path = epi_data_path

# Niftyreg with EPI images
reg_epi(ref, flo_path)

tprint('Finish Registration')


#%% space transformation NAC EPI
# define space matrices
tm_fwd = numpy.loadtxt(py_path + '/UKL_data/tm_epi/tm_trans.txt')
tm_inv = numpy.loadtxt(py_path + '/UKL_data/tm_epi/tm_trans_inv.txt')

tm_nacs = [0]*num_motion_steps

# transform new TM matrices into PET space and save as file and in list
n = 0
for item in num_tm:
    tm_epi = numpy.loadtxt(path_EPI + sorted_alphanumeric(os.listdir(path_EPI))[item])
    tm_nac = tm_epi #tm_inv * tm_epi * tm_fwd
    numpy.savetxt(path_tm + 'tm_' + str(item), tm_nac)
    tm_nacs[n] = tm_nac
    n += 1


#%% resample mu-Map into correct space and transform via invers tm
tprint('Start Resampling mu-Maps')

attn_image = Pet.ImageData(attn_file)
template_image = template_acq_data.create_uniform_image(1.0)

i = 0
for num in num_tm:
    print('Begin resampling mu-Maps: {}'.format(path_EPI + 'tm_epi_' + str(num) + '.txt'))
    
    # read matrix and calculate invers
    matrix = numpy.loadtxt(path_EPI + 'tm_epi_' + str(num) + '.txt')
    matrix2 = numpy.linalg.inv(matrix)
    
    # create affine transformation from numpy array
    tm = Reg.AffineTransformation(matrix2)
    
    resampler = Reg.NiftyResample()
    resampler.set_reference_image(template_image)
    resampler.set_floating_image(attn_image)
    resampler.add_transformation(tm)
    resampler.set_padding_value(0)
    resampler.set_interpolation_type_to_linear()
    attn_image = resampler.forward(attn_image)
    attn_image.write(path_mu + 'stir_mu_map_in_recon_space_' + str(i))

    print('Finish resampling mu-Maps: {}'.format(i))
    
    i += 1

tprint('Finish Resampling mu-Maps')


#%% for-loop, reconstruction of time intervals
path_mu_2 = [f for f in os.listdir(path_mu) if f.endswith(".hv")]

tprint('Start Recon')

for mu, i in zip(sorted_alphanumeric(path_mu_2), range(len(time_intervals)-1)):
    print('Begin reconstruction: Frame {}'.format(i))
    print('Time interval: {} - {}'.format(time_intervals[i], time_intervals[i+1]))

    # listmode-to-sinogram
    lm2sino.set_time_interval(time_intervals[i], time_intervals[i+1])
    lm2sino.set_up()
    lm2sino.process()
    acq_data = lm2sino.get_output()
    acq_data.write(working_folder + '/sino/sino'+str(i))

    # randoms estimate
    randoms = lm2sino.estimate_randoms()
    randoms.write(working_folder + '/rando/rando'+str(i))
    
    # choose correct uMap
    attn_image = Pet.ImageData(path_mu + mu)
    print(mu)
    asm_attn = Pet.AcquisitionSensitivityModel(attn_image, attn_acq_model)

    # reconstruct the data (includes all)
    obj_fun = Pet.make_Poisson_loglikelihood(acq_data)
    asm_attn.set_up(acq_data)
    attn_factors = Pet.AcquisitionData(acq_data)
    attn_factors.fill(1.0)
    asm_attn.unnormalise(attn_factors)
    asm_attn = Pet.AcquisitionSensitivityModel(attn_factors)
    asm = Pet.AcquisitionSensitivityModel(asm_norm, asm_attn)
    acq_model.set_acquisition_sensitivity(asm)
    acq_model.set_background_term(randoms)
    obj_fun.set_acquisition_model(acq_model)
    recon.set_objective_function(obj_fun)
    initial_image = acq_data.create_uniform_image(1.0)
    image = initial_image
    recon.set_up(image)

    recon.set_current_estimate(image)
    recon.process()

    # save recon images
    recon_image = recon.get_output()
    recon_image.write(working_folder + '/recon/recon' + str(i))

    # save Image as .nii
    recon_image = Reg.NiftiImageData(recon.get_output())
    recon_image.write(working_folder + '/floates/recon'+str(i))

    print('Reconstruction successful: Frame {}'.format(i))

tprint('Finish Recon')


#%% convert a array to a SIRF transformation matrix and then resample the float image
# define reference image and float-path
ref_file = working_folder + '/floates/recon0.nii'
ref = Eng_ref.ImageData(ref_file)
flo_path = working_folder + '/floates/'

tprint('Start Resampling')
i = 0
#sorted(os.listdir(path_EPI)),
# for loop, simultaneous matrices and images
for num, image in zip(num_tm, sorted_alphanumeric(os.listdir(flo_path))):
    print('TM: {}, Float-Image: {}'.format('tm_epi_' + str(num) + '.txt', image))

    # read tm-matrix as numpy array and read load float image
    matrix = numpy.loadtxt(path_EPI + 'tm_epi_' + str(num) + '.txt')
    flo = Eng_flo.ImageData(flo_path + image)

    # create affine transformation from numpy array
    tm = Reg.AffineTransformation(matrix)

    # motion information from another source and resample an image with this informations
    print('Begin resampling: {}'.format(image))
    resampler = Reg.NiftyResample()
    resampler.set_reference_image(ref)
    resampler.set_floating_image(flo)
    resampler.add_transformation(tm)
    resampler.set_padding_value(0)
    resampler.set_interpolation_type_to_linear()
    output = resampler.forward(flo)
    output.write(working_folder + '/moco/moco_'+str(i))

    print('Resampling successful: {}'.format(image))

    i += 1

tprint('Finish Resampling')


#%% define RTA method
# define initial image (first image, first frame)
initial = Reg.NiftiImageData(working_folder + '/moco/moco_0.nii')
initial_array = initial.as_array()

# sum over all images (as array)
for image in sorted_alphanumeric(os.listdir(path_moco))[1:]:
    print(image)
    array = Reg.NiftiImageData(path_moco + image).as_array()
    initial_array += array

print('Final array (RTA): {}'.format(initial_array))

# create image
final_image = Reg.NiftiImageData(working_folder + '/moco/moco_0.nii')
print(type(final_image))
final_image.fill(initial_array)
print(type(final_image))
final_image.write('final_image_RTA.nii')

tprint('DONE!')
