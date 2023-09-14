from psychopy import visual, core, event
import numpy as np
import psychopy.iohub as io
from psychopy.hardware import keyboard
import matplotlib.pyplot as plt
from PIL import Image
import os
from dva_to_pix import arcmin_to_px
from numpy.random import choice as randchoice
from numpy.random import random, randint, normal, shuffle, choice as randchoice
from psychopy import sound, gui, visual, core, data, event, logging, clock, colors, layout
from psychopy.constants import (NOT_STARTED, STARTED, PLAYING, PAUSED,
                                STOPPED, FINISHED, PRESSED, RELEASED, FOREVER)
from psychopy.tools import monitorunittools
import sys  # to get file system encoding

from create_conditions import condition_creater
from psychopy import prefs

#prefs.hardware['audioLib'] = ['PTB']
conditions= condition_creater()
# import conditions.npy
#conditions=np.load('conditions.npy')
# Create PsychoPy window covering the whole screen
win = visual.Window(size=(512, 512), fullscr=False, monitor='testMonitor', units='pix', color=[0, 0, 0], useFBO=True)
field_size=(512,512)

#setup screen properties
screen_width=28
screen_height=17
screen_distance=50
### Set the monitor to the correct distance and size
#win.monitor.setSizePix(field_size)
win.monitor.setWidth(screen_width)
win.monitor.setDistance(screen_distance)
mouse = event.Mouse(win=win,visible=False)
# TODO: arcmin to px conversion but open it later
# def arcmin_to_px(arcmin=1,h=19,d=57,r=1080):
#     dva= arcmin/60
#     return(monitorunittools.deg2pix(degrees=dva,monitor=win.monitor))

# Store info about the experiment session
psychopyVersion = '2022.2.4'
expName = 'continous_psych'  # from the Builder filename that created this script
expInfo = {
    'participant': f"{randint(0, 999999):06.0f}",
    'session': '001',
}
expInfo['date'] = data.getDateStr()  # add a simple timestamp
expInfo['expName'] = expName
expInfo['psychopyVersion'] = psychopyVersion

frameRate=win.getActualFrameRate()
print(frameRate)
expInfo['frameRate']=frameRate
if expInfo['frameRate'] != None:
    frameDur = 1.0 / round(expInfo['frameRate'])
else:
    frameDur = 1.0 / 60.0  # could not measure, so guess
defaultKeyboard = keyboard.Keyboard(backend='iohub')
endExpNow = False  # flag for 'escape' or other condition => quit the exp

# Ensure that relative paths start from the same directory as this script
_thisDir = os.path.dirname(os.path.abspath(__file__))
os.chdir(_thisDir)
# Data file name stem = absolute path + name; later add .psyexp, .csv, .log, etc
filename = _thisDir + os.sep + u'data/%s_%s_%s' % (expInfo['participant'], expName, expInfo['date'])


endExpNow = False  # flag for 'escape' or other condition => quit the exp
frameTolerance = 0.001  # how close to onset before 'same' frame
# Noise properties
noise_size = field_size  # Use the window size for noise texture
noise_arcmin = 11  # Standard deviation for pixel noise || noise intensity Adjust to control noise intensity
noise_std =arcmin_to_px(noise_arcmin,h=screen_height,d=screen_distance,r=field_size[0])# convert arcmin to std for Gaussian
# Brownian motion properties
velocity_std = 1.0  # Standard deviation of Gaussian white noise velocities
# Create some handy timers
globalClock = core.Clock()  # to track the time since experiment started
routineTimer = core.Clock()  # to track time remaining of each (possibly non-slip) routine 


########################################################################################################################
# Background noise
def gaussNoise(noise_intensity=1):
    noise = np.random.normal(0, noise_intensity, size=noise_size)
    return noise

# Blob properties
initial_blob_std= arcmin_to_px(arcmin=11,h=screen_height,d=screen_distance,r=field_size[0]) 
total_lum=1*(2*np.pi*initial_blob_std ** 2)*1 # Total luminance of blobm
def generateBlob(space_constant):
    blob_std=arcmin_to_px(arcmin=blob_width,h=screen_height,d=screen_distance,r=field_size[0])
    blob_amplitude =total_lum / (2*np.pi*blob_std ** 2) # Adjust blob amplitude to keep total blob energy constant
    y, x = np.meshgrid(np.arange(noise_size[1]), np.arange(noise_size[0]))
    blob = blob_amplitude * np.exp(-((x - noise_size[0] / 2)**2 + (y - noise_size[1] / 2)**2) / (2 * blob_std**2))
    #blob = (blob - blob.min()) / (blob.max() - blob.min())*2-1 # normalize the blob
    return blob

## create blobs for each of the conditions
blob_widths=[11,13,17,21,25,29]
blobs=[]
for blob_width in blob_widths:
    blobs.append(np.array(generateBlob(arcmin_to_px(arcmin=blob_width,h=screen_height,d=screen_distance,r=field_size[0]))))
blobs=np.array(blobs)
blobs = (blobs - blobs.min()) / (blobs.max() - blobs.min())*2-1
# create dictionary for each blob_width
blob_dict={}
for i in range(len(blob_widths)):
    blob_dict[blob_widths[i]]=blobs[i]


magicNumber=600
expectedFrameRate=60
expectedDuration=20 # in seconds
expectedFrames=expectedFrameRate*expectedDuration

# Blob motion by changing velocity on x and y axis according to the brownian motion
def generateBrownianMotion(field_size, velocity_std, duration):
    # Generate random velocity for each frame
    num_frames = int(duration * expectedFrameRate)+5
    velocies_x=np.random.normal(0, velocity_std, num_frames)
    velocies_y=np.random.normal(0, velocity_std, num_frames)
    # new positions = old position + velocity
    pos_x = np.cumsum(velocies_x)
    pos_y = np.cumsum(velocies_y)
    # clip positions to stay within the field
    pos_x = np.clip(pos_x, -field_size/2, field_size/2 )
    pos_y = np.clip(pos_y, -field_size/2, field_size/2 )
    return (pos_x, pos_y)

# pre-generate noises
noiseQuantity=120
noise_instances=np.empty((noiseQuantity,noise_size[0],noise_size[1]))
for i in range(noiseQuantity):
    noise_instances[i]=gaussNoise()
noise_instances = np.clip(noise_instances, -3, 3)    
noise_instances = (noise_instances - noise_instances.min()) / (noise_instances.max() - noise_instances.min())*2-1 # normalize the noise

    
# Observation pointer for participant tracking of blob
#circle showing mouse position
obs_pointer = visual.Circle(win, radius=10, fillColor=[1, 0, 0], units='pix',size=0.35)
#Mouse properties
mouse.mouseClock = core.Clock()

def noise_visual_func(noiseQuantity):
    noise_obj=[]
    for i in range(noiseQuantity):
        noise_obj.append( visual.GratingStim(win, tex=noise_instances[i], size=noise_size, interpolate=False,units='pix'))
        noise_obj[-1].draw()    # draw noise
    return noise_obj

noise_obj=noise_visual_func(noiseQuantity)

##################### Trials Start Here #####################
# data to save
## blob data
all_blob_x = []
all_blob_y = []
all_blob_v=[]
### mouse data
all_mouse_x = []
all_mouse_y = []
all_mouse_v=[]

def black_to_transparent(img):
    num_channels = img.shape[2]
    dimensions = img.shape[0:2]
    alpha_channel = (np.ones(dimensions) * 255).astype(np.uint8)
    alpha_img = np.dstack((img, alpha_channel ))
    mask = np.sum(alpha_img[:, :, 0:img.shape[2]], axis=2).clip(0,255)
    alpha_img[:,:,num_channels] = mask
    return alpha_img

win.setMouseVisible(False)        
for blob_width in (blob_widths):#conditions:
    # set a timer for the next line record the time spent
    tStart = globalClock.getTime()
    # clear all drawings
    print(blob_width)
    blob_std=arcmin_to_px(arcmin=blob_width,h=screen_height,d=screen_distance,r=field_size[0])
    # create blob
    #blob=generateBlob(blob_std)
    #blob=blobs[int(np.where(np.array(blob_widths)==blob_width)[0])]
    blob=blob_dict[blob_width]

    blob_obj = visual.GratingStim(win, tex=None, mask=blob,  units='pix', interpolate=False)


    # create observation pointer
    obs_pointer.setAutoDraw(False)
    mouse.setPos(newPos=(0,0))
    win.flip(clearBuffer=True)

    mouse.y=[]
    mouse.x=[]
    # create magicNumber frames of stimuli
    # create blob motion positions
    pos_x_obj, pos_y_obj = generateBrownianMotion(field_size=noise_size[0], velocity_std=1, duration=expectedDuration)# pre-generate brownian motion
    pos_x_obj=np.insert(pos_x_obj,0,0)
    pos_y_obj=np.insert(pos_y_obj,0,0)
    # calcuate velocity of blob
    blob_v = np.sqrt(np.diff(pos_x_obj)**2 + np.diff(pos_y_obj)**2)
    blob_v = np.insert(blob_v, 0, 0)
    # save blob velocity
    all_blob_x.append(pos_x_obj)
    all_blob_y.append(pos_y_obj)
    all_blob_v.append(blob_v)

    tEnd = globalClock.getTime()
    win.clearBuffer()
    intensity_profiles = []  # List to store intensity profiles
    continueRoutine = True
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    frameN = -1
    routineTimer.reset()
    t = 0

    # # wait for 1 second before starting the trial
    waiterTime=0
    while  waiterTime <0.5:
        waiterTime = globalClock.getTime() - tStart
        win.flip()

    t = 0
    # draw the stimulus
    obs_pointer.setAutoDraw(True)
    blob_obj.setAutoDraw(True)
    mouse.setVisible(False)  
    ##################### Trial Start #####################
    while continueRoutine:
        t = routineTimer.getTime()
        #tThisFlip = win.getFutureFlipTime(clock=routineTimer)
        #tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN+= 1  # number of completed frames (so 0 is the first frame)

        # draw the stimulus
        obs_pointer.setPos(mouse.getPos())
        noise_obj[int(randint(0,noiseQuantity))].draw()    # draw noise
        blob_obj.setPos((pos_x_obj[frameN], pos_y_obj[frameN]))

        #save mouse position
        mouse.x.append(obs_pointer.pos[0])
        mouse.y.append(obs_pointer.pos[1])
        # flip the window
        win.flip()
        # end the loop after given seconds
        if frameN > expectedFrames-1:
            continueRoutine = False
        # check for quit (typically the Esc key)
        if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
            core.quit()
 
    blob_obj.setAutoDraw(False)

    t = routineTimer.getTime()
    print(t)
    print(tEnd-tStart)
    print(frameN)

    # save positions of target, mouse, and but first calculate velocities
    mouse.x = np.array(mouse.x)
    mouse.y = np.array(mouse.y)
    mouse.v = np.sqrt(np.diff(mouse.x)**2 + np.diff(mouse.y)**2)
    mouse.v = np.insert(mouse.v, 0, 0)
    all_mouse_x.append(mouse.x)
    all_mouse_y.append(mouse.y)
    all_mouse_v.append(mouse.v)
    ###
    # save target velocity


    # save mouse positions

# Close the window
event.waitKeys()
win.close()

# write mouse positions, velocities, and blob velocities and positions and conditions as a matlab data file
# save mouse positions

import scipy.io as sio
conditions=conditions.T
sio.savemat('recorded/all.mat', {'mouse_x': all_mouse_x, 'mouse_y': all_mouse_y, 'mouse_v': all_mouse_v, 'blob_x': all_blob_x, 'blob_y': all_blob_y, 'blob_v': all_blob_v, 'blob_width': conditions})

# transpoze conditions






# ################################# Plotting #############################################################################

# # save intensity profiles based on final_stims array
# intensity_profiles = []  # List to store intensity profiles
# for frameN,frame in enumerate(final_stims):
#     #stim_array = np.array(frame.tex)
#     intensity_profile = frame.tex[noise_size[1]//2,:]
#     normalized_profile = (intensity_profile - intensity_profile.min()) / (intensity_profile.max() - intensity_profile.min())
#     intensity_profiles.append(normalized_profile)

# # function to plot the intensity profile for each frame
# def plot_intensity_profile(intensity_profiles, frameN):
#     plt.figure()
#     plt.plot(intensity_profiles)
#     plt.xlabel("Horizontal Position (pixels)")
#     plt.ylabel("Intensity")
#     plt.title("Cross-Sections of Intensity"+ "  Frame: " + str(frameN))
#     plt.show()

# # plot intensity profile for each frame
# def each_frame_intensities(blob_width, intensity_profiles, plot_intensity_profile):
#     for frameN in range(len(intensity_profiles)):
#         plot_intensity_profile(intensity_profiles[frameN], frameN)
#         #plt.pause(0.0001)
#         if frameN == 0:
#             plt.savefig('recorded/intensity_profile_'+str(blob_width)+'.png')

# #each_frame_intensities(blob_width, intensity_profiles[0:10], plot_intensity_profile)



# # Plot intensity profiles
# def intensities_normalized(noise_size, blob_width, intensity_profiles):
#     plt.figure(figsize=(8, 6))

#     for profile in intensity_profiles:
#         normalized_profile = (profile - profile.min()) / (profile.max() - profile.min())
#         plt.plot(normalized_profile)
# # x label in degrees of visual angle
#     plt.xticks(np.arange(0, noise_size[0], 100), np.arange(-5, 6, 1))
#     plt.xlabel("Horizontal Position (degrees)")
#     plt.ylabel("Normalized Intensity")
#     plt.title("Cross-Sections of Normalized Intensity")
#     plt.ylim(0, 1)
# # legend for each frame
#     plt.legend(np.arange(len(intensity_profiles)))
#     plt.savefig('recorded/intensity_profiles_'+str(blob_width)+'.png')
#     plt.show()


# # intensities_normalized(noise_size, blob_width, intensity_profiles)

# # plot intensity profiles mean
# def mean_profile(intensity_profiles):
#     mean_profile = np.mean(intensity_profiles, axis=0)
#     normalized_mean_profile = (mean_profile - mean_profile.min()) / (mean_profile.max() - mean_profile.min())
#     return normalized_mean_profile

# # function to plot intensity profiles
# def plot_intensity_profiles(intensity_profiles):
#     normalized_mean_profile = mean_profile(intensity_profiles)
#     plt.figure(figsize=(8, 6))
#     # x label in degrees of visual angle
#     plt.xticks(np.arange(0, noise_size[0], 100), np.arange(-5, 6, 1))
#     plt.xlabel("Horizontal Position (degrees)")
#     plt.ylabel("Normalized Intensity")
#     plt.title("Cross-Sections of Normalized Intensity")
#     plt.ylim(0, 1)
#     # legend for each frame
#     plt.legend(np.arange(len(intensity_profiles)))
#     plt.savefig('recorded/intensity_profiles_'+str(blob_width)+'.png')
#     plt.show()

# # save intensity profiles as csv
# np.savetxt("recorded/intensity_profiles.csv", intensity_profiles, delimiter=",")
# # save intensity profiles as numpy array
# np.save("recorded/intensity_profiles.npy", intensity_profiles)
# # save mean intensity profile as numpy array
# np.save("recorded/mean_intensity_profile.npy", mean_profile)