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
# Create PsychoPy window covering the whole screen
win = visual.Window(size=(512, 512), fullscr=False, monitor='testMonitor', units='pix', color=[0, 0, 0], useFBO=True)
field_size=(256,256)

#setup screen properties
screen_width=34
screen_height=21
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
    noise = np.clip(noise, -3*noise_intensity, 3*noise_intensity)    
    noise = (noise - noise.min()) / (noise.max() - noise.min())*2-1 # normalize the noise
    return noise

# Blob properties
initial_blob_std= arcmin_to_px(arcmin=11,h=screen_height,d=screen_distance,r=field_size[0]) 
total_lum=1*(2*np.pi*initial_blob_std ** 2)*1 # Total luminance of blobm
def generateBlob(space_constant):
    #blob_std=arcmin_to_px(arcmin=blob_width,h=screen_height,d=screen_distance,r=field_size[0])
    blob_amplitude =total_lum / (2*np.pi*blob_std ** 2) # Adjust blob amplitude to keep total blob energy constant
    y, x = np.meshgrid(np.arange(noise_size[1]), np.arange(noise_size[0]))
    blob = blob_amplitude * np.exp(-((x - noise_size[0] / 2)**2 + (y - noise_size[1] / 2)**2) / (2 * blob_std**2))
    return blob


###################################################################################################################
"""                         Pregenerate noise and blob instances for each frame """
###################################################################################################################
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
noise_instances=np.empty((int(expectedFrameRate*expectedDuration)+5,noise_size[0],noise_size[1]))
for i in range((expectedFrameRate*expectedDuration)+5):
    noise_instances[i]=gaussNoise()
# ## create visual patches for each frame
def full_stimuli(blob_width,expectedDuration,frameRate):
    """ Pregenerate noise and blob instances for each frame """
    noise= noise_instances[-1]# first noise
    blob=generateBlob(blob_width)    # Create Gaussian blob
    pos_x, pos_y = generateBrownianMotion(field_size=noise_size[0], velocity_std=1, duration=expectedDuration)# pre-generate brownian motion
    stim_arrays=np.empty((magicNumber, noise_size[0], noise_size[0]))
    stim_arrays[0]=blob+noise
    rolled_blob_frames = np.empty((magicNumber, noise_size[0], noise_size[0]))
    final_stims = [None] * ((expectedDuration*frameRate)+1)
    final_stims[0] = visual.GratingStim(win, tex=stim_arrays[0], size=noise_size, interpolate=False,units='pix')
    final_stims[0].draw()  # first draw is slower. So do it now.
    for i in range(1,magicNumber):#int(expectedDuration*frameRate)+5): #assuming program can achieve the actual frame rate of screen
        #noise_instances.append(gaussNoise())# pre-generated noise
        rolled_blob_frames[i]=np.roll(blob, (int(pos_x[(i)]), int(pos_y[(i)])), axis=(1, 0)) # Shift blob according to the brownian motion
        stim_arrays[i]=rolled_blob_frames[i]+noise_instances[i] # Combine noise with the blob
    #stim_arrays = np.clip(stim_arrays, -3*noise_std, 3*noise_std)        #clip the stim array to stay within the range of 3*noise_intensity
    #stim_arrays = (stim_arrays - stim_arrays.min()) / (stim_arrays.max() - stim_arrays.min())*2-1 # normalize the stim array
    for i in range(1,(magicNumber)):
        final_stims[i] = visual.GratingStim(win, tex=stim_arrays[i], size=[512,512], interpolate=False,units='pix',autoLog=False)
        # another way to optimize the creation of visual stimuli
        final_stims[i].draw()  # first draw is slower. So do it now.
    return final_stims, pos_x, pos_y,blob,rolled_blob_frames

def full_stim_for_single_frame(frameOrder,pos_x,pos_y,blob):
    blob_moved=np.roll(blob, (int(pos_x[(frameOrder)]), int(pos_y[(frameOrder)])), axis=(1, 0)) # Shift blob according to the brownian motion
    stim_array=blob_moved+noise_instances[frameOrder] # Combine noise with the blob
    #stim_array = np.clip(stim_array, -3*noise_std, 3*noise_std)        #clip the stim array to stay within the range of 3*noise_intensity
    #stim_array = (stim_array - stim_array.min()) / (stim_array.max() - stim_array.min())*2-1 # normalize the stim array
    final_stim = visual.GratingStim(win, tex=stim_array, size=[512,512], interpolate=False,units='pix')
    #final_stim.draw()  # first draw is slower. So do it now.
    return final_stim


    
# Observation pointer for participant tracking of blob
#circle showing mouse position
obs_pointer = visual.Circle(win, radius=10, fillColor=[1, 0, 0], units='pix',size=0.35)
#Mouse properties
mouse.mouseClock = core.Clock()


##################### Trials Start Here #####################
all_blob_x = []
all_blob_y = []
all_blob_v=[]


all_mouse_x = []
all_mouse_y = []
all_mouse_v=[]

win.setMouseVisible(False)        
for blob_width in conditions:
    # clear all drawings
    print(blob_width)
    blob_std=arcmin_to_px(arcmin=blob_width,h=screen_height,d=screen_distance,r=field_size[0])
    obs_pointer.setAutoDraw(False)
    mouse.setPos(newPos=(0,0))
    win.flip(clearBuffer=True)
    # set a timer for the next line record the time spent
    tStart = globalClock.getTime()
    #mouse.x = np.zeros(expectedDuration*expectedFrameRate)
    #mouse.y = np.zeros(expectedDuration*expectedFrameRate)
    mouse.y=[]
    mouse.x=[]
    # create magicNumber frames of stimuli
    final_stims_obj, pos_x_obj, pos_y_obj, blob_obj, rolledBlobPositions=full_stimuli(blob_std,expectedDuration,expectedFrameRate)
    pos_x_obj=np.insert(pos_x_obj,0,0)
    pos_y_obj=np.insert(pos_y_obj,0,0)
    # calcuate velocity of blob
    blob_vx = np.diff(pos_x_obj)
    blob_vy = np.diff(pos_y_obj)
    blob_v = np.sqrt(blob_vx**2 + blob_vy**2)
    blob_v = np.insert(blob_v, 0, 0)
    # save blob velocity
    all_blob_v.append(blob_v)

    tEnd = globalClock.getTime()
    print(tEnd-tStart)

    win.clearBuffer()
    # update component parameters for each repeat
    # Initialize components for Routine "Trial"
    intensity_profiles = []  # List to store intensity profiles
    continueRoutine = True
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    frameN = -1
    routineTimer.reset()
    t = 0
    # create a blank screen for 0.3 seconds
    obs_pointer.setAutoDraw(True)
    ##################### Trial Start #####################
    while continueRoutine:
        mouse.setVisible(False)  
        t = routineTimer.getTime()
        tThisFlip = win.getFutureFlipTime(clock=routineTimer)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN+= 1  # number of completed frames (so 0 is the first frame)
        # draw the stimulus
        obs_pointer.setPos(mouse.getPos())
        final_stims_obj[frameN].draw()
        #print(pos_x_ob[frameN], pos_y_obj[frameN])
        #save mouse position
        mouse.x.append(obs_pointer.pos[0])
        mouse.y.append(obs_pointer.pos[1])
        #mouse.x[frameN]=obs_pointer.pos[0]
        #mouse.y[frameN]=obs_pointer.pos[1]
        #print(frameN)
        #print("time="+str(t))
        # flip the window
        win.flip()
        # end the loop after given seconds
        if frameN > expectedFrames-1:
            continueRoutine = False
        # check for quit (typically the Esc key)
        if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
            core.quit()
        if frameN < expectedFrames-magicNumber+1:
            blob_moved2nd=np.roll(blob_obj, (int(pos_x_obj[(frameN+magicNumber)]), int(pos_y_obj[(frameN+magicNumber)])), axis=(1, 0)) # Shift blob according to the brownian motion
            stim_array2nd=blob_moved2nd+noise_instances[frameN+magicNumber] # Combine noise with the blob
            final_stims_obj[frameN+magicNumber]=visual.GratingStim(win, tex=stim_array2nd, size=[512,512], interpolate=False,units='pix')
            #win.clearBuffer()
    print(tEnd-tStart)
    print(frameN)
    t = routineTimer.getTime()
    print(t)
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








################################# Plotting #############################################################################

# save intensity profiles based on final_stims array
intensity_profiles = []  # List to store intensity profiles
for frameN,frame in enumerate(final_stims):
    #stim_array = np.array(frame.tex)
    intensity_profile = frame.tex[noise_size[1]//2,:]
    normalized_profile = (intensity_profile - intensity_profile.min()) / (intensity_profile.max() - intensity_profile.min())
    intensity_profiles.append(normalized_profile)

# function to plot the intensity profile for each frame
def plot_intensity_profile(intensity_profiles, frameN):
    plt.figure()
    plt.plot(intensity_profiles)
    plt.xlabel("Horizontal Position (pixels)")
    plt.ylabel("Intensity")
    plt.title("Cross-Sections of Intensity"+ "  Frame: " + str(frameN))
    plt.show()

# plot intensity profile for each frame
def each_frame_intensities(blob_width, intensity_profiles, plot_intensity_profile):
    for frameN in range(len(intensity_profiles)):
        plot_intensity_profile(intensity_profiles[frameN], frameN)
        #plt.pause(0.0001)
        if frameN == 0:
            plt.savefig('recorded/intensity_profile_'+str(blob_width)+'.png')

#each_frame_intensities(blob_width, intensity_profiles[0:10], plot_intensity_profile)



# Plot intensity profiles
def intensities_normalized(noise_size, blob_width, intensity_profiles):
    plt.figure(figsize=(8, 6))

    for profile in intensity_profiles:
        normalized_profile = (profile - profile.min()) / (profile.max() - profile.min())
        plt.plot(normalized_profile)
# x label in degrees of visual angle
    plt.xticks(np.arange(0, noise_size[0], 100), np.arange(-5, 6, 1))
    plt.xlabel("Horizontal Position (degrees)")
    plt.ylabel("Normalized Intensity")
    plt.title("Cross-Sections of Normalized Intensity")
    plt.ylim(0, 1)
# legend for each frame
    plt.legend(np.arange(len(intensity_profiles)))
    plt.savefig('recorded/intensity_profiles_'+str(blob_width)+'.png')
    plt.show()


# intensities_normalized(noise_size, blob_width, intensity_profiles)

# plot intensity profiles mean
def mean_profile(intensity_profiles):
    mean_profile = np.mean(intensity_profiles, axis=0)
    normalized_mean_profile = (mean_profile - mean_profile.min()) / (mean_profile.max() - mean_profile.min())
    return normalized_mean_profile

# function to plot intensity profiles
def plot_intensity_profiles(intensity_profiles):
    normalized_mean_profile = mean_profile(intensity_profiles)
    plt.figure(figsize=(8, 6))
    # x label in degrees of visual angle
    plt.xticks(np.arange(0, noise_size[0], 100), np.arange(-5, 6, 1))
    plt.xlabel("Horizontal Position (degrees)")
    plt.ylabel("Normalized Intensity")
    plt.title("Cross-Sections of Normalized Intensity")
    plt.ylim(0, 1)
    # legend for each frame
    plt.legend(np.arange(len(intensity_profiles)))
    plt.savefig('recorded/intensity_profiles_'+str(blob_width)+'.png')
    plt.show()

# save intensity profiles as csv
np.savetxt("recorded/intensity_profiles.csv", intensity_profiles, delimiter=",")
# save intensity profiles as numpy array
np.save("recorded/intensity_profiles.npy", intensity_profiles)
# save mean intensity profile as numpy array
np.save("recorded/mean_intensity_profile.npy", mean_profile)