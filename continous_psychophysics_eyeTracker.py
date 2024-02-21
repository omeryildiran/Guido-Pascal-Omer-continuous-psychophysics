# optimizing continuous psychophysics experiment  code
# Omer F. Yildiran
# Start date: August 2023
# ENS-PSL LSP, Cognitice Science Masters 2nd year thesis project
# Supervisors: Guido Maiello and Pascal Mamassian


from psychopy import visual, core, event
import numpy as np
import psychopy.iohub as io
from psychopy.hardware import keyboard
import matplotlib.pyplot as plt
from PIL import Image
import os
from dva_to_pix import arcmin_to_px, dva_to_px
from numpy.random import choice as randchoice
from numpy.random import random, randint, normal, shuffle, choice as randchoice
from psychopy import sound, gui, visual, core, data, event, logging, clock, colors, layout
from psychopy.constants import (NOT_STARTED, STARTED, PLAYING, PAUSED,
                                STOPPED, FINISHED, PRESSED, RELEASED, FOREVER)
from psychopy.tools import monitorunittools
import sys  # to get file system encoding

from create_conditions import condition_creater
from psychopy import prefs
from audio_cue import create_stereo_sound, positional_audio
#from psychopy import microphone
from psychopy import visual, core, event

#from psychopy.sound import Sound

### Start eye tracker mouse
from psychopy.iohub import launchHubServer
sizeIs=512

win = visual.Window(size=(sizeIs,sizeIs), fullscr=False, monitor='testMonitor', units='pix', color=[0, 0, 0], useFBO=True,screen=1,colorSpace='rgb')

# --- Setup input devices ---
ioConfig = {}

# Setup eyetracking
ioConfig['eyetracker.hw.mouse.EyeTracker'] = {
    'name': 'tracker',
    'controls': {
        'move': [],
        'blink':('MIDDLE_BUTTON',),
        'saccade_threshold': 0.5,
    }
}

# Setup iohub keyboard
ioConfig['Keyboard'] = dict(use_keymap='psychopy')
ioServer = io.launchHubServer(window=win, **ioConfig)

tracker = ioServer.getDevice('tracker')

#########################################################################
prefs.hardware['audioLib'] = ['pygame']
"""          Experiment INFO Setup"""

# Store info about the experiment session
psychopyVersion = '2022.2.4'
expName = 'continous_psych'  # from the Builder filename that created this script
expInfo = {
    'participant': f"{randint(0, 999999):06.0f}",
    'session': '001',
}
# --- Show participant info dialog --
dlg = gui.DlgFromDict(dictionary=expInfo, sortKeys=False, title=expName)
if dlg.OK == False:
    core.quit()  # user pressed cancel
expInfo['date'] = data.getDateStr()  # add a simple timestamp
expInfo['expName'] = expName
expInfo['psychopyVersion'] = psychopyVersion

# Ensure that relative paths start from the same directory as this script
_thisDir = os.path.dirname(os.path.abspath(__file__))
os.chdir(_thisDir)
# Data file name stem = absolute path + name; later add .psyexp, .csv, .log, etc
filename = _thisDir + os.sep + u'data/%s_%s_%s' % (expInfo['participant'], expName, expInfo['date'])

conditions= condition_creater()
# import conditions.npy
#conditions=np.load('conditions.npy')
# Create PsychoPy window covering the whole screen
sizeIs=512
win = visual.Window(size=(sizeIs,sizeIs), fullscr=False, monitor='testMonitor', units='pix', color=[0, 0, 0], useFBO=True,screen=1,colorSpace='rgb')
field_size=[sizeIs,sizeIs]

frameRate=win.getActualFrameRate()
print(frameRate)
expInfo['frameRate']=frameRate
if expInfo['frameRate'] != None:
    frameDur = 1.0 / round(expInfo['frameRate'])
else:
    frameDur = 1.0 / 60.0  # could not measure, so guess
defaultKeyboard = keyboard.Keyboard(backend='iohub')
space2pass=keyboard.Keyboard()
endExpNow = False  # flag for 'escape' or other condition => quit the exp

#setup screen properties
screen_width=22 # actual size of my screen in cm is 28x17
screen_height=17
screen_distance=57
### Set the monitor to the correct distance and size
#win.monitor.setSizePix(field_size)
win.monitor.setWidth(screen_width)
win.monitor.setDistance(screen_distance)
mouse = event.Mouse(win=win,visible=False)
# TODO: arcmin to px conversion but open it later
# def arcmin_to_px(arcmin=1,h=19,d=57,r=1080):
#     dva= arcmin/60
#     return(monitorunittools.deg2pix(degrees=dva,monitor=win.monitor))
frameTolerance = 0.001  # how close to onset before 'same' frame



########################################################################################################################
"""                     NOISE BACKGROUND           """
# Background noise
# Noise properties
noise_size = field_size  # Use the window size for noise texture
noise_arcmin = 11  # Standard deviation for pixel noise || noise intensity Adjust to control noise intensity
noise_std =arcmin_to_px(noise_arcmin,h=screen_height,d=screen_distance,r=field_size[0])# convert arcmin to std for Gaussian
noiseQuantity=120
class NoiseGenerator:
    def __init__(self, noise_size, noise_intensity=1, noise_quantity=120):
        self.noise_size = noise_size
        self.noise_intensity = noise_intensity
        self.noise_quantity = noise_quantity
        self.noise_instances = self._generate_noises()

    def gauss_noise(self):
        noise = np.random.normal(0, self.noise_intensity, size=self.noise_size)
        return noise

    def _generate_noises(self):
        noise_instances = np.empty((self.noise_quantity, self.noise_size[0], self.noise_size[1]))
        for i in range(self.noise_quantity):
            noise_instances[i] = self.gauss_noise()
        
        noise_instances = np.clip(noise_instances, -3, 3)
        noise_instances = (noise_instances - noise_instances.min()) / (noise_instances.max() - noise_instances.min()) * 2 - 1
        return noise_instances

    def create_visual_objects(self, win):
        noise_objects = []
        for i in range(self.noise_quantity):
            noise_objects.append(visual.GratingStim(win, tex=self.noise_instances[i], size=self.noise_size, interpolate=False, units='pix'))
            noise_objects[-1].draw()  # draw noise
        return noise_objects

########################################################################################################################


"""             Welcome screen to give instructions to the participant         """
# Welcome screen
welcomeText="""
Welcome to the experiment your aim is to 
follow the blob with your mouse.
In the beginning of each trial you will see a fixation cross
Then you will see a blob moving around the screen
You will be asked to follow the blob with your mouse 
"""
welcomeText = visual.TextStim(win, text=welcomeText,
                              color=[1, 1, 1], units='pix', height=20)
welcomeText.draw()
win.flip()
noise_gen = NoiseGenerator(noise_size)
noise_obj = noise_gen.create_visual_objects(win)# Generate and draw noise objects
win.flip()
# draw a text that indicate the participant to press space to continue
pressText=visual.TextStim(win, text='Press any key to continue',
                        pos=(0, -300),
                        color='red',colorSpace='rgb',units='pix', height=20)
welcomeText.draw()
pressText.draw()
win.flip()

event.waitKeys()
win.flip()

# Brownian motion properties
blob_motion_std=2#arcmin_to_px(arcmin=1.32,h=screen_height,d=screen_distance,r=field_size[0])  # Standard deviation of Gaussian white noise velocities
# Create some handy timers
globalClock = core.Clock()  # to track the time since experiment started
routineTimer = core.Clock()  # to track time remaining of each (possibly non-sl€ip) routine 


####################### Blob Generation #######################
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
# create blobs for each of the conditions
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

"""                         BLOB MOTION                """ 
# Blob motion by changing velocity on x and y axis according to the brownian motion
def generateBrownianMotion(field_size, velocity_std, duration):
    # Generate random velocity for each frame
    num_frames = int(duration * expectedFrameRate)-1
    velocies_x=np.random.normal(0, velocity_std, num_frames)
    velocies_y=np.random.normal(0, velocity_std, num_frames)
    
    ## Create UP DOWN motion Calculate the number of frames for each motion
    # up_frames = num_frames // 3
    # down_frames = num_frames - up_frames
    # # Create arrays for each motion"
    # up_motion = np.ones(up_frames)
    # down_motion = -np.ones(down_frames)
    # velocies_y=np.concatenate((up_motion, down_motion))

    # new positions = old position + velocity
    pos_x = np.cumsum(velocies_x)
    pos_y = np.cumsum(velocies_y)
    # clip positions to stay within the field
    pos_x = np.clip(pos_x, -field_size/2, field_size/2 )
    pos_y = np.clip(pos_y, -field_size/2, field_size/2 )
    return (pos_x, pos_y)
    
##            OBSERVER POINTER             """
obs_pointer = visual.Circle(win, radius=10, fillColor='red',colorSpace='rgb', units='pix',size=0.30)
##    fixation cross for the beginning of the trial (before start of the trial )  
fixationCross=visual.TextStim(win, text='+', color=[1, 1, 1], units='pix', height=20)

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


"""                  Have a REST SCREEN       """
trialNum=1
haveRest=False
haveRestText=visual.TextStim(win, text='Press space to continue', color=[1, 1, 1], units='pix', height=20)
haveRestNum=1
#####
sigma_trials=[]
win.setMouseVisible(False)      
# Check for and print any eye tracker events received...
tracker.setRecordingState(True)  
for blob_width in conditions:
    sigma_trials.append(blob_width)
    _space2pass_allKeys = []
    space2pass.keys = []
    space2pass.clearEvents(eventType='keyboard')
    #save conditions
    if trialNum%haveRestNum==0 and trialNum!=1:
        haveRest=True
        while haveRest:
            haveRestText.draw()
            win.flip()
            theseKeys = space2pass.getKeys(keyList=['space'], waitRelease=False)
            _space2pass_allKeys.extend(theseKeys)
            if len(_space2pass_allKeys)>0:
                haveRest=False

    mouse.setPos((0,0))
    # set a timer for the next line record the time spent
    tStart = globalClock.getTime()
    # clear all drawings
    #print(blob_width)
    blob_std=arcmin_to_px(arcmin=blob_width,h=screen_height,d=screen_distance,r=field_size[0])
    # create blob
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
    pos_x_obj, pos_y_obj = generateBrownianMotion(field_size=noise_size[0], velocity_std=blob_motion_std, duration=expectedDuration)# pre-generate brownian motion
    pos_x_obj=np.insert(pos_x_obj,0,0)
    pos_y_obj=np.insert(pos_y_obj,0,0)
    # calcuate velocity of blob
    blob_v = np.sqrt(np.diff(pos_x_obj)**2 + np.diff(pos_y_obj)**2)
    blob_v = np.insert(blob_v, 0, 0)
    # save blob velocity
    all_blob_x.append(pos_x_obj)
    all_blob_y.append(pos_y_obj)
    all_blob_v.append(blob_v)

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
        fixationCross.draw()
        win.flip()
    # draw the stimulus
    blob_obj.setAutoDraw(True)
    obs_pointer.setAutoDraw(True)
    mouse.setVisible(False)  
    response_x=np.empty((expectedFrames))
    response_y=np.empty((expectedFrames))
    mouse_v=np.empty((expectedFrames))
    checkPointX=0
    checkPointY=0

    # Define the colors two shades of green
    green_light = np.array([0.5, 1, 0.5])  # Light green
    green_dark = np.array([0, 0.5, 0])  # Dark green
    red = np.array([1, 0, 0])  # Red
    green = np.array([0, 1, 0])  # Green
    curColor=[1, 1, 1]

    tStart = globalClock.getTime()
    """ TRIAL STARTS HERE """
    """ This is where trial starts"""
    ##################### Trial Start #####################
    while continueRoutine:
        random_noise_index = int(randint(0, noiseQuantity))
        # t = routineTimer.getTime()
        #tThisFlip = win.getFutureFlipTime(clock=routineTimer)
        #tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN+= 1  # number of completed frames (so 0 is the first frame)


        # draw the stimulus
        obs_pointer.setPos(mouse.getPos())
        noise_obj[random_noise_index].draw()    # draw noise
        blob_obj.setPos((pos_x_obj[frameN], pos_y_obj[frameN]))
        # Eye tracker recording
        print(tracker.getPosition())
    
        #save mouse position
        response_x[frameN]=obs_pointer.pos[0]
        response_y[frameN]=obs_pointer.pos[1]
        # flip the window
        win.flip()
        # end the loop after given seconds
        if frameN > expectedFrames-2:
            continueRoutine = False
        # check for quit (typically the Esc key)
        if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
            endExpNow = True
            continueRoutine = False
            #core.quit()

    """ TRIAL ENDS HERE """

    tEnd = globalClock.getTime()
    print("trial lasteD: "+str(tEnd-tStart))
    blob_obj.setAutoDraw(False)
    obs_pointer.setAutoDraw(False)
    trialNum+=1
    

    t = routineTimer.getTime()
    
    print(t)
    print(frameN)

    # save positions of target, mouse, and but first calculate velocities
    response_x = np.array(response_x)
    response_y = np.array(response_y)
    mouse_v = np.sqrt(np.diff(response_x)**2 + np.diff(response_y)**2)
    mouse_v = np.insert(mouse_v, 0, 0)
    all_mouse_x.append(response_x)
    all_mouse_y.append(response_y)
    all_mouse_v.append(mouse_v)
    space2pass.keys = None

    if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
        win.close()
        break

# Create an Ending Screen
endText=visual.TextStim(win, text='Thank you for your participation', color=[1, 1, 1], units='pix', height=20)
endText.draw()
win.flip()
# and wait for participant to press space
event.waitKeys()
win.close()

##

# Close the window
#event.waitKeys()
#win.close()
# write mouse positions, velocities, and blob velocities and positions and conditions as a matlab data file
import scipy.io as sio
conditions=conditions.T
filename = u'data/%s_%s_%s' % (expInfo['participant'], expName, expInfo['date'])
sio.savemat(filename+'.mat', {'mouse_x': all_mouse_x, 'mouse_y': all_mouse_y, 'response': all_mouse_v, 'blob_x': all_blob_x, 'blob_y': all_blob_y, 'target': all_blob_v, 'sigma': sigma_trials})
core.quit()



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