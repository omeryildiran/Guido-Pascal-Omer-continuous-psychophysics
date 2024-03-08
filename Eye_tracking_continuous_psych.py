# optimizing continuous psychophysics experiment  code
# Omer F. Yildiran
# Start date: August 2023
# ENS-PSL LSP, Cognitice Science Masters 2nd year thesis project
# Supervisors: Guido Maiello and Pascal Mamassian


from psychopy import visual, core, event
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import os
from dva_to_pix import arcmin_to_px, dva_to_px,arcmin_to_dva
from numpy.random import choice as randchoice
from numpy.random import random, randint, normal, shuffle, choice as randchoice
from psychopy import sound, gui, visual, core, data, event, logging, clock, colors, layout, iohub, hardware

from psychopy.constants import (NOT_STARTED, STARTED, PLAYING, PAUSED,
                                STOPPED, FINISHED, PRESSED, RELEASED, FOREVER)
import sys  # to get file system encoding

from create_conditions import condition_creater, create_conditions
from psychopy import prefs
from audio_cue import create_stereo_sound, positional_audio
#from psychopy import microphone
from psychopy import visual, core, event
#from psychopy.sound import Sound
prefs.hardware['audioLib'] = ['pygame']
from psychopy.hardware import keyboard
import psychopy.iohub as io
from psychopy.iohub.util import hideWindow, showWindow
from psychopy.tools.monitorunittools import deg2pix, pix2deg
from psychopy import monitors

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

conditions= create_conditions(numOfBlocks=1, blob_widths=[11,17,25], repeats=5)
# import conditions.npy
#conditions=np.load('conditions.npy')
# Create PsychoPy window covering the whole screen
sizeIs=1024
#setup screen properties
screen_width=37 #31 asuSs 14 # actual size of my screen in cm is 28x17
screen_height=23 # 16.5 asus
screen_distance=57
# define monitor
labMonitor=monitors.Monitor('labMon', width=37, distance=57)
labMonitor.setSizePix((sizeIs, sizeIs))
myMon=monitors.Monitor('asusMon', width=31, distance=57)
#labMonitor.setgamma(2.4)
#labMonitor.setwhite([0.95, 0.95, 0.95])  # Example white point with slightly lower values
win = visual.Window(size=(sizeIs,sizeIs),
                     fullscr=True, 
                     monitor=labMonitor, 
                    units='pix', 
                    color=[0, 0, 0],
                      useFBO=True,
                      screen=0,
                      colorSpace='rgb')
field_size=[sizeIs,sizeIs]
### Set the monitor to the correct distance and size
#win.monitor.setSizePix(field_size)
win.monitor.setWidth(screen_width)
win.monitor.setDistance(screen_distance)

frameRate=win.getActualFrameRate()
print(frameRate)
expInfo['frameRate']=frameRate
if expInfo['frameRate'] != None:
    frameDur = 1.0 / round(expInfo['frameRate'])
else:
    frameDur = 1.0 / 60.0  # could not measure, so guess
space2pass=keyboard.Keyboard()
endExpNow = False  # flag for 'escape' or other condition => quit the exp


mouse = event.Mouse(win=win,visible=False)
# TODO: arcmin to px conversion but open it later
# def arcmin_to_px(arcmin=1,h=19,d=57,r=1080):
#     dva= arcmin/60
#     return(monitorunittools.deg2pix(degrees=dva,monitor=win.monitor))
frameTolerance = 0.001  # how close to onset before 'same' frame



# Eye tracker to use ('mouse', 'eyelink', 'gazepoint', or 'tobii')
TRACKER = 'eyelink'
BACKGROUND_COLOR = [0, 0, 0]

devices_config = dict()
eyetracker_config = dict(name='tracker')
if TRACKER == 'mouse':
    eyetracker_config['calibration'] = dict(screen_background_color=BACKGROUND_COLOR,
                                            auto_pace=True,
                                            target_attributes=dict(animate=dict(enable=True, expansion_ratio=1.5,
                                                                                contract_only=False))
                                            )
    devices_config['eyetracker.hw.mouse.EyeTracker'] = eyetracker_config
elif TRACKER == 'eyelink':
    eyetracker_config['model_name'] = 'EYELINK 1000 DESKTOP'
    eyetracker_config['simulation_mode'] = False
    eyetracker_config['runtime_settings'] = dict(sampling_rate=1000, track_eyes='RIGHT')
    eyetracker_config['calibration'] = dict(screen_background_color=BACKGROUND_COLOR, auto_pace=True)
    devices_config['eyetracker.hw.sr_research.eyelink.EyeTracker'] = eyetracker_config
else:
    print("{} is not a valid TRACKER name; please use 'mouse', 'eyelink', 'gazepoint', or 'tobii'.".format(TRACKER))
    core.quit()

# # --- Setup input devices ---
#ioConfig = {}

# # Setup eyetracking
# ioConfig['eyetracker.hw.sr_research.eyelink.EyeTracker'] = {
#     'name': 'tracker',
#     'model_name': 'EYELINK 1000 DESKTOP',
#     'simulation_mode': False,
#     'network_settings': '100.1.1.1',
#     'default_native_data_file_name': 'EXPFILE',
#     'runtime_settings': {
#         'sampling_rate': 1000.0,
#         'track_eyes': 'LEFT',
#         'sample_filtering': {
#             'sample_filtering': 'FILTER_LEVEL_2',
#             'elLiveFiltering': 'FILTER_LEVEL_OFF',
#         },
#         'vog_settings': {
#             'pupil_measure_types': 'PUPIL_AREA',
#             'tracking_mode': 'PUPIL_CR_TRACKING',
#             'pupil_center_algorithm': 'ELLIPSE_FIT',
#         }
#     }
# }

# Setup iohub keyboard
#ioConfig['Keyboard'] = dict(use_keymap='psychopy')

ioSession = '1'
if 'session' in expInfo:
    ioSession = str(expInfo['session'])

#ioServer = io.launchHubServer(window=win, **ioConfig)
ioServer = io.launchHubServer(window=win, **devices_config)

eyetracker = ioServer.getDevice('tracker')
# create a default keyboard (e.g. to check for escape)
defaultKeyboard = keyboard.Keyboard(backend='iohub')
# Display calibration gfx window and run calibration.
result = eyetracker.runSetupProcedure()
print("Calibration returned: ", result)




########################################################################################################################
### -------Initialize componens  for target tracking trial-------------------------
"""                 NOISE BACKGROUND           """
# Background noise
# Noise properties
noise_size = field_size  # Use the window size for noise texture
noise_arcmin = 11  # Standard deviation for pixel noise || noise intensity Adjust to control noise intensity
noise_std =deg2pix(degrees=arcmin_to_dva(noise_arcmin),monitor=win.monitor) #arcmin_to_px(noise_arcmin,h=screen_height,d=screen_distance,r=field_size[0])# convert arcmin to std for Gaussian
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
blob_widths=list(set(conditions))#[11,13,17,21,25,29]
blobs=[]


for blob_width in blob_widths:
    blobs.append(np.array(generateBlob(arcmin_to_px(arcmin=blob_width,h=screen_height,d=screen_distance,r=field_size[0]))))
blobs=np.array(blobs)
blobs = (blobs - blobs.min()) / (blobs.max() - blobs.min())*2-1
# create dictionary for each blob_width
blob_dict={}
for i in range(len(blob_widths)):
    blob_dict[blob_widths[i]]=blobs[i]


expectedFrameRate=60
expectedDuration=20 # in seconds
expectedFrames=expectedFrameRate*expectedDuration

"""                         BLOB MOTION                """ 
# Blob motion by changing velocity on x and y axis according to the brownian motion
def generateBrownianMotion(field_size, velocity_std, duration):
    # Generate random velocity for each frame
    num_frames = int(duration * expectedFrameRate)+1500 # add 600 frames to the duration
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

# data to save
## blob data
all_blob_x = []
all_blob_y = []
all_blob_v=[]
### mouse data
all_mouse_x = []
all_mouse_y = []
all_mouse_v=[]
### stimulus data
all_stim_x = []
all_stim_y = []
all_stim_v=[]
jumped=False


## --- define have a rest screen ###
trialNum=1
haveRest=False
haveRestText=visual.TextStim(win, text='Press space to continue', color=[1, 1, 1], units='pix', height=20)
haveRestNum=1
#####
sigma_trials=[]
win.setMouseVisible(False)    

# create a warning text to indicate the participant to look at the blob if they look 2dva away from the blob
look_warning=visual.TextStim(win, text='Look at the blob', color=[1, 1, 1], units='pix', height=20)

for blob_width in sorted(conditions):
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
    realFrameN=-1

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
    obs_pointer.setAutoDraw(False)
    mouse.setVisible(False)  
    response_x=np.empty((expectedFrames))
    response_y=np.empty((expectedFrames))
    stim_x=np.empty((expectedFrames))
    stim_y=np.empty((expectedFrames))
    mouse_v=np.empty((expectedFrames))
    checkPointX=0
    checkPointY=0

    #ioServer.ClearEvents()
    tStart = globalClock.getTime()
    ##################### Trial Start #####################
    eyetracker.setRecordingState(True)  # start recording of gaze
    while continueRoutine:
        realFrameN=realFrameN+1
        random_noise_index = int(randint(0, noiseQuantity))
        #frameN+= 1  # number of completed frames (so 0 is the first frame)

        # Get the latest gaze position in display coord space.
        gpos = eyetracker.getLastGazePosition()
        valid_gaze_pos = isinstance(gpos, (tuple, list))
        if valid_gaze_pos:
            mX=gpos[0]
            mY=gpos[1]
            jumped=np.sqrt((mX-obs_pointer.pos[0])**2+(mY-obs_pointer.pos[1])**2)>blob_motion_std*20
            looksAway=(np.sqrt((pos_x_obj[realFrameN-1]-mX)**2+(pos_y_obj[realFrameN-1]-mY)**2))>blob_motion_std*30
        if jumped or looksAway or not valid_gaze_pos:
            pass
        else:
            frameN+=1
        # Update stim based on gaze position
        # draw the stimulus
        if valid_gaze_pos and gpos is not None:
            obs_pointer.setPos(gpos)
        noise_obj[random_noise_index].draw()    # draw noise
        blob_obj.setPos((pos_x_obj[realFrameN], pos_y_obj[realFrameN]))

        #save mouse position
        if not jumped or not looksAway:
            stim_x[frameN]=blob_obj.pos[0]
            stim_y[frameN]=blob_obj.pos[1]
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

    # the Routine "trial" was not non-slip safe, so reset the non-slip timer
    eyetracker.setRecordingState(False)  # stop recording of gaze
    tEnd = globalClock.getTime()
    print("trial lasteD: "+str(tEnd-tStart))
    blob_obj.setAutoDraw(False)
    obs_pointer.setAutoDraw(False)
    trialNum+=1
    t = routineTimer.getTime()
    print("time of finish: "+str(t))
    print("Maximum frame achieved"+str(frameN))
    print("Maximum frame achieved"+str(realFrameN))

    # save positions of target, mouse, and but first calculate velocities
    response_x = np.array(response_x)
    response_y = np.array(response_y)
    mouse_v = np.sqrt(np.diff(response_x)**2 + np.diff(response_y)**2)
    mouse_v = np.insert(mouse_v, 0, 0)
    all_mouse_x.append(response_x)
    all_mouse_y.append(response_y)
    all_mouse_v.append(mouse_v)
    stim_x = np.array(stim_x)
    stim_y = np.array(stim_y)
    stim_v = np.sqrt(np.diff(stim_x)**2 + np.diff(stim_y)**2)
    stim_v = np.insert(stim_v, 0, 0)
    # save mouse velocity
    all_stim_x.append(stim_x)
    all_stim_y.append(stim_y)
    all_stim_v.append(stim_v)
    space2pass.keys = None

    if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
        # Create an Ending Screen
        endText=visual.TextStim(win, text='Thank you for your participation', color=[1, 1, 1], units='pix', height=20)
        endText.draw()
        win.flip()
        # and wait for participant to press space
        event.waitKeys()
        win.close()
        eyetracker.setConnectionState(False)
        break

# Create an Ending Screen
if not endExpNow:
    endText=visual.TextStim(win, text='Thank you for your participation press any key to exit', color=[1, 1, 1], units='pix', height=20)
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
# make conditions as numpy
conditions=np.array(conditions)
conditions=conditions.T
filename = u'data/%s_%s_%s' % (expInfo['participant'], expName, expInfo['date'])
sio.savemat(filename+'.mat', {'mouse_x': all_mouse_x, 'mouse_y': all_mouse_y, 'response': all_mouse_v, 'blob_x': all_stim_x, 'blob_y': all_stim_y, 'target': all_stim_v, 'sigma': sigma_trials})
core.quit()
