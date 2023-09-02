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

# Create PsychoPy window covering the whole screen
win = visual.Window(size=(512, 512), fullscr=False, monitor='testMonitor', units='pix', color=[0, 0, 0], useFBO=True)
field_size=(512,512)

#setup screen properties
screen_width=31
screen_height=17.5
screen_distance=50
### Set the monitor to the correct distance and size
#win.monitor.setSizePix(field_size)
# win.mouseVisible = False
win.monitor.setWidth(screen_width)
win.monitor.setDistance(screen_distance)

# Store info about the experiment session
psychopyVersion = '2022.2.4'
expName = 'expectation_shapes_perceived_time'  # from the Builder filename that created this script
expInfo = {
    'participant': f"{randint(0, 999999):06.0f}",
    'session': '001',
}
expInfo['date'] = data.getDateStr()  # add a simple timestamp
expInfo['expName'] = expName
expInfo['psychopyVersion'] = psychopyVersion

frameRate=win.getActualFrameRate()
expInfo['frameRate']=frameRate
if expInfo['frameRate'] != None:
    frameDur = 1.0 / round(expInfo['frameRate'])
else:
    frameDur = 1.0 / 60.0  # could not measure, so guess
defaultKeyboard = keyboard.Keyboard(backend='iohub')
endExpNow = False  # flag for 'escape' or other condition => quit the exp


# Noise properties
noise_size = field_size  # Use the window size for noise texture
noise_arcmin = 11  # Standard deviation for pixel noise || noise intensity Adjust to control noise intensity
noise_std =arcmin_to_px(noise_arcmin,h=screen_height,d=screen_distance,r=field_size[0])# convert arcmin to std for Gaussian

# Brownian motion properties
velocity_std = 1.0  # Standard deviation of Gaussian white noise velocities

# Create some handy timers
globalClock = core.Clock()  # to track the time since experiment started
routineTimer = core.Clock()  # to track time remaining of each (possibly non-slip) routine 

# Background noise
def gaussNoise(noise_intensity=1):
    noise = np.random.normal(0, noise_intensity, size=noise_size)
    noise = np.clip(noise, -3*noise_intensity, 3*noise_intensity)    
    noise = (noise - noise.min()) / (noise.max() - noise.min())*2-1 # normalize the noise
    return noise

# Blob properties
blob_width=11 # in arcmins
initial_blob_std= arcmin_to_px(arcmin=11,h=screen_height,d=screen_distance,r=field_size[0]) 
# blob properties conversion
#blob_std = space_constant / (2 * np.sqrt(2 * np.log(2))) # Convert arcmin to std for Gaussian
blob_std=arcmin_to_px(arcmin=blob_width,h=screen_height,d=screen_distance,r=field_size[0])
total_lum=1*(2*np.pi*initial_blob_std ** 2)*1 # Total luminance of blobm
blob_amplitude =total_lum / (2*np.pi*blob_std ** 2) # Adjust blob amplitude to keep total blob energy constant

# Blob generator
def generateBlob(space_constant):
    y, x = np.meshgrid(np.arange(noise_size[1]), np.arange(noise_size[0]))
    blob = blob_amplitude * np.exp(-((x - noise_size[0] / 2)**2 + (y - noise_size[1] / 2)**2) / (2 * blob_std**2))
    return blob

# Blob motion by changing velocity on x and y axis according to the brownian motion
def generateBrownianMotion(field_size, velocity_std, duration):
    # Generate random velocity for each frame
    num_frames = int(duration * frameRate)+10
    velocies_x=np.random.normal(0, velocity_std, num_frames)
    velocies_y=np.random.normal(0, velocity_std, num_frames)
    # new positions = old position + velocity
    pos_x = np.cumsum(velocies_x)
    pos_y = np.cumsum(velocies_y)
    # clip positions to stay within the field
    pos_x = np.clip(pos_x, -field_size[0]/2, field_size[0]/2 )
    pos_y = np.clip(pos_y, -field_size[0]/2, field_size[0]/2 )
    return (pos_x, pos_y)

print(frameRate)

""" Pregenerate noise and blob instances for each frame """
duration=10 # in seconds
noise_instances=[]
noise_instances.append(gaussNoise())# first noise
blob=generateBlob(blob_width)    # Create Gaussian blob
pos_x, pos_y = generateBrownianMotion(field_size=noise_size[0], velocity_std=1, duration=duration)# pre-generate brownian motion
blob_rolleds=[]
blob_rolleds.append(blob) # first blob instance

## create visual patches for each frame
final_stims=[]
for i in range(int(duration*frameRate)+10): #assuming program can achieve the actual frame rate of screen
    noise_instances.append(gaussNoise())# pre-generate noise

    blob_rolleds.append(np.roll(blob, (int(pos_x[(i)]), int(pos_y[(i)])), axis=(1, 0))) # Shift blob according to the brownian motion

    stim_array=blob_rolleds[i]+noise_instances[i]
    #clip the stim array to stay within the range of 3*noise_intensity
    stim_array = np.clip(stim_array, -3*noise_std, 3*noise_std)
    stim_array = (stim_array - stim_array.min()) / (stim_array.max() - stim_array.min())*2-1 # normalize the stim array
    final_stim = visual.PatchStim(win, tex=stim_array, size=noise_size, interpolate=False,units='pix')
    final_stims.append(final_stim)
    final_stims[-1].draw()  # first draw is slower. So do it now.


##################### Loop Start #####################
intensity_profiles = []  # List to store intensity profiles
continueRoutine = True
_timeToFirstFrame = win.getFutureFlipTime(clock="now")
frameN = -1
win.clearBuffer()
routineTimer.reset()
t = 0
while continueRoutine:
    t = routineTimer.getTime()
    tThisFlip = win.getFutureFlipTime(clock=routineTimer)
    tThisFlipGlobal = win.getFutureFlipTime(clock=None)
    frameN+= 1  # number of completed frames (so 0 is the first frame)

    # show frame number on top right corner
    #frame_text = visual.TextStim(win, text="Frame: "+str(frameN), pos=(field_size[0]/2-150,field_size[1]/2-100 ), color=[1,0,0], height=50, alignHoriz='center', alignVert='center', units='pix')
    #frame_text.draw()

    # draw the stimulus
    final_stims[frameN].draw()
    print(frameN)
    #print("time="+str(t))

    # flip the window
    win.flip(clearBuffer=True)




    # end the loop after given seconds
    if t > duration:
        continueRoutine = False
    # check for quit (typically the Esc key)
    if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
        core.quit()
    
    
    # save the intensity profile for each frame
    # stim_array = np.array(final_stims[frameN].tex)
    # intensity_profile = stim_array[noise_size[1]//2,:]
    # morm_intensity_profile = stim_array[noise_size[1]//2,:]
    # # normalize the intensity profile
    # normalized_profile = (intensity_profile - intensity_profile.min()) / (intensity_profile.max() - intensity_profile.min())
    # intensity_profiles.append(normalized_profile)
    # #intensity_profiles.append(intensity_profile)
    # # plot the intensity profile for each frame
    # #plt.plot(intensity_profiles)
    # plt.plot(intensity_profiles[-1]) 
    # plt.xlabel("Horizontal Position (pixels)")
    # plt.ylabel("Intensity")
    # plt.title("Cross-Sections of Intensity"+ "  Frame: " + str(frameN))
    # plt.ylim(0, 1)
    # plt.pause(0.0001)
    # if frameN == 0:
    #     plt.savefig('recorded/intensity_profile_'+str(blob_width)+'.png')
    # plt.clf()
# Close the window
event.waitKeys()
win.close()

# save intensity profiles based on final_stims array
intensity_profiles = []  # List to store intensity profiles
for frame in final_stims:
    stim_array = np.array(frame.tex)
    intensity_profile = stim_array[noise_size[1]//2,:]
    intensity_profiles.append(intensity_profile)
# plot the intensity profile for each frame
plt.plot(intensity_profiles[-1])
plt.xlabel("Horizontal Position (pixels)")
plt.ylabel("Intensity")
plt.title("Cross-Sections of Intensity"+ "  Frame: " + str(frameN))
plt.ylim(0, 1)
plt.show()
# if frameN == 0:
#     plt.savefig('recorded/intensity_profile_'+str(blob_width)+'.png')
# plt.clf()



# Plot intensity profiles
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

# plot intensity profiles mean
plt.figure(figsize=(8, 6))
mean_profile = np.mean(intensity_profiles, axis=0)
normalized_mean_profile = (mean_profile - mean_profile.min()) / (mean_profile.max() - mean_profile.min())
plt.plot(mean_profile)
# x label in degrees of visual angle
plt.xticks(np.arange(0, noise_size[0], 100), np.arange(-5, 6, 1))
plt.xlabel("Horizontal Position (degrees)")
plt.ylabel("Intensity")
plt.title("Mean of Intensity")
plt.ylim(0, 1)
plt.savefig('recorded/mean_intensity_profile_'+str(blob_width)+'.png')
plt.show()

# save intensity profiles as csv
np.savetxt("recorded/intensity_profiles.csv", intensity_profiles, delimiter=",")
# save intensity profiles as numpy array
np.save("recorded/intensity_profiles.npy", intensity_profiles)
# save mean intensity profile as numpy array
np.save("recorded/mean_intensity_profile.npy", mean_profile)