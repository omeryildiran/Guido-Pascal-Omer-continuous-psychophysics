from psychopy import visual, core, event
import numpy as np
from psychopy.visual import filters
import psychopy.iohub as io
from psychopy.hardware import keyboard
import matplotlib.pyplot as plt
from PIL import Image
import os
# Create PsychoPy window covering the whole screen
win = visual.Window(size=(1024, 1024), fullscr=False, monitor='testMonitor', units='pix', color=[0, 0, 0])#, useFBO=True)



defaultKeyboard = keyboard.Keyboard(backend='iohub')
endExpNow = False  # flag for 'escape' or other condition => quit the exp

# Noise properties
noise_duration = 10  # Duration of noise presentation in seconds
noise_size = win.size  # Use the window size for noise texture

noise_arcmin = 11  # Standard deviation for pixel noise || noise intensity Adjust to control noise intensity
# convert arcmin to std for Gaussian
noise_std = noise_arcmin / (2 * np.sqrt(2 * np.log(2)))


# Create some handy timers
globalClock = core.Clock()  # to track the time since experiment started
routineTimer = core.Clock()  # to track time remaining of each (possibly non-slip) routine 

# Background noise
def gaussNoise(noise_intensity=noise_std):
    noise = np.random.normal(0, 3, size=noise_size)
    noise = np.clip(noise, -3*noise_std, 3 * noise_std)    
    noise = (noise - noise.min()) / (noise.max() - noise.min()) * 2 - 1
    #noise_stim = visual.ImageStim(win, image=noise, size=noise_size, units='pix')
    return noise


initial_blob_std= 11 / (2 * np.sqrt(2 * np.log(2)))  # Convert arcmin to std for Gaussian
total_lum=1*(2*np.pi*initial_blob_std ** 2)*1 # Total luminance of blob
# Blob generator
def generateBlob(space_constant):
    y, x = np.meshgrid(np.arange(noise_size[1]), np.arange(noise_size[0]))

    # blob properties conversion
    blob_std = space_constant / (2 * np.sqrt(2 * np.log(2)))
    blob_amplitude =total_lum / (2*np.pi*blob_std ** 2) # Adjust blob amplitude to keep total blob energy constant
    blob = blob_amplitude * np.exp(-((x - noise_size[1] / 2)**2 + (y - noise_size[0] / 2)**2) / (2 * blob_std**2))
    return blob



intensity_profiles = []  # List to store intensity profiles
##################### Loop Start #####################
continueRoutine = True
t = 0
_timeToFirstFrame = win.getFutureFlipTime(clock="now")
frameN = -1
blob_width=11 # in arcmins
while continueRoutine:
    t = routineTimer.getTime()
    tThisFlip = win.getFutureFlipTime(clock=routineTimer)
    tThisFlipGlobal = win.getFutureFlipTime(clock=None)
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)

    noise=gaussNoise() # noise instance
    blob=generateBlob(blob_width)    # Create Gaussian blob

    # combine blob and noise
    stim_array=blob + noise
    # create a new psychopy image with the blob and noise
    stim = visual.ImageStim(win, image=stim_array, size=noise_size, units='pix')
    # draw the stim
    stim.draw()
    # flip the window
    win.flip()

    # save the intensity profile for each frame
    intensity_profile = stim_array[noise_size[1]//2,:]
    # normalize the intensity profile
    normalized_profile = (intensity_profile - intensity_profile.min()) / (intensity_profile.max() - intensity_profile.min())
    intensity_profiles.append(normalized_profile)
    #intensity_profiles.append(intensity_profile)
    # plot the intensity profile for each frame
    #plt.plot(intensity_profiles)
    plt.plot(intensity_profiles[-1]) 
    plt.xlabel("Horizontal Position (pixels)")
    plt.ylabel("Intensity")
    plt.title("Cross-Sections of Intensity"+ "  Frame: " + str(frameN))
    plt.ylim(0, 1)
    plt.pause(0.0001)
    if frameN == 0:
        plt.savefig('recorded/intensity_profile_'+str(blob_width)+'.png')
    plt.clf()

    # end the loop after given seconds
    if t > 10:
        continueRoutine = False
    # check for quit (typically the Esc key)
    if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
        core.quit()
    
# Close the window
event.waitKeys()
win.close()



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