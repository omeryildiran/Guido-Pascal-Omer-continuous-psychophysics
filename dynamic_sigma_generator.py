
import numpy as np
from psychopy import visual
""" --------------- SIGMA DRIFT ----------------------         """
# Constants
screen_height = 1024
screen_distance = 57
field_size = (screen_height, screen_height)
expectedFrameRate = 60
toleranceTrialN = 300
expectedDuration=20 # in seconds
expectedFrames=expectedFrameRate*expectedDuration

class SigmaGenerator:
    def __init__(self, meanSigma=15, velocity_std=1, duration=20, minSigma=5, maxSigma=35):
        self.meanSigma = meanSigma
        self.velocity_std = velocity_std
        self.duration = duration
        self.minSigma = minSigma
        self.maxSigma = maxSigma
        self.sigma_width = self._generate_walk()

    
    def _generateBlob(blob_width):
        blob_std=arcmin_to_px(arcmin=blob_width,h=screen_height,d=screen_distance,r=field_size[0])
        blob_amplitude =total_lum / (2*np.pi*blob_std ** 2) # Adjust blob amplitude to keep total blob energy constant
        y, x = np.meshgrid(np.arange(noise_size[1]), np.arange(noise_size[0]))
        blob = blob_amplitude * np.exp(-((x - noise_size[0] / 2)**2 + (y - noise_size[1] / 2)**2) / (2* blob_std**2))
        #blob = (blob - blob.min()) / (blob.max() - blob.min())*2-1 # normalize the blob
        return blob


    def _generate_walk(self):
        num_frames = int(self.duration * expectedFrameRate) + toleranceTrialN
        increments = np.array([-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1])
        sigma_width = np.zeros(num_frames)
        sigma_width[0] = self.meanSigma  # Start from the meanSigma

        for i in range(1, num_frames):
            # Filter increments to only those that do not push sigma_width out of bounds
            valid_increments = increments[(sigma_width[i-1] + increments >= self.minSigma) & (sigma_width[i-1] + increments <= self.maxSigma)]
            if valid_increments.size == 0:
                # If no valid increments, maintain current position
                sigma_width[i] = sigma_width[i-1]
            else:
                # Randomly choose from the valid increments
                sigma_width[i] = sigma_width[i-1] + np.random.choice(valid_increments)
        return sigma_width

    def get_sigma_width(self):
        return self.sigma_width

    def get_sigma_width_at(self, frame):
        return self.sigma_width[frame]

    def get_sigma_width_at_time(self, time):
        return self.sigma_width[int(time * expectedFrameRate)]

    def get_sigma_width_at_time(self, time):
        return self.sigma_width[int(time * expectedFrameRate)]

    def get_sigma_width_at_frame(self, frame):
        return self.sigma_width[frame]

    def get_sigma_width_at_frame(self, frame):
        return self.sigma_width[frame]

    def get_sigma_width_at_time(self, time):
        return self.sigma_width[int(time * expectedFrameRate)]

    def get_sigma_width_at_time(self, time):
        return self.sigma_width[int(time * expectedFrameRate)]

    def get_sigma_width_at_frame(self, frame):
        return self.sigma_width[frame]

    def create_visual_blobs(self, win):
        blobs = []
        for sigma in self.sigma_width:
            blobs.append(self._generate_blob(sigma))
        return blobs
        

def generateWalkofSigmaDifficulty(meanSigma=15, velocity_std=1, duration=20, minSigma=5, maxSigma=35):
    num_frames = int(duration * expectedFrameRate) + toleranceTrialN
    increments = np.array([-1, -0.75, -0.5, -0.25, 0,0.25, 0.5, 0.75, 1])
    print(increments.std())
    sigma_width = np.zeros(num_frames)
    sigma_width[0] = meanSigma  # Start from the meanSigma

    for i in range(1, num_frames):
        # Filter increments to only those that do not push sigma_width out of bounds
        valid_increments = increments[(sigma_width[i-1] + increments >= minSigma) & (sigma_width[i-1] + increments <= maxSigma)]
        if valid_increments.size == 0:
            # If no valid increments, maintain current position
            sigma_width[i] = sigma_width[i-1]
        else:
            # Randomly choose from the valid increments
            sigma_width[i] = sigma_width[i-1] + np.random.choice(valid_increments)
    return sigma_width

# Generate sigma drift for all possible sigmas in the experiment
minBlobSigma=10
maxBlobSigma=28
initial_blob_std= arcmin_to_px(arcmin=minBlobSigma,h=screen_height,d=screen_distance,r=field_size[0]) 
total_lum=200*(2*np.pi*initial_blob_std ** 2)*1 # Total luminance of blobm

blobs=[]
allPossibleSigma=np.arange(minBlobSigma, maxBlobSigma+0.25, 0.25)
for sigma in allPossibleSigma:
    blobs.append(np.array(generateBlob(arcmin_to_px(arcmin=sigma,h=screen_height,d=screen_distance,r=field_size[0]))))
blobs=np.array(blobs)
blobs = (blobs - blobs.min()) / (blobs.max() - blobs.min())*2-1

blobVisualObjsDict={}
for i,sigma in enumerate(allPossibleSigma):
    blobMask=blobs[i]
    blobVisualObjsDict[sigma]=visual.GratingStim(win, tex=None, mask=blobs[i],  units='pix', interpolate=False)

