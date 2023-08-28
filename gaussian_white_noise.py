from psychopy import filters
import numpy as np
from psychopy import visual, event

def makeFilteredNoise(res, radius, shape='gauss'):
    noise = np.random.random([res, res])
    kernel = filters.makeMask(res, shape=shape, radius=radius)
    filteredNoise = filters.conv2d(kernel, noise)
    filteredNoise = (filteredNoise-filteredNoise.min())/(filteredNoise.max()-filteredNoise.min())*2-1
    return filteredNoise
    
def arcmin2deg(arcmin):
    return arcmin/60.0

radius = arcmin2deg(.05)


filteredNoise = makeFilteredNoise(256, radius)

win = visual.Window([800,800], monitor='testMonitor')
stim = visual.ImageStim(win, image = filteredNoise, mask=None)
stim.draw()
win.flip()
event.waitKeys()
