from psychopy import visual, core
import numpy as np
import scipy.ndimage

# Create PsychoPy window covering the whole screen
win = visual.Window(size=(1920, 1080), fullscr=True, monitor='testMonitor', units='pix', color=[0, 0, 0])

# Noise properties
noise_duration = 10  # Duration of noise presentation in seconds
frame_rate = 60  # Frames per second
noise_std = 0.1  # Standard deviation for pixel noise
noise_size = win.size  # Use the window size for noise texture

# Blob properties
blob_width_arcmin = 11  # Space constant of the Gaussian blob (standard deviation) in arcmin
blob_amplitude = 1.0  # Amplitude for the Gaussian blob
blob_std = blob_width_arcmin / (2 * np.sqrt(2 * np.log(2)))  # Convert arcmin to std for Gaussian

# Generate and display dynamic Gaussian noise with embedded blob
for frame in range(int(noise_duration * frame_rate)):
    noise = np.random.normal(0, noise_std, size=noise_size)
    
    # Create Gaussian blob
    y, x = np.meshgrid(np.arange(noise_size[1]), np.arange(noise_size[0]))
    blob = blob_amplitude * np.exp(-((x - noise_size[1] / 2)**2 + (y - noise_size[0] / 2)**2) / (2 * blob_std**2))
    
    # Combine noise with the blob
    final_image = blob + noise
    
    # Clip noise values at three standard deviations
    final_image = np.clip(final_image, -3 * noise_std, 3 * noise_std)
    
    # Set the maximum value to match the maximum monitor output
    final_image = final_image / (3 * noise_std)  # Normalize to [-1, 1] range
    final_image = final_image * 255  # Scale to [0, 255]
    
    # Display the combined image
    noise_stim = visual.ImageStim(win, image=final_image, units='pix', size=noise_size)
    noise_stim.draw()
    win.flip()

# Close the window
win.close()
