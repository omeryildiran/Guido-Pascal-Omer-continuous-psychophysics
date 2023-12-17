from psychopy import sound, core

# Initialize the sound objects
right_sound = sound.Sound('right.wav')
left_sound = sound.Sound('left.wav')

# Initialize a clock
clock = core.Clock()

# Create a loop to continuously play the sounds
while True:
    right_sound.play()  # Play right.wav
    while clock.getTime() < right_sound.getDuration():
        pass
    clock.reset()  # Reset the clock
    right_sound.stop()  # Stop the sound
    left_sound.play()  # Play left.wav
    while clock.getTime() < left_sound.getDuration():
        pass
    left_sound.stop()  # Stop the sound
    clock.reset()  # Reset the clock
