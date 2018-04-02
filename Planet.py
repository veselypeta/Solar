# Planet Class
import numpy as np

class Planet(object):

    def __init__(self, name, mass, position, velocity, radius, colour):
        self.name = name            # a name for the planet
        self.mass = mass            #
        self.position = position    # numpy 2D array to represent a 2D vector with [x,y] coordinates
        self.velocity = velocity    # numpy 2D array to represent a 2D velocity vecotor with [x,y] velocities
        self.radius = radius        # Integer value used to scale the planets on screen
        self.colour = colour        # Is either a character or a triple to represent a colour in the animation

        # initialise the previous and current acceleration to zero
        self.currentAcceleration = np.array([0,0])
        self.previousAcceleration = np.array([0,0])