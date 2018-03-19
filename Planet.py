# Planet Class
import numpy as np

class Planet(object):

    def __init__(self, name, mass, position, velocity, radius):
        self.name = name
        self.mass = mass
        self.position = position
        self.velocity = velocity
        self.radius = radius

        # initialise the previous acceleration to zero
        self.currentAcceleration = np.array([0,0])
        self.previousAcceleration = np.array([0,0])