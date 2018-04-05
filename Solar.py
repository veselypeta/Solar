# Solar class
import csv
import numpy as np
from Planet import Planet
import scipy.constants as k
import matplotlib.pyplot as plt
import matplotlib.animation as f
import math



class Solar(object):
    # declare G as a constant -- static/class variable
    G = k.gravitational_constant

    # a class method to convert degrees to radians
    @classmethod
    def degreesToRadians(cls, degrees):
        r = (degrees/180.0) * math.pi
        return r

    # a class method to convert radians to degrees
    @classmethod
    def radiansToDegrees(cls, rad):
        d = rad/math.pi * 180
        return d

    # method takes a filename for planet data and a value for timestep (kwarg) - default value of 100,000
    def __init__(self, filename, timestep=100000):
        self.planets = []
        # read in file to create planets
        self.loadPlanets(filename)
        # Initialize previous acceleration initialy for each planet to be the same value as the first tiemstep
        self.initialisePreviousAccelerations()
        self.timestep = timestep
        self.time = 0

    # this method will accept if a planet object or a planet name is passed in as a parameter
    # This method will be used to get a value for when a satellite can be launched to mars
    def getAngleBetween(self, planet1, planet2):
        # check if object or planet name is passed in
        sun = self.getPlanet('Sun')
        if type(planet1) == str:
            planet1 = self.getPlanet(planet1)
            planet2 = self.getPlanet(planet1)
        # get a vector from the sun to the each planet
        p1Vec = planet1.position - sun.position
        p1Mag = np.linalg.norm(p1Vec)
        p2Vec = planet2.position - sun.position
        p2Mag = np.linalg.norm(p2Vec)
        # using the dot product formula to calculate angle between the two planets
        theta = math.acos((np.dot(p1Vec, p2Vec))/(p1Mag*p2Mag))
        return theta

    # function which converts seconds to days on earth
    def convertSecToDays(self, i):
        days = i/(60 * 60 * 24)
        return days

    # Calculate the period of orbit of a planet using keplers laws
    def calculateOrbitalPeriods(self, planet):
        # magnitude of the planets position and velocity
        r = np.linalg.norm(planet.position)
        v = np.linalg.norm(planet.velocity)
        # this is a check to make sure the mag of velocity isn't zero - i.e. the sun
        if v != 0:
            period = (2 * math.pi * r)/v
        else:
            period = 0
        return period

    # given two planets, returns the distance between them
    def calculateDistance(self, planet1, planet2):
        s = planet2.position - planet1.position
        dist = np.linalg.norm(s)
        return dist

    # given a name, this method returns a Planet object with that name
    def getPlanet(self, name):
        for p in self.planets:
            if p.name == name:
                return p
        return "No such planet exists with that name"

    # This reads in planets from a csv file and adds the planets to the self.planets list
    def loadPlanets(self, filename):
        with open(filename) as planetFile:
            csvReader = csv.reader(planetFile)
            for planet in csvReader:
                name = planet[0]
                mass = float(planet[1])
                position = np.array([float(planet[2]), float(planet[3])])
                velocity = np.array([float(planet[4]), float(planet[5])])
                radius = float(planet[6])
                colour = str(planet[7])
                p = Planet(name, mass, position, velocity, radius, colour)
                self.planets.append(p)

    # A function used to return a normalised vector of a vecotr passsed in
    # i.e. return a unit vector in the same direction as the vector passed in
    def normaliseVecotr(self, v):
        mag = np.linalg.norm(v)
        if mag == 0:
            return v
        else:
            return v/mag

    # Calculates force applied to a planet by summing all individual forces applied to it
    # from all other planets
    def calculateForceApplied(self, planet):
        # create a zero force initially
        totalForce = np.array([0,0])
        # go through every planet to calculate force applied from that planet
        for p in self.planets:
            # skip for itself
            if p == planet:
                continue
            else:
                s = p.position - planet.position
                forceDirection = self.normaliseVecotr(s)
                magSqr = np.linalg.norm(s) * np.linalg.norm(s)
                # avoid division by zero i.e. for the sun
                if magSqr != 0:
                    f = ((Solar.G *p.mass*planet.mass)/(magSqr)) * forceDirection
                else:
                    f = np.zeros(2)

            totalForce = totalForce + f

        return totalForce

    # this method calculates the next position of a planet at the next timestep - using the beeman algorithm
    def calculateNextPosition(self, planet):
        # calculate current acceleration applied to the planet from force applied
        # by all other planets from their positions
        planet.currentAcceleration = self.calculateForceApplied(planet)/ planet.mass
        # beeman algorithm
        newPosition = planet.position + (planet.velocity*self.timestep) + \
                      (1.0/6.0)*((4*planet.currentAcceleration) - planet.previousAcceleration)\
                      *(self.timestep**2)
        return newPosition

    # calculate next velocity at a timestep for a given planet
    def calculateNextVelocity(self, planet):
        # this is the new acceleration for next timestep because all positions have been updated
        nextAcceleration = self.calculateForceApplied(planet)/planet.mass
        # beeman algorithm
        nextVelocity = planet.velocity + (1.0/6.0)\
                       *((2.0*nextAcceleration) + (5.0*planet.currentAcceleration) - planet.previousAcceleration)*self.timestep
        # change values of previous and current acceleration for the next timestep
        planet.previousAcceleration = planet.currentAcceleration
        planet.currentAcceleration = nextAcceleration
        return nextVelocity

    # completes and updates all velocities and position for a complete timestep of the simulation
    def runTimestep(self):
        # increment time)
        self.time += self.timestep
        # log total energy to file
        # First Calculate all new Positions for all planets, and update their positions
        newPositions = []
        for planet in self.planets:
            newPos = self.calculateNextPosition(planet)
            newPositions.append(newPos)

        # update all positions at once
        for i in range(len(self.planets)):
            self.planets[i].position = newPositions[i]

        # Now we can calculate new velocities based on new positions - all objects have new values for positions
        newVelocities = []
        for p in self.planets:
            newV = self.calculateNextVelocity(p)
            p.velocity = newV

    # method used to move the planets in the animation function -- i parameter not used
    def movePlanets(self, i, patches):
        # run a timestep
        self.runTimestep()
        for i in range(len(patches)):
            patches[i].center = (self.planets[i].position[0], self.planets[i].position[1])

    # This function just returns a max magnitude of all position vectors
    # used for scaling the image so everything is visible in the plot
    def getAnimationRange(self):
        max = 0
        for i in range(len(self.planets)):
            posMag = np.linalg.norm(self.planets[i].position)
            if max < posMag:
                max = posMag

        return max

    # i use this method to initialise accelerations for when t=0
    def initialisePreviousAccelerations(self):
        # set the previous acceleration of the previous timestep for when time = 0 - timestep;
        for p in self.planets:
            initA = self.calculateForceApplied(p)/p.mass
            p.previousAcceleration = initA

    def runAnimation(self):
        # create a figure and axes
        fig = plt.figure()
        ax = plt.axes()

        ax.axis('scaled')
        # set the range of the plot
        range = self.getAnimationRange() * 1.5
        ax.set_xlim(-range, range)
        ax.set_ylim(-range, range)

        # create a patches list to hold all patches
        patches = []
        # add patches to the list as planets, with their positions as their centre
        radiusScaleFactor = 30000
        for p in self.planets:
            c = plt.Circle((p.position[0], p.position[1]), p.radius * radiusScaleFactor, color=p.colour)
            ax.add_patch(c)
            patches.append(c)

        # run the animation and show the plot
        animation = f.FuncAnimation(fig, self.movePlanets, fargs=(patches,))
        plt.show()

# This is the code to run the animation -- i've commented it out so it doesn't take too long to run all experiments
# def __main__():
#     s = Solar('solardata.csv')
#     s.runAnimation()
#
# __main__()