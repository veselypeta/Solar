# Solar class
import csv
import numpy as np
from Planet import Planet
import scipy.constants as k
import matplotlib.pyplot as plt
import matplotlib.animation as f
import math
# declare G as a constant
G = k.gravitational_constant


class Solar(object):

    def __init__(self, filename):
        self.planets = []
        # read in file to create planets
        self.loadPlanets(filename)

        self.initialisePreviousAccelerations()
        self.timestep = 100
        self.time = 0

    def getAngleBetween(self, planet1, planet2):
        sun = self.getPlanet('Sun')
        p1Vec = planet1.position - sun.position
        p1Mag = np.linalg.norm(p1Vec)
        p2Vec = planet2.position - sun.position
        p2Mag = np.linalg.norm(p2Vec)
        theta = math.acos((p1Vec.dot(p2Vec))/(p1Mag*p2Mag))
        return theta

    def convertSecToDays(self, i):
        days = i/(60 * 60 * 24)
        return days

    def calculateOrbitalPeriods(self, planet):
        r = np.linalg.norm(planet.position)
        v = np.linalg.norm(planet.velocity)
        if v != 0:
            period = (2 * math.pi * r)/v
        else:
            period = 0
        return period

    def relativePeriods(self):
        earth = self.getPlanet('Earth')
        eathPeriod = self.calculateOrbitalPeriods(earth)
        relativePeriods = {}
        for i in self.planets:
            relativePeriods[i.name] = self.calculateOrbitalPeriods(i)/eathPeriod
        return relativePeriods

    def logEnergy(self):
        with open('energyConservation.csv', 'a') as csvFile:
            writer = csv.writer(csvFile, quoting=csv.QUOTE_MINIMAL, lineterminator='\n')
            writer.writerow([self.calculateTotalEnergy()])

    def calculateDistance(self, planet1, planet2):
        s = planet2.position - planet1.position
        dist = np.linalg.norm(s)
        return dist

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
            if p == planet:
                continue
            else:
                s = p.position - planet.position
                forceDirection = self.normaliseVecotr(s)
                magSqr = np.linalg.norm(s) * np.linalg.norm(s)
                if magSqr != 0:
                    f = ((G *p.mass*planet.mass)/(magSqr)) * forceDirection
                else:
                    f = np.zeros(2)

            totalForce = totalForce + f

        return totalForce

    def calculateNextPosition(self, planet):
        # calculate current acceleration applied to the planet from force applied
        # by all other planets from their positions
        planet.currentAcceleration = self.calculateForceApplied(planet)/ planet.mass
        newPosition = planet.position + (planet.velocity*self.timestep) + \
                      (1.0/6.0)*((4*planet.currentAcceleration) - planet.previousAcceleration)\
                      *(self.timestep**2)
        return newPosition

    def calculateNextVelocity(self, planet):
        # this is the new acceleration for next timestep because all positons have been updated
        nextAcceleration = self.calculateForceApplied(planet)/planet.mass

        nextVelocity = planet.velocity + (1.0/6.0)\
                       *((2.0*nextAcceleration) + (5.0*planet.currentAcceleration) - planet.previousAcceleration)*self.timestep
        # change values of previous and current acceleration for the next timestep
        planet.previousAcceleration = planet.currentAcceleration
        planet.currentAcceleration = nextAcceleration
        return nextVelocity

    def runTimestep(self):
        # increment time
        self.time += self.timestep
        # First Calculate all new Positions for all planets, and update their positions
        newPositions = []
        for planet in self.planets:
            newPos = self.calculateNextPosition(planet)
            newPositions.append(newPos)

        # update all positions at once
        for i in range(len(self.planets)):
            self.planets[i].position = newPositions[i]

        # Now we can calculate new velocities based on new positions
        newVelocities = []
        for p in self.planets:
            newV = self.calculateNextVelocity(p)
            p.velocity = newV

    # To calculate Gravitational Potential Energy, just use the formula E = mgh
    # where m is mass of planet and g is gravitational acceleration and h is distance from sun
    def calculateGravitationalPotentialEnergy(self, planet):
        # Since energy is a scalar quantity then we just need the magnitude of the vectors involved.
        f = np.linalg.norm(self.calculateForceApplied(planet))
        g = f/planet.mass
        h = np.linalg.norm(planet.position - self.planets[0].position)
        e = planet.mass * g * h
        return e

    def calculateTotalEnergy(self):
        totalEnergy = 0
        for p in self.planets:
            magVel = np.linalg.norm(p.velocity)
            kinetic_energy = 0.5 * p.mass * magVel * magVel
            gravit_energy = self.calculateGravitationalPotentialEnergy(p)
            totalEnergy += (kinetic_energy + gravit_energy)
        return totalEnergy

    # method used to move the planets in the animation function -- i parameter not used
    def movePlanets(self, i, patches):
        self.runTimestep()
        for i in range(len(patches)):
            patches[i].center = (self.planets[i].position[0], self.planets[i].position[1])

    # This function just returns a max value from all positions
    # used for scaling the image so everything is visible
    def getAnimationRange(self):
        max = 0
        for i in range(len(self.planets)):
            for j in self.planets[i].position:
                if max < j:
                    max = j

        return max

    def initialisePreviousAccelerations(self):
        # set the previous acceleration of the previous timestep for when time = 0 - timestep;
        for p in self.planets:
            initA = self.calculateForceApplied(p)/p.mass
            p.previousAcceleration = initA

    def runAnimation(self):
        fig = plt.figure()
        ax = plt.axes()

        ax.axis('scaled')
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

        animation = f.FuncAnimation(fig, self.movePlanets, fargs=(patches,))
        plt.show()


