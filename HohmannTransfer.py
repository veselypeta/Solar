from Solar import Solar
from Planet import Planet
import numpy as np
import scipy.constants as k
import math

G = k.gravitational_constant

class HohmannTransfer():

    # initialise a simulation
    def __init__(self):
        self.sim = Solar('solardata.csv')

    # Works!
    def semiMajorAxis(self, p1Name, p2Name):
        sun = self.sim.getPlanet('Sun')
        p1 = self.sim.getPlanet(p1Name)
        p2 = self.sim.getPlanet(p2Name)
        r1 = np.linalg.norm(p1.position - sun.position)
        r2 = np.linalg.norm(p2.position - sun.position)
        a = (r1 + r2) / 2
        return a

    # This function calculates the orbital period of the eliptical transfer orbit
    # Calculated Using Keplers 3'rd law -- returns full period of elliptical orbit, 1/2 of this period is then used
    # to calculate phase angle of launch
    def getOrbitalPeriod(self, p1Name, p2Name):
        sun = self.sim.getPlanet('Sun')
        k = (4*math.pow(math.pi, 2))/(sun.mass * G)
        p = math.sqrt(math.pow(self.semiMajorAxis(p1Name, p2Name), 3)* k)
        return p

    # Take two planet names and from which you get a value for velocity to get
    # an eliptical orbit transfer from one to another
    def transferVelocity(self, p1Name, p2Name):
        sun = self.sim.getPlanet('Sun')
        M = sun.mass
        p1 = self.sim.getPlanet(p1Name)
        p2 = self.sim.getPlanet(p2Name)
        v_direction = p1.velocity/np.linalg.norm(p1.velocity)
        r1 = np.linalg.norm(p1.position - sun.position)
        r2 = np.linalg.norm(p2.position - sun.position)
        v_trans = math.sqrt(2*G*M*(r2/(r1*(r1+r2)))) * v_direction
        return v_trans

    # Given two names of planets returns the phase angle at launch to get
    def transferAngle(self, p1Name, p2Name):
        p2 = self.sim.getPlanet(p2Name)
        # we take half of the whole orbital transfer period since we are only concered with it's position at apoapse
        theta = math.pi*(1-2*((self.getOrbitalPeriod(p1Name, p2Name)/2)/self.sim.calculateOrbitalPeriods(p2)))
        return theta

    # give it a phase angle to launch for -- then it get the orbital transfer velocity
    # and return a value for the closest approach to the planet
    def launchSatelite(self, satelliteName, launchPlanetName, targetPlanetName, angle):
        launchPlanet = self.sim.getPlanet(launchPlanetName)
        targetPlanet = self.sim.getPlanet(targetPlanetName)
        # run the simulation until the angle is reached.
        angle = round(angle, 3)

        # When planet is in position of the angle - launch the satellite
        while True:
            self.sim.runTimestep()
            currAngle = round(self.sim.getAngleBetween(launchPlanetName, targetPlanetName), 3)
            # print(currAngle)
            if currAngle == angle:
                launchVelocity = self.transferVelocity(launchPlanetName, targetPlanetName)
                # print(launchVelocity)
                # set the position so that there no effect from the gravity of the launchPlanet.
                satellite = Planet(satelliteName, 100, launchPlanet.position , launchVelocity, launchPlanet.radius, 'b')
                self.sim.planets.append(satellite)
                break

        self.sim.runAnimation()
        # period = self.getOrbitalPeriod(launchPlanetName, targetPlanetName) + self.sim.time
        # minDistance = float('inf')
        # while self.sim.time < (period):
        #     self.sim.runTimestep()
        #     dist = self.sim.calculateDistance(targetPlanet, self.sim.getPlanet(targetPlanetName))
        #     if (dist < minDistance):
        #         minDistance = dist
        #
        # return minDistance






h = HohmannTransfer()

