from Solar import Solar
from Planet import Planet
import numpy as np
import scipy.constants as k
import math


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
        k = (4*math.pow(math.pi, 2))/(sun.mass * Solar.G)
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

    def launchProbleToMars(self, angle):
        earth = self.sim.getPlanet('Earth')
        mars = self.sim.getPlanet('Mars')
        closest_approach = float('inf')
        # theoretical_angle = self.transferAngle('Earth', 'Mars')
        # range = Solar.degreesToRadians(10)
        # test_angles = np.arange(theoretical_angle-range, theoretical_angle+range, Solar.degreesToRadians(1))

        # wait 2 years initially to have the planets not all in alignment
        two_years = self.sim.calculateOrbitalPeriods(earth) * 0.5
        while self.sim.time < two_years:
            self.sim.runTimestep()

        while True:
            self.sim.runTimestep()
            angle = round(angle, 3)
            earth_mars_angle = round(self.sim.getAngleBetween(earth, mars),3)
            if angle == earth_mars_angle:
                # launch Sattelite


                # run the simulation for 1 period of the eliptical orbit of the probe
                simulation_time = self.sim.time + self.getOrbitalPeriod('Earth', 'Mars')
                probe_velocity = self.transferVelocity('Earth', 'Mars')
                mars_probe = Planet('Mars Probe', 100, earth.position, probe_velocity, earth.radius, (0.6, 0.6, 0.6))
                self.sim.planets.append(mars_probe)
                self.sim.runAnimation()
                while self.sim.time < simulation_time:
                    self.sim.runTimestep()
                    dist = self.sim.calculateDistance(self.sim.getPlanet('Mars Probe'), self.sim.getPlanet('Mars'))
                    if dist < closest_approach:
                        closest_approach = dist
                return closest_approach




h = HohmannTransfer()
r = 44/180.0 * math.pi
print(h.getOrbitalPeriod('Earth', 'Mars'))
d = h.launchProbleToMars(r)
print(d)
