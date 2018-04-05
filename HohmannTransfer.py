from Solar import Solar
from Planet import Planet
import numpy as np
import math
import csv


class HohmannTransfer():

    # initialise a simulation
    def __init__(self):
        self.sim = Solar('solardata.csv', 10000)
        self.logfile = "marsExperiment.txt"
        self.rawData = "marsExperimentData.csv"

    # this function returns the clockwise angle of a planet and the x-axis
    def getClockwiselockwiseAngle(self, planet):
        sun = self.sim.getPlanet('Sun')
        vec = planet.position - sun.position
        a = math.atan2(vec[1], vec[0])
        if a < 0:
            a = 2*math.pi + a
        return a

    # return the semi-major axis of the ellipse of the Hohmann transfer between two planets.
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
    # an elliptical orbit transfer from one to another
    def transferVelocity(self, p1Name, p2Name):
        sun = self.sim.getPlanet('Sun')
        M = sun.mass
        p1 = self.sim.getPlanet(p1Name)
        p2 = self.sim.getPlanet(p2Name)
        v_direction = p1.velocity/np.linalg.norm(p1.velocity)
        r1 = np.linalg.norm(p1.position - sun.position)
        r2 = np.linalg.norm(p2.position - sun.position)
        v_trans = math.sqrt(2*Solar.G*M*(r2/(r1*(r1+r2)))) * v_direction
        return v_trans

    # Given two names of planets returns the phase angle at launch to get - this is the ideal angle - theoretically
    def transferAngle(self, p1Name, p2Name):
        p2 = self.sim.getPlanet(p2Name)
        # we take half of the whole orbital transfer period since we are only concered with it's position at apoapse
        theta = math.pi*(1-2*((self.getOrbitalPeriod(p1Name, p2Name)/2)/self.sim.calculateOrbitalPeriods(p2)))
        return theta

    def launchProbe(self, launchAngle):
        # shift the angle by 180 degrees since i'm launching the satellite from the opposite side
        launchAngle = round(launchAngle + math.pi, 3)
        earth = self.sim.getPlanet('Earth')
        mars = self.sim.getPlanet('Mars')
        # run the simulation for two years
        two_years = 2 * self.sim.calculateOrbitalPeriods(earth)
        while self.sim.time < two_years:
            self.sim.runTimestep()

        # I'm rounding the angles to 3 d.p. because otherwise they would be too precise to equate
        marsAngle = round(self.getClockwiselockwiseAngle(mars) - self.getClockwiselockwiseAngle(earth), 3)
        while marsAngle != launchAngle:
            self.sim.runTimestep()
            marsAngle = round(self.getClockwiselockwiseAngle(mars) - self.getClockwiselockwiseAngle(earth), 3)

        # I've placed the probe on the opposite side of the earth, just to ensure that the gravity from the
        # earth does not affect the trajectory
        probe_velocity = -self.transferVelocity('Earth','Mars')
        probe_position = -earth.position
        probe = Planet('Mars Probe', 100, probe_position, probe_velocity, 10000, 'b')
        self.sim.planets.append(probe)
        closestApproach = float('inf')
        # i'm running the simulation for 1 complete orbit of the probe
        simulation_time = self.sim.time + self.getOrbitalPeriod('Earth', 'Mars')
        while self.sim.time < simulation_time:
            self.sim.runTimestep()
            a = self.sim.getPlanet('Mars Probe')
            b = self.sim.getPlanet('Mars')
            distance = self.sim.calculateDistance(a, b)
            if distance < closestApproach:
                closestApproach = distance
        self.logClosestApproach(closestApproach, Solar.radiansToDegrees(marsAngle - math.pi))
        self.logData(closestApproach, Solar.radiansToDegrees(marsAngle - math.pi))

    # write closest approach to file
    def logClosestApproach(self, distance, angle):
        with open(self.logfile, 'a') as myFile:
            message = "The probe was launched to mars at an angle of " + str(angle) \
                      + " and was at a closest approach of " + str(distance) + " m." + '\n'
            myFile.write(message)

    # Each row in the csv file has two fields - angle, distance
    def logData(self, distance, angle):
        with open(self.rawData, 'a') as myFile:
            writer = csv.writer(myFile, quoting=csv.QUOTE_MINIMAL, lineterminator='\n')
            writer.writerow([angle, distance])

# display the orbital period of the Hohmann transfer
def orbitalPeriodOfTransfer():
    h = HohmannTransfer()
    t = h.getOrbitalPeriod('Earth', 'Mars')
    t = h.sim.convertSecToDays(t)
    print("Orbital Period of probe" + str(t))

# Print the theoretical angle of the hohmann transfer.
def idealAngle():
    h = HohmannTransfer()
    a = h.transferAngle('Earth', 'Mars')
    print("Ideas launch angle: " + str(Solar.radiansToDegrees(a)))



def main():
    angles = np.arange(30, 50, 1)
    for i in range(len(angles)):
        h = HohmannTransfer()
        h.launchProbe(Solar.degreesToRadians(angles[i]))

orbitalPeriodOfTransfer()
idealAngle()
main()


