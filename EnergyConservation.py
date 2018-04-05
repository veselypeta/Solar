from Solar import Solar
import numpy as np
import csv
# here we will check if energy is conserved throughout the simulation of planets
class EnergyConservation(object):

    def __init__(self):
        self.sim = Solar('solardata.csv')
        self.logfile = 'energyConservation.csv'
    # To calculate Gravitational Potential Energy, just use the formula E = mgh
    # where m is mass of planet and g is gravitational acceleration and h is distance from sun
    def calculateGravitationalPotentialEnergy(self, planet):
        # Since energy is a scalar quantity then we just need the magnitude of the vectors involved.
        f = np.linalg.norm(self.sim.calculateForceApplied(planet))
        g = f / planet.mass
        h = np.linalg.norm(planet.position - self.sim.planets[0].position)
        e = planet.mass * g * h
        return e

    def calculateTotalEnergy(self):
        totalEnergy = 0
        for p in self.sim.planets:
            magVel = np.linalg.norm(p.velocity)
            kinetic_energy = 0.5 * p.mass * magVel * magVel
            gravit_energy = self.calculateGravitationalPotentialEnergy(p)
            # I subtract the gravitational energy since it is negative.
            totalEnergy += (kinetic_energy - gravit_energy)
        return totalEnergy

    # write results to log file
    def logEnergy(self):
        with open(self.logfile, 'a') as csvFile:
            writer = csv.writer(csvFile, quoting=csv.QUOTE_MINIMAL, lineterminator='\n')
            writer.writerow([self.calculateTotalEnergy()])

    # Run the simulation
    def runSimulation(self):
        earth = self.sim.getPlanet('Earth')
        year = self.sim.calculateOrbitalPeriods(earth)
        while self.sim.time < year:
            self.sim.runTimestep()
            self.logEnergy()

def __main__():
    a = EnergyConservation()
    a.runSimulation()

__main__()