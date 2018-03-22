from Solar import Solar
from Planet import Planet
import numpy as np
import matplotlib.pyplot as plt

class HohmannTransfer():

    def __init__(self, multiple):
        # Initialise a simulation of the solar system
        self.multiple = multiple
        self.closestApproach = 0
        sim = Solar('solardata.csv')
        earth = sim.getPlanet('Earth')
        # run the simulation for 2 years
        two_years = sim.calculateOrbitalPeriods(earth) * 2
        while sim.time < two_years:
            sim.runTimestep()

        # create a satellite with the same position as earth, but a velocity which is a multiple of earths.
        mySatellite = Planet('MarsProbe', 10, earth.position, earth.velocity * multiple, earth.radius, 'b')
        sim.planets.append(mySatellite)


        sim.time = 0
        minDistance = float('inf')
        # run for another two years
        while sim.time < two_years:
            dist = sim.calculateDistance(sim.getPlanet('Mars'), sim.getPlanet('MarsProbe'))
            if ( dist < minDistance):
                minDistance = dist
            sim.runTimestep()
            # self.angles.append(sim.getAngleBetween(sim.getPlanet('Earth'), sim.getPlanet('Mars')))
        self.closestApproach = minDistance


# a = HohmannTransfer(1.1)
# print(a.closestApproach)
b = {}
# for i in np.arange(1, 1.2, 0.01):
#     a = HohmannTransfer(i)
#     b[i] = a.closestApproach




