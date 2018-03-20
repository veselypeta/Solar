from Solar import Solar
from Planet import Planet

class HohmannTransfer():

    def __init__(self):
        # Initialise a simulation of the solar system
        sim = Solar('solardata.csv')
        earth = sim.getPlanet('Earth')
        # run the simulation for 2 years
        two_years = sim.calculateOrbitalPeriods(earth) * 2
        while sim.time < two_years:
            sim.runTimestep()
        mySatellite = Planet('MarsProbe', 10, earth.position, earth.velocity * 1.1, earth.radius, 'b')
        sim.planets.append(mySatellite)
        sim.runAnimation()

a = HohmannTransfer()

