from Planet import Planet
import random as rand
import numpy as np
from Solar import Solar


class Asteroid():

    def __init__(self):
        self.sim = Solar('solardata.csv')

    # This function runs the simulation for a given number of years
    def runAsteroidSimulation(self, years):
        earth = self.sim.getPlanet('Earth')
        time = years * self.sim.calculateOrbitalPeriods(earth)
        while self.sim.time < time:
            self.sim.runTimestep()

    # This function checks if the asteroid is near the earth -- i.e. will there be an impact
    def closeEncounter(self, asteroid):
        encounter = False
        earth = self.sim.getPlanet('Earth')
        distance = self.sim.calculateDistance(asteroid, earth)
        # a closest approach is considered around between 10 - 100 thousand km
        if distance < 100*(10**6):
            encounter = True
        return encounter


    # This function will generate a random position for an asteroid that is reasonable.
    # i.e. not too close to the sun and not miles away
    def getRandomPosition(self):
        # generate an integer value that is the maximum value the asteroid can spawn
        max_range = int(self.sim.getAnimationRange())
        rand_x = rand.randrange(-max_range,max_range, 1)
        rand_y = rand.randrange(-max_range,max_range, 1)
        random_position = np.array([rand_x, rand_y])
        return random_position

    # This function will generate a random velocity for an asteroid -- also within a certain range.
    def getRandomVelocity(self):
        max_velocity = 50000
        vel_x = rand.randrange(-max_velocity, max_velocity, 1)
        vel_y = rand.randrange(-max_velocity, max_velocity, 1)
        velocity = np.array([vel_x, vel_y])
        return velocity

    # This function adds an asteroid randomly to the simulation with random velocity and position
    def addAsteroid(self):
        # asteroids typically weigh between 2.8-3.2 *10^21 kg
        asteroid_mass = rand.randrange(280, 320, 1) * 10**19
        asteroid_position = self.getRandomPosition()
        asteroid_velocity = self.getRandomVelocity()
        grey = (0.5, 0.5, 0.5)
        asteroid_radius = 10000
        asteroid = Planet('asteroid', asteroid_mass, asteroid_position, asteroid_velocity, asteroid_radius, grey)
        self.sim.planets.append(asteroid)



i = Asteroid()
i.runAsteroidSimulation(0.01)