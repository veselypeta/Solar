from Planet import Planet
import random as rand
import numpy as np
from Solar import Solar


class Asteroid():

    def __init__(self):
        self.sim = Solar('solardata.csv')
        self.numberOfAsteroids = 0
        self.log = 'asteroidExperimentLog.txt'

    # This function runs the simulation for a given number of years
    def runAsteroidSimulation(self, years, probability):
        self.logString("New Experiment Started for " + str(years) + " Years!")
        earth = self.sim.getPlanet('Earth')
        time = years * self.sim.calculateOrbitalPeriods(earth)
        while self.sim.time < time:
            self.sim.runTimestep()
            if rand.random() < probability:
                self.addAsteroid()
            # check if there is a close encounter and log it
            for obj in self.sim.planets:
                if obj.name == 'asteroid':
                    if self.closeEncounter(obj):
                        self.logCloseEncouter(obj)


    def logCloseEncouter(self, asteroid):
        with open(self.log, 'a') as myFile:
            earth = self.sim.getPlanet('Earth')
            distance = self.sim.calculateDistance(earth, asteroid)
            myFile.write("There was a close encouter with an asteroid at time: " + str(self.sim.time) +
                         " and the asteroid had a mass of " + str(asteroid.mass)+ ". The asteroid was " +
                         str(distance) + "m  from earth" + '\n')

    def logAsteroid(self, asteroid):
        with open(self.log, 'a') as myFile:
            myFile.write("An asteroid was randomly created with a position of " + str(asteroid.position) +
                         " and a velocity of " + str(asteroid.velocity) + ". There are now " +
                         str(self.numberOfAsteroids) +  '\n')

    def logString(self, string):
        with open(self.log, 'a') as myFile:
            myFile.write(string + '\n')


    # This function checks if the asteroid is near the earth -- i.e. will there be an impact
    def closeEncounter(self, asteroid):
        encounter = False
        earth = self.sim.getPlanet('Earth')
        distance = self.sim.calculateDistance(asteroid, earth)
        # a closest approach is considered around between 10 - 100 thousand km
        # had to increase the distance so that we can log how far away it is.
        if distance < 100*(10**8):
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
        self.numberOfAsteroids += 1
        # asteroids typically weigh between 2.8-3.2 *10^21 kg
        asteroid_mass = rand.randrange(280, 320, 1) * 10**19
        asteroid_position = self.getRandomPosition()
        asteroid_velocity = self.getRandomVelocity()
        grey = (0.5, 0.5, 0.5)
        asteroid_radius = 10000
        asteroid = Planet('asteroid', asteroid_mass, asteroid_position, asteroid_velocity, asteroid_radius, grey)
        self.logAsteroid(asteroid)
        self.sim.planets.append(asteroid)

def __main__()
    i = Asteroid()
    i.runAsteroidSimulation(4, 0.1)

__main__()