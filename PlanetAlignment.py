from Solar import Solar


class PlanetAlignment(object):

    def __init__(self, tolerance = 10, logfile = 'planetAlignment.txt'):
        self.sim = Solar('solardata.csv')
        self.logfile = logfile
        self.tolerance = tolerance

    # this function cycles through all the plants and calculates an angle from each, to every other (skipping the sun)
    # Then if the angle is less than the tolerance - it continues checking the others - before finally returning True
    # if and only if all the angles are below the tolerance
    def isAligned(self):
        aligned = False
        sun = self.sim.getPlanet('Sun')
        for p1 in self.sim.planets:
            # skip the sun
            if p1 == sun:
                continue
            for p2 in self.sim.planets:
                if p2 == sun:
                    continue
                if p1 != p2:
                    angle = self.sim.getAngleBetween(p2, p1)
                    # if the angle is less than some value, the they are aligned
                    # otherwise not
                    if angle < Solar.degreesToRadians(self.tolerance):
                        aligned = True
                    else:
                        return False
        return aligned

    # Works -- but quite slow
    def runSimulation(self,):
        # run simulation for 1 year before checking if planets are aligned - so they're initially misaligned
        earth = self.sim.getPlanet('Earth')
        year = self.sim.calculateOrbitalPeriods(earth)
        while self.sim.time < year:
            self.sim.runTimestep()

        while not self.isAligned():
            self.sim.runTimestep()
        self.sim.time = 0

        while self.sim.time < year:
            self.sim.runTimestep()

        while not self.isAligned():
            self.sim.runTimestep()

        time_of_alignment = self.sim.time / year
        # the time until the alignment happens is then logged.
        self.logResult(time_of_alignment)

    def logResult(self, time):
        with open(self.logfile, 'a') as myFile:
            message = "The Planets Have Aligned to within " + str(self.tolerance) + " degrees after " + \
                      str(time) + " years." + '\n'
            myFile.write(message)

def __main__():
    for s in range(10, 20):
        a = PlanetAlignment(s)
        a.runSimulation()

__main__()

