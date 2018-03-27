from Solar import Solar


class PlanetAlignment(object):

    def __init__(self):
        self.sim = Solar('solardata.csv')

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
                    if angle < Solar.degreesToRadians(5):
                        aligned = True
                    else:
                        return False
        return aligned

    # Works -- but quite slow
    def runSimulation(self):
        # run simulation for 1 year before checking if planets are aligned
        earth = self.sim.getPlanet('Earth')
        year = self.sim.calculateOrbitalPeriods(earth)
        while self.sim.time < year:
            self.sim.runTimestep()

        while not self.isAligned():
            self.sim.runTimestep()

        print(Solar.degreesToRadians(5))
        merc = self.sim.getPlanet('Mercury')
        print(self.sim.getAngleBetween(earth, merc))
        time_of_alignment = self.sim.time / year
        return time_of_alignment





