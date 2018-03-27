from Solar import Solar

# this class is to create a log file to outplut the relative periods of orbits
# as fractions of an earth year
class OrbitalPeriods(object):

    #initialise a simulation and the filename of the logfile
    def __init__(self):
        self.sim = Solar('solardata.csv')
        self.logfile = 'relativeOrbitalPeriods.txt'

    # use the method to calculate the periods and then divide by the period of earth
    # return a dictionary file with keys being the planet names and values as the relative periods
    def relativePeriods(self):
        earth = self.sim.getPlanet('Earth')
        eathPeriod = self.sim.calculateOrbitalPeriods(earth)
        relativePeriods = {}
        for i in self.sim.planets:
            relativePeriods[i.name] = self.sim.calculateOrbitalPeriods(i)/eathPeriod
        return relativePeriods

    # write the data to a file
    def logRelativePeriods(self):
        rel_periods = self.relativePeriods()
        with open(self.logfile, 'w') as myFile:
            for val in rel_periods:
                myFile.write(str(val) + " has period " + str(rel_periods[val]) + " Years" + '\n')

# run code
a = OrbitalPeriods()
a.logRelativePeriods()