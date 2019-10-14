import os
import sys
import ROOT
import matplotlib.pyplot as plt
import datetime

class Sample:
    """
    This is a radioactive sample, that decays over time
    """
    def __init__(self, atoms = 0, decayConstant = 0):
        self._atoms = atoms # number of atoms contained in the sample
        self._decayConstant = decayConstant # decay constant of the sample [s**(-1)]
    
    @property
    def atoms(self):
        # print("Getting atom number")
        return self._atoms
    
    @atoms.setter
    def atoms(self, value):
        # print("Setting number of atoms")
        if value < 0:
            print("ERROR: setting number of atoms below zero!")
            sys.exit(1)
        else:
            self._atoms = value
    
    @property
    def decayConstant(self):
        # print("Getting decayConstant")
        return self._decayConstant
    
    @decayConstant.setter
    def decayConstant(self, value):
        # print("Setting decayConstant")
        if value < 0:
            print("ERROR: setting decayConstant below or equal zero!")
            sys.exit(1)
        else:
            self._decayConstant = value

class Experiment:
    """
    This is an experiment that simulates the radioactive decay of samples
    """
    def __init__(self, samples, timeStep=1, timeEnd=20):
        self._samples = samples # list of samples
        self._timeStep = timeStep # time step for the simulation (s)
        self._timeEnd = timeEnd # time at which the simulation ends (s)
        self._rng = ROOT.TRandom3(datetime.datetime.now().microsecond) # generate a random seed extracted by the system timestamp
        self._simulationData = {} # dictionary containing the simulated data for each dataset
        self._theoryData = {} # dictionary containing the theoretical prediction at each time step, for each dataset
        
        # compute the probabilities that an atom of a sample decays within a timeStep
        stepProbabilities = [sample.decayConstant * timeStep for sample in samples] 
        
        # check that the condition decayConstant * timeStep << 1 (< 0.25) is false and, in that case, set the timeStep to satisfy it
        if max(stepProbabilities ) >= 0.25: 
            maxIndex = stepProbabilities.index(max(stepProbabilities)) # position of the maximum
            self._timeStep = 0.25 / samples[maxIndex].decayConstant
            print('Resetting time step to {timeStep:.2f} s, to satisfy the decayConstant * timeStep << 1 condition'.format(timeStep=self._timeStep))
            
        # time points at which the number of remaining atoms will be simulated
        self._timePoints = [step * self._timeStep for step in range(int(self._timeEnd / self._timeStep) + 1)] 

    def simulate(self,timeEnd=0): # simulate the decay of the samples
        if timeEnd == 0:
            timeEnd = self._timeEnd
        timePoints = [t for t in self._timePoints if t <= timeEnd]
        for sample in self._samples:
            initialAtoms = sample.atoms
            measurements = [] # list containing the simulation results
            
            # list containing the theoretical prediction
            theoryMeasurements = [sample.atoms * ROOT.TMath.Exp(- sample.decayConstant * timePoint) for timePoint in timePoints] 
            
            for timePoint in self._timePoints:
                measurements.append(sample.atoms)
                for atom in range(sample.atoms): 
                    
                    # check if each atom decays and if it does, decrease the size of the sample by 1 (if there are still atoms)
                    if (self._rng.Rndm() < self._timeStep * sample.decayConstant): 
                        if sample.atoms > 0:
                            sample.atoms -= 1
            
            self._simulationData[sample] = measurements
            self._theoryData[sample] = theoryMeasurements
            
            sample.atoms = initialAtoms # restore the number of atoms the sample had at the beginning after the simulation

    # plot the number of atoms vs. time
    def PlotDecayVsTime(self,ax=None,semilogy=False): 
        if not ax:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        
        for sample in self._samples:
            if semilogy:
                plt.semilogy(self._timePoints, self._simulationData[sample], 'o' , label = 'Simulation - Decay constant: '+str(sample.decayConstant)+' 1/s')
                plt.semilogy(self._timePoints, self._theoryData[sample], label = 'Theory - Decay constant: '+str(sample.decayConstant)+' 1/s')
            else:
                plt.plot(self._timePoints, self._simulationData[sample], 'o' , label = 'Simulation - Decay constant: '+str(sample.decayConstant)+' 1/s')
                plt.plot(self._timePoints, self._theoryData[sample], label = 'Theory - Decay constant: '+str(sample.decayConstant)+' 1/s')

        ax.legend()
        ax.set_xlabel('time (s)')
        ax.set_ylabel('Number of atoms')
        ax.set_title('Radioactive Decay Experiment')
        plt.show()
        
    # plot the number of decays per sample within timeInterval for N (= iterations) simulations
    def PlotDecaysInInterval(self,timeInterval,iterations,ax=None): 
        if not ax:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        
        if timeInterval > self._timeEnd:
            print('The time interval must be shorter than the full simulation time!')
            sys.exit(1)
        data = {sample:[] for sample in self._samples}

        for i in range(iterations):
            if i % 50 == 0: # print only once every 50 simulations
                print('Running simulation number: '+str(i)+'/'+str(iterations))
            self.simulate(timeInterval)
            for sample in self._samples:
                # index of the last measurement before timeInterval
                lastMeasurementIndex = [index for index,timePoint in enumerate(self._timePoints) if timePoint < timeInterval][-1] 
                # number of decays happened before timeInterval
                nDecays = self._simulationData[sample][0] - self._simulationData[sample][lastMeasurementIndex] 
                data[sample].append(nDecays)
        
        # last bin of the plot (5 times the most probable number of decays)
        maxDecays = 5 * round(max([timeInterval * sample.atoms * sample.decayConstant for sample in self._samples]))
        
        for sample in self._samples:
            plt.hist(data[sample],label='Decay constant: '+str(sample.decayConstant)+' 1/s',histtype='step', bins=range(maxDecays))
        
        ax.legend()
        ax.set_xlabel('Number of decays')
        ax.set_title('Number of decays in '+str(timeInterval)+' seconds')
        plt.show()

if __name__ == '__main__':
    samples = [
    Sample(500,4e-5), 
    Sample(500,2e-4)
    ]

    exp = Experiment(samples, 1, 2500)
    exp.simulate()
    exp.PlotDecayVsTime(semilogy=False)
    exp.PlotDecaysInInterval(100,1000)