"""
Abstraction of signal and field objects. 
A signal is a variable depending on time.
A field is a 3d variable depending on time 
defined for a set of points.
Copyright Andres Agudelo-Toro (https://sites.google.com/site/aagudelotoro/)
"""

import numpy as np

class Signal:
    """
    A general class to represent a signal in time.
    """

    def __init__(self,t):
        self.t = t
        self.s = np.zeros_like(self.t)
       
    def load_from_file(self,filename):
        """
        Load time trajectory from a specified
        comma separated file.
        """
        sraw = np.loadtxt(filename,delimiter=',')
        self.s[0:len(sraw)] = sraw[:,1]
    
    def set_raw(self,sraw):
        """
        Make a sequence of values a time signal
        """
        
        #self.s[:len(sraw)] = sraw
        self.s = sraw[:len(self.s)]

    def invert(self):
        self.s = -self.s
    
    def multiply(self,a):
        self.s = a*self.s

    def delay(self,s0,delay):
        """
        Set this signal to be a delayed version of other
        signal. delay is in sample steps.
        """     
        #print self.s.shape,self.s[delay:].shape,(len(self.s)-delay)
        
        #self.s[delay:delay+len(s0.s)] = s0.s
        self.s[delay:] = s0.s[0:len(self.s)-delay]

class StepSignal(Signal):
    """
    Step signal. tstart is how long it takes for the signal 
    to start, that is, be different from 0. 
    The signal ends at tend. 
    t is a regularly spaced time vector. 
    Resulting s is the same size of the time vector.
    """

    def __init__(self, A, t, tstart, tend):
        Signal.__init__(self, t)
        
        #Allocate data
        self.s = np.empty(len(self.t))
        self.set(A, tstart, tend)
        
    """
    Update period, amplitude, start and end
    """
    def set(self, A, tstart, tend):
        dt = self.t[1]
        istart = np.ceil(tstart/dt)
        iend = np.ceil(tend/dt)
        self.s[:] = 0.
        self.s[istart:iend] = A

class BipolarSignal(Signal):
    """
    Square bipolar signal. tstart is how long it takes for the signal 
    to start, that is, be different from 0.
    At tswitch the signal changes sign.
    The signal ends at tend. 
    t is a regularly spaced time vector. 
    Resulting s is the same size of the time vector.
    """

    def __init__(self, A, t, tstart, tswitch, tend):
        Signal.__init__(self, t)
        
        #Allocate data
        self.s = np.empty(len(self.t))
        self.set(A, tstart, tswitch, tend)
        
    """
    Update period, amplitude, start and end
    """
    def set(self, A, tstart, tswitch, tend):
        dt = self.t[1]
        istart = np.ceil(tstart/dt)
        iswitch = np.ceil(tswitch/dt)
        iend = np.ceil(tend/dt)
        self.s[:] = 0.
        self.s[istart:iswitch] = A
        self.s[iswitch:iend] = -A


class CosSignal(Signal):
    """
    Cosinosoindal signal. tstart is not a phase shift but 
    how long it takes for the signal to start, that is, 
    be different from 0. 
    The signal ends at tend. 
    t is a regularly spaced time vector. 
    Resulting s is the same size of the time vector.
    """

    def __init__(self, A, period, t, tstart, tend):
        Signal.__init__(self,t)
        
        #Allocate data
        self.s = np.empty(len(self.t))
        self.set(A,period,tstart,tend)
        
    """
    Update period, amplitude, start and end
    """
    def set(self, A, period, tstart, tend):
        dt = self.t[1]
        istart = np.ceil(tstart/dt)
        iend = np.ceil(tend/dt)
        self.s[:] = 0.
        self.s[istart:iend] = A*np.cos(
            2.*np.pi*(1/period)*self.t[:iend-istart])

class RLCSignal(Signal):
    """
    Electric field induced by the inductor in a Resistor, 
    Inductor, Capacitor circuit. It emulates the electric
    field produced by a transcraneal magnetic stimulation 
    coil.
    tstart is how long it takes for the signal to start.
    The signal ends at tend. 
    t is a regularly spaced time vector. 
    Resulting s is the same size of the time vector.
    """

    def __init__(self, t, tstart, tend):
        Signal.__init__(self,t)
        self.set(tstart,tend)
        
    """
    Update period, amplitude, start and end
    """
    def set(self, tstart, tend):
        dt = self.t[1] - self.t[0]
        istart = np.ceil(tstart/dt)
        iend = np.ceil(tend/dt)
        print istart,iend
        self.s[:] = 0.0
        print self.s[istart:iend].shape
        print self.didt(self.t[:iend-istart]).shape
        self.s[istart:iend] = self.didt(self.t[:iend-istart])
            
    def creepin(self,t):
      """
      A function to represent a more realistic rise trajectory
      of the electric field signal      
      """
      risetime = 5e-6 #Risetime in seconds
      return 1.0 - np.exp(-t/risetime)
     
    def didt(self,t):
      """
      Expresion for the time derivative of the current in the
      inductor
      """
      v0 = 5000.0 #volts
      i0 = 0.0 #amps
      capacitance = 110e-6 #Farad
      resistance = 0.08 #Ohm
      inductance = 12e-6 #Henry
      alpha = resistance/(2*inductance)
      beta = np.sqrt(1.0/(inductance*capacitance) - alpha**2)

      return self.creepin(t)*0.5*(np.sign(t) + 1.0)*0.5*(np.sign(-t + 2.0*np.pi*np.sqrt(inductance*capacitance)) + 1.0)*(-np.exp(-alpha*t)*alpha*np.cos(beta*t)*i0
        -np.exp(-alpha*t)*beta*np.sin(beta*t)*i0
        +np.exp(-alpha*t)*beta*np.cos(beta*t)*
        (-((alpha*i0)/beta) + v0/(inductance*beta)) 
        -np.exp(-alpha*t)*alpha*np.sin(beta*t)*
        (-((alpha*i0)/beta) + v0/(inductance*beta)))

class Field(object):
    """
    A general class to represent a field.
    """
    def __init__(self):
        """
        """    
        #The amplitude factor
        self.A = 1.0

    def amplitude(self,A=None):
        """
        Set a new amplitude factor.
        """
        if A is not None:
            self.A = A
        return self.A

class HomogField(Field):
    """
    Class that produces a spatially homogenous field for a given
    time changing signal. The field is multiplied by the scale
    factor.
    """
    def __init__(self, signal, v=np.array([1.0,0.0,0.0])):
        """
        signal -- the time trajectory
        v -- array with the 3 components of the field (preferibly normalized)
        """
        Field.__init__(self)
        self.s = signal
        self.v = v
        
    def calculate_field(self, X, Y, Z):
        """
        Calculate field for a series of points.
        X,Y,Z must be a 1D arrays forming 3D coordinates.
        It is assumed coordinates are in meters.
        Results are stored in the F variable.
        """
        self.X = X
        self.Y = Y
        self.Z = Z
        
        #Allocate data
        F = np.zeros((len(self.s.t),len(self.X),3))
        
        #Set X component of the field in every point to 
        #have the signal value
        As = self.A*self.s.s
        for i in range(0,len(self.X)):
            F[:,i,0] = As*self.v[0]
            F[:,i,1] = As*self.v[1]
            F[:,i,2] = As*self.v[2]
        
        return F
