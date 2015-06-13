"""
Implements the NEURON calls of different
methods for extra and some intracellular
stimulation.

The mechanism objects take a list of neuronal 
sections and assigns the right NEURON mechanisms
for stimulation.
Copyright Andres Agudelo-Toro (https://sites.google.com/site/aagudelotoro/)
"""

class StimulationMechanism(object):
    """
    """
    def __init__(self):
        """
        """

from numpy import array,dot,sqrt,c_
from pylab import plot,figure,title

class ActivatingFunction(StimulationMechanism):
    """
    For stimulation, the amplitude of the electric field 
    needs to be converted to the equivalent current injection. 
    This depends on the axial resistance ra.
    The current injection at the 
    ends is proportional to the electric field and not 
    the potential. This is due to the sealed ends condition.
    A current clamp is used to simulate this.
    """
    
    def __init__(self,h,sections):
        StimulationMechanism.__init__(self)
        
        self.h = h
        self.nrnclamp = []
        self.nrnvec = []
        self.sections = sections
        #Allocate NEURON Vectors for each section to store current
        #values and create current clamps for initial and final ends
        for s in self.sections.list:
            self.nrnvec.append((self.h.Vector(),self.h.Vector()))

            stim0 = self.h.IClamp(s.nrnsec(0.0))
            stim1 = self.h.IClamp(s.nrnsec(1.0))
            
            #For some weird reason NEURON requires
            #the clamps to have their variables set
            #even tough these values will not be used
            stim0.amp = 0
            stim0.dur = 1e9
            stim1.amp = 0
            stim1.dur = 1e9

            #The reference to the stimulators needs to be
            #kept to avoid garbage collection
            self.nrnclamp.append([stim0, stim1])
            
            #DEBUG
            #print 'Clamps created'

    def setstimuli(self,E,dt):
        """
        Set electric field.
        E should be in units of V/m.
        There must be a field value for each section. 
        This assumes the sections are small enough that
        only one field value corresponds to them.
        The field is a 3D array where each row corresponds 
        to time, each column to an independent segment, and 
        the 3rd dimension to the Ex,Ey,Ez components.
        dt is in milliseconds.
        """
        
        #Dummy electric field, for debuging
        #E = array([
        #    [[0,0,0],[0,0,0]],
        #    [[280,0,0],[280,0,0]],
        #    [[280,0,0],[280,0,0]],
        #    [[280,0,0],[280,0,0]],
        #    [[0,0,0],[0,0,0]]])
        
        #Get number of timesteps
        nt = len(E)
        
        #Prepare NEURON Vectors of each section to store field
        #values
        for j,s in enumerate(self.sections.list):
            self.nrnvec[j][0].resize(nt)
            self.nrnvec[j][1].resize(nt)
             
        #For each time step find the field projection
        #on each neuronal section
        for i in range(nt):
            for j,s in enumerate(self.sections.list):

                #Find projection of the field over the sections's unitary
                #vector u
                p0,p1 = s.points()
                u = p1 - p0
                u = u/sqrt(dot(u,u))
                Eu = dot(E[i,j,:],u)
                
                #DEBUG
                if i==int(nt*0.25):
                    print '%s %.3g %.3g %.3g %.3g %.3g'%(s.name(), Eu, s.ra(), s.nrnsec.Ra, s.d(), 1e9*Eu/s.ra())
                
                #Store current values in NEURON Vector objects
                #for future simulation. Current is converted 
                #to nano amps first
                ra = s.ra()  #Here ra might not be correct.
                             #It could include the extracellular resistance!
                             #See Monai et al. 2010 Biophys. J.
                self.nrnvec[j][0].x[i] = -1e9*Eu/ra
                self.nrnvec[j][1].x[i] = +1e9*Eu/ra
                
        #Set current clamps for both section ends to represent 
        #the effect of the activating function.
        for j,s in enumerate(self.sections.list):
            #This is the way current clamps are given a time trajectory 
            #in NEURON (See Vector.play())
            
            self.nrnvec[j][0].play(self.nrnclamp[j][0]._ref_amp,dt)
            self.nrnvec[j][1].play(self.nrnclamp[j][1]._ref_amp,dt)
                        
            #DEBUG
            #print 'Play added'
        
        #DEBUG: see neuron first vector
        #f=figure(id(self))
        #plot(c_[self.nrnvec[0][0],self.nrnvec[0][1]],hold=False)
        #title('Ix0(t),Iy0(t)')
        
class ExtracellularMechanism(StimulationMechanism):
    """
    Uses NEURON's extracellular mechanism to set the extra-
    cellular potential for a set of sections.
    Each NEURON section is divided in segments.
    For every segment on every section, a potential value over
    time must be provided. This is required by NEURON because 
    the Vector.play method is used to drive the e_extracellular
    variable of each segment. For further reference see the 
    NEURONS mechanisms documentation and these forum posts:
    
    https://www.neuron.yale.edu/phpBB2/viewtopic.php?t=168
    http://www.neuron.yale.edu/phpBB/viewtopic.php?t=212
    http://www.neuron.yale.edu/phpBB/viewtopic.php?t=1814
    """
    
    def __init__(self,h,sections):
        StimulationMechanism.__init__(self)
        
        self.h = h
        self.secvec = []
        self.sections = sections
        #Allocate NEURON Vectors for each segment to store field
        #values. Also prepare for stimulation.
        for s in self.sections.list:
            segvec = []
            
            for i in range(0,s.nseg()):
                segvec.append(self.h.Vector())

            self.secvec.append(segvec)
            
            #Insert extracellular mechanism
            s.nrnsec.insert('extracellular')


    def setstimuli(self,phi,dt):
        """
        Set potential field.
        There must be a field value for each section. 
        The field is a list of 2D arrays, each corresponding
        to a section. Each row of the 2D array corresponds to 
        time and each column to each segment. 
        dt is in milliseconds.
        """
   
        #Get number of timesteps
        nt = len(phi[0][:,0])
        
        #For the HOC-Vector.play trick
        self.h('objref vect')
        
        #Fill NEURON Vectors of each section to store field
        #values
        for i,s in enumerate(self.sections.list):
            ns = s.nseg()
            #Spatial differential for the section
            ds = 1./ns
            
            #For each segment
            for j in range(0,ns):
                self.secvec[i][j].resize(nt)
                
                #For each time step
                for k in range(0,nt):
                    self.secvec[i][j].x[k] = phi[i][k,j]
                    
                #Prepare for vector assignment during simulation
                #self.secvec[i][j].play(s.e_extracellular[ds/2 + j*ds]._ref_amp,
                #    dt*1e3)
                
                #I couldn't find a way to get the extracellular mechanism 
                #varible references, so the play method can't be set from
                #python. Next code uses a trick to set the play from HOC.
                
                #Set as the currently accessed section
                s.nrnsec.push()
                                
                self.h.vect = self.secvec[i][j]
                
                self.h('vect.play(&e_extracellular('+
                    str(ds/2 + j*ds)+'),'+str(dt)+')')
                
                #Remove this section from the HOC pile
                self.h.pop_section()
            
class ElectricClamp(StimulationMechanism):
    """
    Amplitude is given in terms of electric field not current
    making it dependent on axial resistance.
    For stimulation, the amplitude of the electric field 
    needs to be converted to the equivalent current injection. 
    This depends on the axial resistance ra.
    A current clamp is used.
    Only the X component of E is considered.
    This is useful for single or clamp tests.
    """
    
    def __init__(self,h,sections,locations):
        StimulationMechanism.__init__(self)
        
        self.h = h
        self.nrnclamp = []
        self.nrnvec = []
        self.sections = sections
        self.locations = locations
        
        #Allocate NEURON Vectors for each section to store field
        #values and create current clamps for each one
        for s,l in zip(self.sections.list,self.locations):
            self.nrnvec.append(self.h.Vector())

            stim = self.h.IClamp(s.nrnsec(l))
            #For some weird reason NEURON requires
            #the clamps to have their variables set
            #even tough these values will not be used
            stim.amp = 0
            stim.dur = 1e9

            #The reference to the stimulators needs to be
            #kept to avoid garbage collection
            self.nrnclamp.append(stim)

    def setstimuli(self,E,dt):
        """
        Set electric field.
        There must be a field value for each section. 
        The field is a 3D array where each row corresponds 
        to time, each column to an independent segment, and 
        the 3rd dimension to the Ex,Ey,Ez components.
        dt is in milliseconds.
        """
        
        #Get number of timesteps
        nt = len(E)
        
        #Prepare NEURON Vectors of each section to store field
        #values
        for j,s in enumerate(self.sections.list):
            self.nrnvec[j].resize(nt)
       
        #For each time step find set the current
        #on each neuronal section
        for i in range(nt):
            for j,s in enumerate(self.sections.list):

                #Store current values in NEURON Vector objects
                #for future simulation. Current is converted 
                #to nano amps
                self.nrnvec[j].x[i] = 1e9*E[i,j,0]/s.ra()

        #Set current clamps
        for j,s in enumerate(self.sections.list):
            #This is the way current clamps are given a time trajectory 
            #in NEURON (See Vector.play())
            self.nrnvec[j].play(self.nrnclamp[j]._ref_amp,dt)

class CurrentClamp(StimulationMechanism):
    """
    A Python wrapper to NEURON's IClamp mechanism. This
    implementation only inserts a current clamp on every 
    section at the given range value in list p in [0,1]. 
    An arbitrary current trajectory can be given.
    """
    
    def __init__(self,h,sections,p):
        StimulationMechanism.__init__(self)
        
        self.h = h
        self.nrnclamp = []
        self.nrnvec = []
        self.sections = sections
        #Allocate NEURON Vectors for each section to store current
        #values and create a current clamp in position p
        for i,s in enumerate(self.sections.list):
            self.nrnvec.append(self.h.Vector())

            stim = self.h.IClamp(s.nrnsec(p[i]))
            
            #For some weird reason NEURON requires
            #the clamps to have their variables set
            #even tough these values will not be used
            stim.amp = 0
            stim.dur = 1e9

            #The reference to the stimulators needs to be
            #kept to avoid garbage collection
            self.nrnclamp.append(stim)

    def setstimuli(self,I,dt):
        """
        Specify current trajectory in time for each clamp. I must 
        be a 2D array were the row index corresponds to time.
        Current is in nanoamperes, dt is in milliseconds.
        """
        
        #Get number of timesteps
        nt = len(I[:,0])
        
        #Prepare the NEURON Vectors of each section that
        #will store current values
        for j,s in enumerate(self.sections.list):
            self.nrnvec[j].resize(nt)
       
        #For each time step assign current value on each section
        for i in range(nt):
            for j,s in enumerate(self.sections.list):

                #Store current values in NEURON Vector objects
                #for future simulation
                self.nrnvec[j].x[i] = I[i,j]

        #Set current clamps for each section.
        for j,s in enumerate(self.sections.list):
            #This is the way current clamps are given a time trajectory 
            #in NEURON (See Vector.play())
            self.nrnvec[j].play(self.nrnclamp[j]._ref_amp,dt)

