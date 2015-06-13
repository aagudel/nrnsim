#Section.py: Two classes that help manipulating NEURON geometries
#based on the NEURON Python interface. These clases are simply 
#wrapers around NEURON Python functions.
#Copyright Andres Agudelo-Toro (https://sites.google.com/site/aagudelotoro/)


import numpy as np
from numpy import array,pi,cos,sin,dot,empty

class Section():
    """
    A Class to wrap NEURON sections providing a cleaner
    interface to obtain information about them
    """
    
    def __init__(self,h,nrnsec):
        self.nrnsec = nrnsec

        #Spatial differential for the non dimensional section
        self.ds = 1./self.nrnsec.nseg

        #Record voltage at extremes
        self.nrnV0 = h.Vector()
        self.nrnV0.record(nrnsec(0.0)._ref_v)
        #self.nrnV0.record(nrnsec(self.ds/2)._ref_v)
        self.nrnV1 = h.Vector()
        self.nrnV1.record(nrnsec(1.0)._ref_v)
        #self.nrnV1.record(nrnsec(1.-self.ds/2)._ref_v)
        self.h = h
        
    def recordvoltage(self):
        """
        Enable membrane voltage recording. Create a neuron
        vector for each segment. It uses NEURON vectors.
        Values can be accessed in the nrnV attribute.
        """
        self.nrnV = []

        #For each segment
        for i in range(0,self.nrnsec.nseg):
            v = self.h.Vector()
            v.record(self.nrnsec(self.ds/2. + i*self.ds)._ref_v)
            #DEBUG
            print i,self.ds/2. + i*self.ds
            self.nrnV.append(v)
    
    def recordcurrent(self):
        """
        Enable membrane's current recording. Create a neuron
        vector for each segment. This requires the extracellular
        mechanism.
        """
        self.nrnI = []

        #Insert extracellular mechanism just in case
        self.nrnsec.insert('extracellular')

        #For the HOC Vector.record trick
        self.h('objref rvec')

        #Set the currently accessed section
        self.nrnsec.push()

        #For each segment
        for i in range(0,self.nrnsec.nseg):
            self.h.rvec = self.h.Vector() 
            
            #I couldn't find a way to get the extracellular mechanism 
            #varible references so the play method can't be set from
            #python. Next code uses a trick to set the play from HOC.
            #This fails:
            #rvec.record(self.nrnsec(ds/2 + i*ds)._ref_i_membrane)

            self.h('rvec.record(&i_membrane('+str(self.ds/2 + i*self.ds)+'))')
            
            self.nrnI.append(self.h.rvec)
                
        #Remove this section from the HOC pile
        self.h.pop_section()
    
    def currentrec(self):
        """
        Return recorded current as a numpy array
        """
                
        nt = len(self.nrnI[0])
        ns = self.nrnsec.nseg
        
        I = empty((nt,ns))
        
        for i in range(0,nt):
            for j in range(0,ns):
                I[i,j] = self.nrnI[j].x[i]
        
        return I

    def extvolt(self):
        """
        Return recorded voltages at the extremes as two numpy arrays
        """
        return np.array(self.nrnV0), np.array(self.nrnV1)
            
    def name(self):
        return self.nrnsec.name()

    def points(self):
        """
        Retrieve initial and final 3D coordinates 
        of the section from NEURON. Return is in meters.
        TODO: use centimeters for consistency
        """

        #Set as the currently accessed section
        self.nrnsec.push()
        
        #Obtain the index of the last segment end
        m = int(self.h.n3d())-1
        
        #We only care about the initial and last. 
        #NEURON uses microns
        p0 = array([self.h.x3d(0),self.h.y3d(0),self.h.z3d(0)])/1.0e6
        p1 = array([self.h.x3d(m),self.h.y3d(m),self.h.z3d(m)])/1.0e6
        
        #WARNING: Ignore Z
        #p0 = array([self.h.x3d(0),self.h.y3d(0),0.0])/1.0e6
        #p1 = array([self.h.x3d(m),self.h.y3d(m),0.0])/1.0e6
        
        #WARNING: Ignore Y and Z
        #p0 = array([self.h.x3d(0),0.0,0.0])/1.0e6
        #p1 = array([self.h.x3d(m),0.0,0.0])/1.0e6
        
        #Remove this section from the HOC pile
        self.h.pop_section()
        
        return p0,p1
    
    def segments(self):
        """
        Retrieve all the pt3d defined segments in an array
        where each row is a point and a diameter.
        Return is in millimeters.
        """

        #Set as the currently accessed section
        self.nrnsec.push()
        
        #Obtain the number of segments for this section
        m = int(self.h.n3d())
        
        s = np.empty([m,4])
        
        for i in range(m):
            s[i,:] = np.array([self.h.x3d(i),self.h.y3d(i),self.h.z3d(i),self.h.diam3d(i)])
        
        #Remove this section from the HOC pile
        self.h.pop_section()
        
        return s
    
    def nseg(self):
        """
        Retrieve number of segments 
        of the section from NEURON
        """
        return self.nrnsec.nseg

    def dx(self):
        """
        Calculate spatial differential
        of the section using NEURON parameters 
        and convert it from microns to cm
        """
        return 1e-4*self.nrnsec.L/self.nrnsec.nseg

    def L(self,val=None):
        """
        Set or get section's length in microns
        """
        
        if val == None:
          return self.nrnsec.L
        else:
          self.nrnsec.L = val
    
    def d(self):
        """
        Get section's diameter in microns as specified using
        pt3dadd.
        Several implementations have been tried trying to
        resolve the issue of which diameter should be
        considered when a section changes in diameter. Using
        two diameters (one for each end) is dangerous because
        if the next section uses the same diameter the total
        current injected for this section's end will be
        cancelled by the next. Returning either the start or
        end diameters can be guaranteed to work always as
        trees can be connected disorderly.
        """        
        #Diameter for standard neuron sections.
        #When the section has pt3dadd() segments, this
        #returns an erroneous value (500) on initialization
        #d = self.nrnsec.diam
        
        #Diameter from pt3dadd() defined sections

        #Set as the currently accessed section
        self.nrnsec.push()
        
        #Obtain the index of the last segment end
        m = int(self.h.n3d())-1
        
        d0 = self.h.diam3d(0)
        d1 = self.h.diam3d(m)
    
        d = (d0+d1)/2.0
    
        #Remove this section from the HOC pile
        self.h.pop_section()
        
        return d

    def setpoints(self,p0,p1):
        """
        Change the section coordinates. Use meters.
        """
        #Set as the currently accessed section
        self.nrnsec.push()
        
        #Obtain the index of the last segment end
        m = int(self.h.n3d())-1
        
        #import pdb; pdb.set_trace()
        
        #Set initial and final values
        self.h.pt3dchange(0,1e6*float(p0[0]),1e6*float(p0[1]),
            1e6*float(p0[2]),self.h.diam3d(0))
        self.h.pt3dchange(m,1e6*float(p1[0]),1e6*float(p1[1]),
            1e6*float(p1[2]),self.h.diam3d(m))
        
        #Remove this section from the HOC pile
        self.h.pop_section()
    
    def ra(self):
        """
        TODO: change to Ohms/cm
        Obtain the section's axial resistance in Ohm/m.
        For stimulation, the amplitude of the electric field 
        needs to be converted to the equivalent current injection. 
        This depends on the axial resistance ra, that can be 
        calculated from the section specific axial resistance Ra
        and diameter. 
        The specific resistance is specified in Ohm*cm and 
        the diameter in microns. The first is converted to Ohm*m
        and the second to meters.
        """
        #As the section could have variable diameter (when specified
        #with pt3dadd) two values could be returned for each end
        #d0,d1 = self.d()
        #ra0 = (4.0*(self.nrnsec.Ra/100.0))/(pi*(d0/1.0e6)**2)
        #ra1 = (4.0*(self.nrnsec.Ra/100.0))/(pi*(d1/1.0e6)**2)
        #return ra0,ra1
        
        return (4.0*(self.nrnsec.Ra/100.0))/(pi*(self.d()/1.0e6)**2)

    def rm(self):
        """
        Obtain the section's membrane resistance in Ohm*cm.
        This depends on the specific axial resistance Rm and
        diameter. The specific resistance is specified in Ohm*cm^2 and 
        the diameter in microns and is converted to cm.
        """
        return self.Rm()/(pi*(self.d()/1.0e4))

    def Ra(self,val=None):
        """
        Set or get section's axial resistance in Ohm*cm.
        """
        
        if val == None:
          return self.nrnsec.Ra
        else:
          self.nrnsec.Ra = val

    def Rm(self,val=None):
        """
        Set or get section's membrane resistance in Ohm*cm^2.
        """
        
        if val == None:
          return 1./self.nrnsec.g_pas
        else:
          self.nrnsec.g_pas = 1./val

class SectionList():
    """
    Object that holds a list of NEURON Sections out of 
    a NEURON SectionList object for easy handling within Python.
    SectionLists should usually represent a morphological part of a 
    neuron (e.g. branch, apical dendrite, axon)
    """
    def __init__(self,h,name):
        
        #Retrieve the section group by name as created by NEURON's
        #cell builder 
        h('obfunc returnobject(){ return '+name+' }')
        self.nrnlist = h.returnobject()
        self.list = []
        for nrnsec in self.nrnlist:
            sec = Section(h,nrnsec)
            self.list.append(sec)
            #DEBUG
            #print h.secname()
            
        #Recort the time vector of each simulation
        self.nrnT = h.Vector()
        self.nrnT.record(h._ref_t)
    
    def time(self):
        """
        Return a numpy array of the recorded time
        """
        return np.array(self.nrnT)
        
    def segments(self):
        """
        Return a list with the array for each section's
        segments
        """
        l = []
        for s in self.list:
            l.append(s.segments())
        return l

    def byname(self,name):
        """
        Return a section from the list with the given name
        """
        for s in self.list:
            if s.name() == name:
                return s
        return None

    def rotate(self,dtheta):
        #Rotation matrix
        R = array([[cos(dtheta),-sin(dtheta),0],
            [sin(dtheta), cos(dtheta),0],
            [          0,           0,1]])
    
        for s in self.list:
            p0,p1 = s.points()
            p0 = dot(R,p0)
            p1 = dot(R,p1)
            s.setpoints(p0,p1)

    def translate(self,dx,dy,dz):
        
        for s in self.list:
            p0,p1 = s.points()
            p0 = p0 + array([dx,dy,dz])
            p1 = p1 + array([dx,dy,dz])
            s.setpoints(p0,p1)
            
    def meanpoint(self):
      """
      For each section find its mean point.
      """
      #Create a list to store the mean points of each section
      X = empty((len(self.list),3))

      for i,s in enumerate(self.list):
        p0,p1 = s.points()
        #Calculate the mean point
        X[i] = (p1 + p0)/2.0
      
      return X

    def recordcurrent(self):       
        for s in self.list:
            s.recordcurrent()

    def Ra(self,val):
      for s in self.list:
        s.Ra(val)
    
    def Rm(self,val):
      for s in self.list:
        s.Rm(val)
