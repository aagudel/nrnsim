NRNSIM USER GUIDE
=================

Versions:

- V 0.7 19.7.2013
- V 0.6 23.7.2012
- V 0.3 14.4.2010
- V 0.1 8.12.2009

NRNSIM or NeuRoN SIMulation is a Python tool that facilitates the execution of simulations with the Neuron environment (http://www.neuron.yale.edu). NRNSIM uses Neuron's Python interface to provide access to Neuron. NRNSIM however facilitates the simulation of neurons by providing easy to use Python objects. NRNSIM is designed to extend Neuron with Python easily. NRNSIM provides one example application and this is for the simulation of the extracellular stimulation of neurons. This application can be configured with semi-arbitrary descriptions of electric fields changing in space and time and it will configure Neuron in such a way, that the current neuron geometry will respond as if it was inside this field. The electric fields and the type of Neuron mechanism for representation of the field in Neuron can be chosen as a Python object. At the moment, stimulation with Neuron IClamps and the 'extracellular' mechanism can be configured.

TOOL COMPONENTS
===============

NRNSIM is composed by a collection of Python scripts. In the current version these scripts are:

- nsrun.py:  This is the initialization script. This should be executed for 
  every simulation. In this version the script should be manually 
  edited to load the respective 'application' scripts. This script 
  should always call nrnsim.py.

- nrnsim.py: The main script. This contains the basic functions most 
  'application' script should call, such as those to initialize 
  Neuron and to set the default time parameters. The main 
  functions are:

  - settime(newdt,newtmax): Set time parameters in milliseconds.
  
  - nrngeom(name): Load the Neuron geometry given (.hoc file)

  - nrninit(): Prepare Neuron for simulation

  - updatemech(): Should be called if the Neuron mechanisms changed

  - sim(): initiates the Neuron simulation.

- section.py: Provides classes to wrap Neuron sections and section list 
  providing a cleaner interface to obtain information about
  their parameters. The main functions are:

  - Section.points(): returns a pair of 3D points with the start 
    and end of a neuron section.

  - Section.segments(): Retrieve all the pt3d defined segments in an
    array where each row is a point and a diameter.

  - Section.dx(): Calculate the spatial differential of the section
    using the Neuron parameters of length and number of segments.
    
  - Section.setpoints(p0,p1): Change the section initial and
    end coordinates.

  - Section.ra(): Obtain the section's axial resistance in Ohm/m.

  - Section.rm(): Obtain the section's membrane resistance in 
    Ohm*cm.

  - Section.Ra(val): Set or get section's axial resistance in 
    Ohm*cm.

  - Section.Rm(val): Set or get section's membrane resistance in
    Ohm*cm^2.

- smech.py:  Implements the NEURON calls of different methods for extra and 
            some intracellular stimulation. Supported methods as Python 
            objects are:

  - ActivatingFunction: Implements the activation function [1] as a
    pair of stimulating IClamp objects in the two tips of each 
    section. Current values are calculated as in [2,3].

  - ExtracellularMechanism: Uses Neuron's extracellular mechanism 
    to set the extracellular potential for a set of sections. Each 
    Neuron section is divided in segments. For every segment on 
    every section, a potential value over time must be provided.

  - CurrentClamp:  A simple to use wrapper for Neuron's IClamp
    mechanism.
     
- fields.py: Provides objects for signals and fields. A signal is a variable
            depending on time. A field is a 3d variable depending on time
            defined for a set of points. Classes provided are:
    
  - StepSignal: Given a time vector and a start and end generates
    a time series in a vector that resembles a step function.

  - BipolarSignal: A square bipolar function where a start, 
    switching time and end time should be provided.

  - CosSignal: A cosine signal with and amplitude, period, start 
    of the signal and end of the signal.

  - RLCSignal: Emulates the electric field induced by the inductor 
    in a Resistor, Inductor, Capacitor circuit. It emulates the 
    electric field produced by a transcranial magnetic stimulation 
    coil.

  - Field: A general class to represent a field. Field objects can 
    take as input a signal object to scale them over time.

  - HomogField: produces a spatially homogeneous field but variable 
    in time field for a given time changing signal.


EXAMPLE APPLICATIONS
====================

One example application script is provided with NRNSIM: 

- multi.py:  Simulates multiple types of 1D neural geometries, under 
  geometrical transformations and multiple stimulation methods. 
  The file multi.py can be used as a template for other 
  applications. multi.py provides a graphical user interface to 
  change the magnitude of and duration of the stimulation field. 
  On execution of nsrun.py, parameters can be given in the 
  command line to load an specific hoc file and an extra Python 
  script. This script can modify the default amplitude, duration, 
  orientation of the field and timesteps.

TOOL USAGE
==========

NRNSIM requires that Neuron is installed with Python support and that Neuron is visible to the Python interpreter, this is explained in 
http://www.neuron.yale.edu/neuron/download

For execution the regular Python interpreter can be used but IPython is recommended. To run the default multi.py application with IPython < 0.11 use the command:

    ipython -pylab -i -c '%run -i nsrun.py file'

Where file is the root name of a set of files:
- file.hoc:     a Neuron geometry. Cell parameters can also be set.
- file.ses:     a Neuron session file.
- file.py:      a configuration file with the stimulation parameters.
- file.cin.hoc: alternatively, a file with hoc's init() function. This
                can be used to set arbitrary initial potentials 
                for different cells.
                         
IPython 0.11 and superior changed the interaction of IPython with QT gui (which is required by the tool). To use IPython then using a separate kernel is recommended to do this call instead:

    ipython console --pylab qt -i -c "%run -i nsrun.py file"

alternatively simple IPython can be used (no kernel mode) calling:

    ipython -i -c "%run -i nsrun.py file"

but in this case IPython's plotting routines cannot be used.

REFERENCES
==========

[1] Frank Rattay. Analysis of models for external stimulation 
    of axons. IEEE Transactions on Biomedical Engineering,
    33(10):974-977, 1986.

[2] S. S. Nagarajan, D. M. Durand, and E. N. Warman.
    Efects of induced electric filds on finite neuronal 
    structures: a simulation study. IEEE transactions on
    bio-medical engineering, 40(11):1175-1188, 1993.

[3] Assaf Rotem and Elisha Moses. Magnetic stimulation of
    one-dimensional neuronal cultures. Biophysical Journal,
    94(12):5065-5078, 2008.

COPYRIGHT
=========

Copyright Andres Agudelo-Toro (https://sites.google.com/site/aagudelotoro/). All source code can be used for scientific purposes. The author holds no liability.

