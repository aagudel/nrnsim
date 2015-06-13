"""
nrnsim.py Provide common routines for interaction with NEURON
Copyright Andres Agudelo-Toro (https://sites.google.com/site/aagudelotoro/)
"""

import numpy as np
import simplepanel as sp
import section
import os.path as pt

"""
def plot_sect():
  p = empty((len(secs.list),3))
  for i,s in enumerate(secs.list):
    p[i] = s.points()[0]
  figure(id(secs))
  plot(p[:,0],p[:,2])
"""
    
def updatemech():
  """
  For each section find its average point in space.
  Calculate the stimulating field at each of these points.
  Inform the stimulation mechanism about the new field values.
  """        
  X = secs.meanpoint()    

  E = stim.calculate_field(X[:,0],X[:,1],X[:,2])
  #DEBUG
  print 'The field was calculated for %d points'%X.shape[0]
  #print X

  print dt
  print 'Setting stimuli ...'
  mech.setstimuli(E,dt)
  #DEBUG
  print 'Mechanism\'s field updated'
  #print E

def sim():
  """
  Run the simulation
  """
  h.init()
  h.run(tmax)

  #f = np.figure(1)
  #a1 = f.add_subplot(111)
  #a1.plot(t,np.array(secs.list[0].nrnV0)[:-1],'g')
  #a1.set_xlabel('t [us]')
  #a1.set_ylabel('Vm [mV]',color='g')
  
  #a2 = a1.twinx()
  #a2.plot(1e6*t,array(sm.nrnvec[0][0]),'b')
  #a2.plot(1e6*t,-stimfield.E[:,0,0])
  #a2.set_ylabel('dEx/dx [V/m^2]',color='b')    
    
def nrninit():
  """
  Prepare NEURON for simulation
  """
  from neuron import h
  from neuron import gui
  global h,mech,stim

  #Fillable mechanism
  mech = None
  
  #Fillable stimulating field
  stim = None

def nrngeom(name):
  
  global secs
  
  #Load geometry
  h.load_file(name+'.hoc')
  secs = section.SectionList(h,'all')

def nrnready(name):
  """
  Call some final Neuron routines that couldn't
  be called at initialization.
  """

  #Is there a session file?, load it
  if pt.isfile(name+'.ses'):
    h.load_file(name+'.ses')
  else:
    print 'Session could not be loaded. Loading default'
    h.load_file('default.ses')

  #Is there a custom init file?, load it
  if pt.isfile(name+'.cin.hoc'):
    h.load_file(name+'.cin.hoc')
  else:
    print 'No custom init specified'
  
  #The session usually changes the time parameters.
  #Make sure they are set again.
  nrnsettime()
  
def nrnsettime():
  h.dt = dt
  h.steps_per_ms = 1.0/h.dt
  h.tstop = tmax
  h.v_init = 0.0

def settime(newdt,newtmax):
  """
  Set time parameters in milliseconds.
  """
  global dt,tmax,t
  
  dt = newdt
  tmax = newtmax
  t = np.arange(0,tmax,dt)

  nrnsettime()

def panel():
  """
  Provide basic functionality through a panel
  """
  global pan
  pan = sp.SimplePanel(globals())
  pan.move(0,200)
  pan.addcom('sim()')

def init():
  """
  Start up setting a default global time trajectory
  in milliseconds and initialize GUI
  """
  nrninit()
  settime(0.010, 6.0)
  panel()
  
  """
  #Enable ipython event loop support in IPython 0.12
  #See: https://github.com/ipython/ipython/blob/0.12.1/docs/examples/lib/gui-qt.py
  try:
      from IPython.lib.guisupport import start_event_loop_qt4
      start_event_loop_qt4(sp.app)
  except ImportError:
      print "Import error"
      #sp.app.exec_()
  """


