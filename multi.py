"""
multi.py Simulate multiple types of 1D neural geometries, under
geometrical transformations and multiple stimulation methods.
Copyright Andres Agudelo-Toro (https://sites.google.com/site/aagudelotoro/)
"""
import numpy as np
import pylab as pl
import simplepanel as sp
import nrnsim as ns
import smech as sm
import fields as fi
import analytical as an
import argparse as ap

def update():
  """
  Set the stimulus, maximum amplitude of the field and
  update the mechanism accordingly
  """
  global signal
  signal = fi.StepSignal(1.0,ns.t,20*ns.dt,durat+20*ns.dt)
  ns.stim = fi.HomogField(signal,np.array(orient))
  if ns.mech is None:
    ns.mech = sm.ActivatingFunction(ns.h,ns.secs)
  ns.stim.amplitude(A)
  ns.updatemech()
    
  #DEBUG
  print 'Field and mechanism updated. Amplitude value is now',A

def hfield():
  """
  Set the field to be horizontal
  """
  global orient
  orient = [1.0,0.0,0.0]
  update()

def vfield():
  """
  Set the field to be vertical
  """
  global orient
  orient = [0.0,1.0,0.0]
  update()

def const():
  """
  Print some useful constants for all secs
  """
  for s in ns.secs.list:
    an.sigmai = 1000.0/s.Ra()
    an.Rm = s.Rm()
    an.L = 1e-4*s.L()
    an.d = 1e-4*s.d()
    print s.name(),an.constants()

def panel():
  global pan

  pan = sp.SimplePanel(globals())
  pan.move(0,400)

  pan.addvar('A',-2000.0,2000.0,step=50.0)
  pan.addvar('durat',0.1,durat,step=0.1)
  pan.addcom('update()')
  pan.addcom('hfield()')
  pan.addcom('vfield()')
  pan.addcom('const()')

def nrninit():
  """
  Neuron part of the initialization
  Next actions are performed:
  -Initialize nrnsim NEURON settings
  -Create a signal object to drive the stimulating field
  -Create extracellular field based on signal
  -Select stimulation mechanism
  """
  global signal, A, orient, durat, dt, tmax
  
  #Default values
  A = 1000.0
  durat = ns.tmax
  dt = 0.001
  tmax = 2.5
  orient = [1.0,0.0,0.0]
  
  parser = ap.ArgumentParser(description=__doc__)
  parser.add_argument('geom', action="store", help='Base name of the geometry to be used.')
  param = parser.parse_args()
  print 'Parameters to be used:', param
  
  geom = param.geom

  #Set this simulation geometry and time trajectory in milliseconds
  #geom = '1d_short_vs_long_3_12'; ns.settime(0.0005,2.0)
  #geom = '/home/agdtoro/lab-mpi/results/stim_simplistic_cable_3_12/multi_dendrite'; ns.settime(0.001,2.5)
  #geom = '/home/agdtoro/lab-mpi/results/stim_simplistic_cable_3_12/asymmetric'; ns.settime(0.001,2.5)
  #geom = '/home/agdtoro/lab-mpi/results/stim_simplistic_cable_3_12/dendrite_soma_axon'; ns.settime(0.0001,0.01)
  #geom = '/home/agdtoro/lab-mpi/results/stim_simplistic_cable_3_12/axon_dendrite'; ns.settime(0.0001,0.7)
  #geom = '/home/agdtoro/lab-mpi/results/soma_representation_4_12/nrnsim'; ns.settime(0.000010,0.001400)

  #geom = 'trac3'
  #geom = 'analytic_3_6_11'
  #geom = 'rotemaxon'
  #geom = 'rotemaxon_multi'

  #Load Neuron geometry  
  ns.nrngeom(geom)
  
  #Load parameter and configuration file
  execfile(geom+'.py',globals())
  ns.settime(dt, tmax)

  update()

  #signal = fi.CosSignal(1.0,0.240,ns.t,0.010,0.200)
  #signal = fi.StepSignal(1.0,ns.t,ns.dt,ns.tmax)
  #signal = fi.StepSignal(1.0,ns.t,20*ns.dt,durat+20*ns.dt)
  #signal = fi.StepSignal(1.0,ns.t,20*ns.dt,0.5)
  #signal = fi.BipolarSignal(1.0,ns.t,ns.dt,1.0,2.0)

  #ns.stim = fi.HomogField(signal,np.array(orient))
  #ns.mech = sm.ActivatingFunction(ns.h,ns.secs)
  
  #geom = '3d_better_synapse_1_12_2'
  #ns.nrninit(geom)
  #signal = fi.StepSignal(1.0,ns.t,ns.dt,ns.tmax)
  #ns.stim = fi.HomogField(signal)
  #ns.mech = sm.ActivatingFunction(ns.h,ns.secs)

  #geom = '3d_better_synapse_1_12_1'
  #ns.nrninit(geom)
  #signal = fi.StepSignal(1.0,ns.t,ns.dt,ns.tmax)
  #ns.stim = fi.HomogField(signal)
  #ns.mech = sm.ActivatingFunction(ns.h,ns.secs)

  #geom = '3d_synapse_12_11'
  #ns.nrninit(geom)
  #signal = fi.StepSignal(1.0,ns.t,ns.dt,ns.tmax)
  #ns.stim = fi.HomogField(signal)
  #ns.mech = sm.ActivatingFunction(ns.h,ns.secs)
  
  #geom = 'analytic_3_6_11'
  #ns.nrninit(geom)
  #signal = fi.StepSignal(1.0,ns.t,ns.dt,ns.tmax)
  #ns.stim = fi.HomogField(signal)
  #ns.mech = sm.ActivatingFunction(ns.h,ns.secs)
 
  #geom = 'trac1'
  #ns.nrninit(geom)
  #signal = fi.CosSignal(A,0.240,ns.t,0.010,0.200)
  #ns.stim = fi.HomogField(signal)
  #ns.mech = sm.ActivatingFunction(ns.h,ns.secs)
  
  #geom = 'trac2'
  #ns.nrninit(geom)
  #signal = fi.Signal(ns.t)
  #rawsignal = np.ones_like(ns.t); rawsignal[0] = 0.0; signal.set_raw(rawsignal)
  #ns.stim = fi.HomogField(signal)
  #ns.mech = sm.ActivatingFunction(ns.h,ns.secs)
  
  #tt = np.arange(0,1e-3,0.001e-3)
  #signal = fi.RLCSignal(tt, 0.010e-3, 0.900e-3)
  #pl.figure(10)
  #pl.plot(tt,signal.s)

  #Update the mechanism with the stim field for the first time
  
  #ns.stim.amplitude(A)
  #ns.updatemech()

  ns.nrnready(geom)

def init():
  
  #Amplitude of the electric field in V/m
  global A
  A = 1000.0

  nrninit()
  panel()

  print 'WARNING: check Y and/or Z are not ignored'
