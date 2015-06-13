#Scan parameter module
#Copyright Andres Agudelo-Toro (https://sites.google.com/site/aagudelotoro/)

import pylab as pl
import simulationpanel
import nrnsim as ns
import smech
from scipy.optimize import leastsq

def panel():
  global Ras,Rms,Ie,Ies
  global tstart
  global pan,t0,t1

  #Pulse start time
  tstart = 0.1

  #Ras = [1.,500.,20.] #Ohm*cm
  #Rms = [200.,10000.,200.] #Ohm*cm^2
  
  Ras = [5.0, 500.0, 50.0] #Ohm*cm
  Rms = [200.0, 10000.0, 1000.0] #Ohm*cm^2
  
  Ie = 4. #nA
  Ies = pl.c_[pl.ones_like(ns.t)]
  Ies[:pl.ceil(tstart/ns.dt)] = 0.

  t0 = 0.101
  t1 = 0.6

  pan = simulationpanel.SimulationPanel()
  pan.Move((640,600))
  pan.setdict(globals())
  pan.addcommand('sim()')
  pan.addcommand('plottao()')
  pan.addvar('Ras')
  pan.addvar('Rms')
  pan.addvar('Ie')
  pan.addvar('t0')
  pan.addvar('t1')

def sim():
  global V,Vlin,tao_e,Rar,Rmr
  Rar = pl.arange(Ras[0],Ras[1],Ras[2])
  Rmr = pl.arange(Rms[0],Rms[1],Rms[2])
  
  ns.mech.setcurrent(Ie*Ies,ns.dt)
  
  li = len(Rmr)
  lj = len(Rar)
  tao_e = pl.empty((li,lj))
  tao_l = pl.empty((li,lj))
  tao_n = pl.empty((li,lj))
  
  for i in range(li):
    for j in range(lj):
      
      #Special conditions
      if Rar[j] < 10.:
        sec.L(15000.)
      else:
        sec.L(7000.)
        
      if Rmr[i] > 5000.:
        ns.h.tstop = 50.
      else:
        ns.h.tstop = 20.
      
      sec.Rm(Rmr[i])
      sec.Ra(Rar[j])
      print Rmr[i],Rar[j]
      ns.sim()
    
      #Obtain voltage, steady state voltage, normalize and 
      #get logarithmic values
      t = ns.t
      Vinf = sec.nrnV0[-1]
      V = 1 - pl.array(sec.nrnV0)[:-1]/Vinf
      Vlin = pl.log(V)
      print Vinf
      
      #Estimate the time constant finding the 
      #point at witch the voltage reaches the
      #value 1/e
      nz, = pl.nonzero(V>(1/pl.e))
      #The time where V ~ 1/e is the point 
      #right after the last nz
      tao_e[i,j] = t[nz[-1]+1] - tstart
      print 'tao_e',tao_e[i,j]
        
      #Define least squares data interval and
      #make the pulse starting time to be zero
      i0 = int(t0/ns.dt)
      i1 = int(t1/ns.dt)
      t01 = t[:i1-i0]
      V01 = V[i0:i1]
      Vlin01 = Vlin[i0:i1]
      
      #Linear least squares
      A = pl.c_[t01,pl.ones_like(t01)]
      m, c = pl.lstsq(A, Vlin01.copy())[0]
      tao_l[i,j] = -1./m - tstart
      print 'tao_l',tao_l[i,j],'(',m, c, pl.exp(c),')'
    
      #Parametric function: v is the parameter vector and
      #x the independent varible
      fp = lambda p, t: p[0]*pl.exp(p[1]*t)
      #fp = lambda p, t: p[0]*pl.exp(p[1]*t) + p[2]*pl.exp(p[3]*t)
      #fp = lambda p, t: pl.exp(p[0]*t)

      #Error function
      e = lambda p, t, V: (fp(p,t)-V)

      #Initial parameter guess
      p0 = [1., -5.]
      #p0 = [1., -5., 1., -1.]
      #p0 = [-5.]

      #Fitting
      p, success = leastsq(e, p0, args=(t01,V01), maxfev=10000)
      
      tao_n[i,j] = -1./p[1] - tstart
      print 'tao_n',tao_n[i,j],'(',p,success,')'
  
  """
  f = pl.figure(0)
  pl.clf()
  a = f.add_subplot(111)
  a.plot(t-t0,Vlin)
  a.plot([t01[0], t01[-1]],[Vlin01[0], Vlin01[-1]],'kx')
  a.plot(t01,m*t01 + c)
  pl.title('linear least squares')
  a.set_xlabel('t [ms]')
  a.set_ylabel('log(Vm)')
  pl.legend(('data','lsq interv.','linear'))

  f = pl.figure(1)
  pl.clf()
  a = f.add_subplot(111)
  a.plot(t-t0,V)
  a.plot([t01[0], t01[-1]],[V01[0], V01[-1]],'kx')
  a.plot(t,pl.exp(m*t + c))
  #a.plot(t,fp([p],t))
  a.plot(t,fp(p,t))
  pl.title('linear vs non-linear lsq')
  a.set_xlabel('t [ms]')
  a.set_ylabel('Vm [mV]')
  pl.legend(('data','lsq interv.','linear','non-linear'))
  """

def plottao():
  f = pl.figure(2)
  pl.clf()
  pl.pcolor(Rar,Rmr,tao_e)
  #pl.contourf(Rar,Rmr,tao_e,10)
  cb = pl.colorbar()
  pl.contour(Rar,Rmr,tao_e,10,colors='r')
  pl.ylim(Rms[0],Rms[1])
  cb.set_label('tao_e at x = 0 [ms]')
  pl.xlabel('Ra [Ohm*cm]')
  pl.ylabel('Rm [Ohm*cm^2]')
 
  pl.twinx()
  #Show real tao equivalent values right
  pl.ylim(Rms[0]*1e-3,Rms[1]*1e-3)  
  pl.ylabel('tao [ms]')
  
def nrninit():
  global sec
  geom = 'infaxon'
  ns.nrninit(geom)
  
  #Prepare stimulator
  ns.mech = smech.CurrentClamp(
    ns.h,ns.sections,[0.0])

  sec = ns.sections.list[0]

  ns.loadsession(geom)

def init():
        
  ns.settime(0.005,50.) #milliseconds
  panel()
  nrninit()
