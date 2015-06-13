#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
An utility to present a simple gui panel and to allow easy adding of
buttons and text fields for interactive management of variables
in a python application. Can be used with matplotlib and ipython.

It is based on the idea first used in simulationpanel.py and
the tutorial code from http://zetcode.com/tutorials/pyqt4/

To use with ipython and matplotlib see:
https://github.com/ipython/ipython/blob/0.12.1/docs/examples/lib/gui-qt.py

Ideas for the future:
-Session saving: remember panel position, buttons added, etc.
-Re-callable routines: interactivelly edited, session saved, code you write in the gui itself
-Individual panels represent individual functionalities, modules
-Simulation steps: for visualization in time, a object interface that is 
called every timestep and something is done. To integrate neuron
like simulations or coil simulations some asynchronies must be introduced.
Could be considered a data source
-Plot displays that where data can be added and will be later recalled

Copyright Andres Agudelo-Toro (https://sites.google.com/site/aagudelotoro/)
"""

import sys
from PyQt4 import QtCore, QtGui

app = QtCore.QCoreApplication.instance()

if app is None:
    app = QtGui.QApplication([])
    print "Qt Application didn't exist"
else:
    print "Qt Application already existed"

class SimplePanel(QtGui.QMainWindow):
    
    def __init__(self, symb):
        
        super(SimplePanel, self).__init__()
        self.symb = symb
        self.init()
        
    def init(self):
        
        center = QtGui.QWidget(self)
        self.vbox = QtGui.QVBoxLayout(center)
        self.setCentralWidget(center)

        self.setWindowTitle('SimplePanel')
        self.statusBar()

        self.show()

    def addvar(self, name, vmin, vmax, step=None):
        """
        Add a variable to modify
        """
        hbox = QtGui.QHBoxLayout()
        label = QtGui.QLabel(name)
        spin = QtGui.QDoubleSpinBox()
        spin.setRange(vmin, vmax)
        spin.setAccelerated(True)
        if step is None: step = (vmax-vmin)/20.0
        spin.setSingleStep(step)

        hbox = QtGui.QHBoxLayout()
        hbox.addWidget(label)
        hbox.addWidget(spin)
        
        spin.setValue(self.symb[name])
        
        self.vbox.addLayout(hbox)
        
        spin.setProperty('name',name)
        spin.valueChanged.connect(self.chavar)
            
    def chavar(self, value):
        """
        Change variable through an event
        """
        name = str(self.sender().property('name').toString())
        self.symb[name] = value
        self.statusBar().showMessage(name + ' set')

    def addcom(self, command):
        """
        Add a python command, examples: 'calculate()', 'print var0'
        """
        button = QtGui.QPushButton(command)
        button.clicked.connect(self.docom)
        self.vbox.addWidget(button)
        
    def docom(self):
        self.statusBar().showMessage('...')
        command = str(self.sender().text())
        try:
            exec(command, self.symb)
        except:
            self.statusBar().showMessage('Error in '+command)
            raise
        else:
            self.statusBar().showMessage(command+' done')

if __name__ == '__main__':
    global test0
    test0 = 0
    sp = SimplePanel(globals())
    sp.move(100,100)
    sp.addvar('test0',0,100)
    sp.addcom('print test0')
    sys.exit(app.exec_())

