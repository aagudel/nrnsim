import nrnsim as ns
#import scanparam as sp
import multi as mu

#Run helper script
#See example usage in readme file
#Copyright Andres Agudelo-Toro (https://sites.google.com/site/aagudelotoro/)

param = None

def init():
    """
    Python's reload commands allow to reload a python 
    module while keeping variables on a persisten memory 
    workspace. When an application's code is changed and
    reloaded, some initalization commands do not need to 
    be re-run, for example some gui's, or slow data load 
    functions. 
    Function calls within this function are only run once 
    at the first execution.
    """
    global __persistent
    
    try:
        print __persistent
    except NameError:
        ns.init()
        #sp.init()
        mu.init()
        __persistent = 'Code reloaded'

def restart():
    global __persistent
    del __persistent
    init()

if __name__ == '__main__':
    
    #reload(ns)
    #reload(sp)
    #reload(mu)
    
    init()


