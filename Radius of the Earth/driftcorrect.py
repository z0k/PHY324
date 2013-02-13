#Electron Spin Resonance
#Iris Ma and Rohit Batra
#Tuesday Oct. 26, 2010


from scipy.optimize import leastsq
from numpy import *
from pylab import *

#Loading and Defining the Data                
data=loadtxt("dgdt3.txt")
dt=data[:,0] #frequency originally in mHz converted to Hz.
dg=data[:,1] #Current. Divide by 2 because the aquired current was peak to peak.

#Constants

#Equations

#Initial Parameter Estimates
pini=zeros(2)
pini[0]=0.1
pini[1]=0

#Residuals and Line of best Fit Functions
def residuals(p,dg,x_val):
    err = dg-peval(x_val,p)
    return err

def peval(x_val,p):
    return p[0]*x_val+p[1]

plsq = leastsq(residuals, pini, args=(dg, dt))

#Plotting
title('Calculating the Drift')
xlabel('Delta Time')
ylabel('Delta Gravity')
grid('on')
plot (dt,dg, label='Experminetal', lw=2)
plot (dt,peval(dt,plsq[0]), "--", label='Best fit', lw=3)
legend(('Experimental', 'Best Fit'), loc='upper left', shadow='true')
show()
