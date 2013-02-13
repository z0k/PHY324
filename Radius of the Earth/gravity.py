from numpy import *
from pylab import *
from scipy.optimize import leastsq


data_ascending = loadtxt("data_floor_ascending.txt", dtype='float', skiprows=1)
#Constants
G = 6.67 * 10 ** -11
MASS_OF_SUN = 2.0 * 10 ** 30  # in kilograms
DISTANCE_OF_SUN = 1.5 * 10 ** 8  # in kilometers
FLOOR_HEIGHT = 3.95  # in meters
FLOOR_MASS = 10. ** 6  # in kilograms
MG_CONSTANT = 0.10155
BURTON_G = 9.804253  # in meters per second square


#Data
time = data_ascending[3:, 0]
floor = data_ascending[3:, 1]  # Floor number
#g = abs(data_ascending[:, 2] - 576.5) * MG_CONSTANT / 100.
print floor

g = data_ascending[3:, 2]
print g
g = abs((0.18 * time) + g) * MG_CONSTANT / 100.
print g


#rel_floor = data_ascending[:, 0] - 8  # Floor relative to floor 8, negative
#rel_g = g[:] - g[8]  # Delta g with respect to floor 8.o
y = g / BURTON_G
basement_distance = floor * FLOOR_HEIGHT
print basement_distance


#del_g = zeros(len(g))
#del_g[0] = 0.
#for i in range(1, len(g)):
    #del_g[i] = abs(g[i] - g[i - 1])
#print del_g

p = zeros(2)
p[0] = 6.
p[1] = 1.


def peval(basement_distance, p):
    return -2. / p[0] * basement_distance + p[1]


def residuals(p, y, basement_distance):
    return y - peval(basement_distance, p)


p_final, cov_x, info, mesg, success = leastsq(residuals,
                                              p, args=(y, basement_distance),
                                              full_output=True)


title("Ascending Floors")
plot(basement_distance, y, linestyle='None', marker='o')
plot(basement_distance, peval(basement_distance, p_final))
xlabel("Floor Height from Basement")
ylabel("del g / BURTON G")
show()

#title("Distance from Reference Station")
#plot(rel_floor, rel_g, linestyle='None', marker='o')
#show()

print p_final[0]




#voltage = data[:, 0]
#current = data[:, 1]
#diameter = data[:, 2]

#voltage_error = ones(len(voltage)) * 0.05
#current_error = ones(len(current)) * 0.0005
#diameter_error = ones(len(diameter)) * 0.5
#voltage_root_error = 0.5 * (voltage) ** (-.5) * 0.05

##Convert to radius and convert to meters.
#radius = diameter / 200.
#radius_error = diameter_error / 200.
##Y value is computed as this expression to check linear relationship with
##a changing current as the independent variable.
#y = sqrt(voltage) / radius
##Error propogation in expression for y.
#y_error = (sqrt(voltage) / radius) * sqrt((voltage_root_error / voltage) ** 2 + (0.0025 / radius) ** 2)

##Fit parameters initialized with guess of 1.
#p = zeros(2)
#p[0] = 1.
#p[1] = 1.

##Constants; mu is described in the document. n is the number of turns in
##each coil, and R is the radius of the coils.
#mu = 4. * pi * 10. ** -7
#n = 130.
#R = 0.161
##Expression for the constant k in equation 8 from document.
#k = (mu * n / R) * (1 / sqrt(2.)) * (4. / 5.) ** (3. / 2.)

##Off axis correction ratio.
#Bz_ratio = zeros(len(radius))
#Bz_ratio = 1. - (radius ** 4) / ((R ** 4) * (0.6583 + 0.29 * (radius ** 2 / R ** 2)) ** 2)
#adj_radius_error = radius / Bz_ratio - radius


##The function representing equation 8 from the document.
#def eight(current, p):
    #return p[0] * k * (current - p[1])


#def residuals(p, y, current):
    #return y - eight(current, p)


#p_final, cov_x, info, mesg, success = leastsq(residuals, p,
#args=(y, current), full_output=True)

#y_final = eight(current, p_final)
#chi2 = sum((y - y_final) ** 2 / abs(y_error ** 2))
##Degrees of freedom.
#dof = len(current) - len(p_final)
##Reduced chi-squared value is calculated.
#print "RMS of residuals (sqrt(chisquared/dof))", sqrt(chi2 / dof)
##Output of fit statistics.
#for i in range(len(p_final)):
    #print "p[%d] =%12.4e +/- %.5e" % (i, p_final[i], sqrt(cov_x[i, i]))

#errorbar(current, y, y_error, fmt='r+')

##p[0] is squared in order to obtain the charge mass ratio.
#em_ratio = p_final[0] ** 2
##Error propogation for charge mass ratio.
#em_ratio_error = 2. * p[0] * sqrt(cov_x[0, 0])
#print "\nThe calculated value of (e/m) is%12.4e +/- %.5e." % (em_ratio, em_ratio_error)

#title("Charge Mass Ratio of an Electron")
#plot(current, y, linestyle='None', marker='o')
#plot(current, eight(current, p_final))
#ylabel(r'$\frac{\sqrt{V}}{r}$', fontsize=18)
#xlabel(r'$I$', fontsize=14)
#show()

##we tabulate some radius error - adjusted radius error data
#print "\n%7s%28s%31s" % ("Radius", "Radius Error", "Adjusted Radius Error")
#for i in range(len(radius)):
    #print "%10.3e%22.3e%22.3e" % (radius[i], radius_error[i], adj_radius_error[i])
