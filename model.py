import math
import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import ode, odeint
from math import sin, cos, tan, pi


def odes(t, x, args):
    # altitude, lon, lat, flight path angle, heading angle
    r, theta, phi, V, gamma, psi = x

    omega = 1.06e-7 # angular velocity of the planet

    sigma = 0.0 # bank angle

    m = 3257.0 # vehicle mass
    g = 4.282837e13 / (r)**2 # planet gravity

    rho = 0.015 * np.exp(-(r-3389e3)/12e3) # density at altitude
    
    cd = 0.47 # cd for a sphere

    A = pi * (2.25**2) # reference area for a sphere of radius 2.25m

    LD = 0.32 # lift drag ratio
    D = 0.5*rho*cd*A*(V**2) # drag
    L = 0.0 # zero lift (ballistic entry)

    epsilon = 0.0 # first thrust vector angle
    zeta = 0.0 # second thrust vector angle
    T = 0.0 # thrust


    r_dot = V*sin(gamma)


    theta_dot = (V*cos(gamma)*cos(psi))/(r*cos(phi))


    phi_dot = (V*cos(gamma)*sin(psi))/r


    V_dot = (T/m)*(cos(zeta)*cos(epsilon)) -(D/m) -g*sin(gamma) \
            + r*(omega**2)*cos(phi)*(cos(phi)*sin(gamma) - sin(phi)*sin(psi)*cos(gamma))


    gamma_dot = (1.0/V)*( (T/m)*(sin(zeta)*sin(sigma) + cos(zeta)*sin(epsilon)*cos(sigma)) \
                        + (L/m)*cos(sigma) -g*cos(gamma) \
                        + ((V**2)/(r))*cos(gamma) + 2.0*V*omega*cos(phi)*cos(psi) \
                        + r*(omega**2)*cos(phi)*(cos(phi)*cos(gamma) + sin(phi)*sin(psi)*sin(gamma)) )


    psi_dot = (1.0/V)*( (1.0 / (m*cos(gamma)))*(T*(cos(zeta)*sin(epsilon)*sin(sigma) - sin(zeta)*cos(sigma)) + L*sin(sigma)) \
                        - ((V**2)/r)*cos(gamma)*cos(psi)*tan(phi) + 2.0*V*omega*(sin(psi)*cos(phi)*tan(gamma) - sin(phi)) \
                        - ((r*(omega**2))/cos(gamma))*sin(phi)*cos(phi)*cos(psi) )


    return [r_dot, theta_dot, phi_dot, V_dot, gamma_dot, psi_dot]



# altitude, lon, lat, flight path angle, heading angle
x0 = [140e3 + 3389e3, 0, 0, 5.6e3, np.radians(-15.5), np.radians(0.0)]

r = ode(odes, jac=None)
r.set_initial_value(x0, 0.0).set_f_params(6.0)

t1 = 5 * 60
dt = 0.1

x = []
current_altitude = 140e3 + 3389e3

while r.successful() and current_altitude > 3389e3:
    step = r.integrate(r.t + dt)
    x.append(step)
    current_altitude = step[0]
    print(r.t + dt, step)

x = np.array(x)
x[:, 0] -= 3389e3

"""
t = np.linspace(0, 3*60, 100)
x = odeint(odes, x0, t)

for i in range(100):
    print(x[:, 3][i], x[:, 0][i], np.degrees(x[:, 4][i]))

"""

def do_plot(xlabel, x, ylabel, y, label, title):
    # Basic utility function to simplify plotting
    plt.plot(x, y, label=label)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.xlim([0, 6.0])
    plt.ylabel(ylabel)
    plt.ylim([0, 140.0])
    plt.legend()
    plt.show()


f5 = plt.figure(figsize=(20, 12))

do_plot(
"velocity (km/s)", x[:, 3] / 1e3,
"altitude (km)", x[:, 0] / 1e3,
"test", "velocity/altitude"
)