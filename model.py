import math
import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import ode, odeint
from math import sin, cos, tan, pi, abs


class Parachute:

    def __init__(self, diameter, cd):
        self.diameter = diameter
        self.reference_area = 0.25 * pi * diameter
        self.cd = cd
        # for supersonic inflation of viking type DGB parachute
        self.inflation_time = 0.02 * self.diameter
        self.deployment_time = 0.0
        self.deployment_initial_speed = 0.0
        self.n = 2.0 # inflation curve exponent, 2.0 for DGB parachutes

        self.is_deployed = False

    def deploy(self, t, V):
        if self.is_deployed is False:
            self.deployment_time = t
            self.deployment_initial_speed = V
            self.is_deployed = True

    def inflation(self, t):
        tf = self.deployment_time + self.inflation_time
        return abs(((t / self.deployment_time) / (tf / self.deployment_time))**self.n)

    def get_drag(self, rho, V, t):
        if self.is_deployed == True:
            return 0.5 * rho * self.cd * self.reference_area * (V**2) * self.inflation(t)
        else:
            return 0.0


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

    parachute_cd = 0.8 # from wind tunnel simulation of MSL DGB parachute
    parachute_diameter = 21.5 # MSL parachute diameter

    parachute = Parachute(parachute_diameter, parachute_cd)

    if (r - 3389e3) < 10e3 and parachute.is_deployed == False: # open parachute when altitude is 10km
        parachute.deploy(t, V)

    parachute_drag = parachute.get_drag(rho, V, t)

    LD = 0.24 # lift drag ratio
    D = 0.5 * rho * cd * A * (V**2) + parachute_drag # drag

    if parachute_drag == 0.0:
        L = LD * D # lifting entry
    else:
        L = 0.0 # zero lift when parachute opens since angle of attack is zero

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

# Initial conditions

entry_interface = 140e3
mars_radius = 3389e3
longitude = np.radians(137.42)
latitude = np.radians(-4.49)
entry_velocity = 5.6e3
flight_path_angle = np.radians(-15.5)
heading_angle = np.radians(0.0)


x0 = [entry_interface + mars_radius, longitude, latitude, entry_velocity, flight_path_angle, heading_angle]


r = ode(odes)
r.set_initial_value(x0, 0.0).set_f_params(6.0)

dt = 0.1

x = []
current_altitude = 140e3

while r.successful() and current_altitude > 1.8e3:
    step = r.integrate(r.t + dt)
    x.append(step)
    current_altitude = step[0] - mars_radius

    #print(r.t + dt, step)

x = np.array(x)
x[:, 0] -= mars_radius


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


f5 = plt.figure(figsize=(10, 6))

do_plot(
    "velocity (km/s)", x[:, 3] / 1e3,
    "altitude (km)", x[:, 0] / 1e3,
    "trajectory", "velocity/altitude"
)

def do_3d_plot(title, xlabel, x, ylabel, y, zlabel, z):
    fig = plt.figure(figsize=(10, 6))
    ax = plt.axes(projection="3d")

    plt.title(title)

    ax.set_xlabel(xlabel)
    #plt.xlim([0, 6.0])
    ax.set_ylabel(ylabel)
    #plt.ylim([0, 140.0])
    ax.set_zlabel(zlabel)

    ax.plot3D(x, y, z, "gray")

    # Data for three-dimensional scattered points
    """
    zdata = 15 * np.random.random(100)
    xdata = np.sin(zdata) + 0.1 * np.random.randn(100)
    ydata = np.cos(zdata) + 0.1 * np.random.randn(100)

    ax.scatter3D(xdata, ydata, zdata, c=zdata, cmap="Greens")
    """

    plt.show()

"""
# Data for a three-dimensional line
zline = np.linspace(0, 15, 1000)
xline = np.sin(zline)
yline = np.cos(zline)

do_3d_plot("x", xline, "y", yline, "z", zline)
"""

def to_cart(r, theta, phi):
    x = r*sin(phi)*cos(theta)
    y = r*sin(phi)*cos(theta)
    z = r*cos(phi)
    return [x, y, z]

do_3d_plot("trajectory", "lon", np.degrees(x[:, 1]), "lat", np.degrees(x[:, 2]), "z", x[:, 0])
x = x[x[:, 0] < 25e3]
do_3d_plot("trajectory after parachute deployment", "lon", np.degrees(x[:, 1]), "lat", np.degrees(x[:, 2]), "z", x[:, 0])
