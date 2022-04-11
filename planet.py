from math import exp, sqrt
from mars_atmosphere import Mars_Atmosphere


class Planet:

    def __init__(self, mass, radius, omega, atmosphere=None):
        self.mass = mass # kg
        self.radius = radius # m
        self.omega = omega # angular velocity (rad/s)
        self.atmosphere = atmosphere

    def get_gravity(self, r):
        return self.mass * 6.67408e-11 / r**2


class Mars(Planet):

    def __init__(self):
        Planet.__init__(self, 6.39e23, 3389e3, 7.094834e-5, Mars_Atmosphere())
