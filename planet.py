from math import exp, sqrt
from mars_atmosphere import Mars_Atmosphere


class Planet:

    def __init__(self, mass, radius, omega):
        self.mass = mass
        self.radius = radius
        self.omega = omega # angular velocity (rad)

    def get_gravity(self, r):
        return self.mass * 6.67408e-11 / r**2


class Mars(Planet):

    def __init__(self):
        Planet.__init__(self, 6.39e23, 3389e3, 7.094834e-5)
        self.atmosphere = Mars_Atmosphere()