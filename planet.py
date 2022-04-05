from math import exp, sqrt
from mars_atmosphere import Mars_Atmosphere


class Planet:

    def __init__(self, mass, radius, omega):
        self.mass = mass
        self.radius = radius
        self.omega = omega

    def get_atmosphere(self, r, longitude, latitude):
        pass

    def get_speed_of_sound(self, T):
        pass

    def get_gravity(self, r):
        return self.mass * 6.67408e-11 / r**2


class Mars(Planet):

    def __init__(self):
        Planet.__init__(self, 6.39e23, 3389e3, 1.06e-7)
        self.atmosphere = Mars_Atmosphere()

    def get_atmosphere(self, r, longitude, latitude):
        h = r - self.radius # altitude
        return self.atmosphere.get(h, longitude, latitude)

    def get_speed_of_sound(self, T):
        return sqrt(1.29 * 191.8 * T)

    # ratio of specific heats
    def gamma(self):
        return 1.33