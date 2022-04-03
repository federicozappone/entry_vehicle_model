from math import exp, sqrt


class Planet:

    def __init__(self, mass, radius, omega):
        self.mass = mass
        self.radius = radius
        self.omega = omega

    def get_atmosphere(self, r):
        pass

    def get_speed_of_sound(self, T):
        pass

    def get_gravity(self, r):
        return self.mass * 6.67408e-11 / self.radius**2


class Mars(Planet):

    def __init__(self):
        Planet.__init__(self, 6.39e23, 3389e3, 1.06e-7)

    # mars atmospheric model for h > 7000m
    def get_atmosphere(self, r):
        h = r - self.radius                  # altitude
        T = (-23.4 - 0.00222 * h) + 273.15   # temperature
        p = 0.699 * exp(-0.00009 * h)        # pressure
        rho = p / (0.1921 * T)               # density at altitude
        return T, p, rho

    def get_speed_of_sound(self, T):
        return sqrt(1.29 * 191.8 * T)

    # ratio of specific heats
    def gamma(self):
        return 1.33