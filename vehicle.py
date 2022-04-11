import math
import numpy as np

from scipy.integrate import ode
from math import pi, sin, cos, tan, asin


class Vehicle:

    # for a generic entry "capsule" vehicle the length (L) is the base radius and the reference length (c) is equal to L
    def __init__(self, mass, L, c, planet, I=[1.0, 1.0, 1.0, 0.0, 0.0, 0.0]):
        self.mass = mass         # vehicle mass (kg)
        self.L = L               # length (m)
        self.c = c               # reference length (m)
        self.A = pi * (c**2)     # reference area (m2)
        self.I = I               # moments of inertia (kg-m2)
        self.planet = planet
        self.r = ode(self.odes)

    def set_initial_values(self, x0):
        self.r.set_initial_value(x0, 0.0).set_f_params(12.0)

    def step(self, dt):
        return self.r.integrate(self.r.t + dt)

    def get_aero_coefficients(self, Ma, V, p_inf, rho_inf, alpha, beta):
        pass

    def get_input_aero_coefficients(self, Ma, V, p_inf, rho_inf, alpha, beta):
        pass

    def odes(self, t, x, args):
        # radius, lon, lat, velocity, flight path angle, heading angle, angular velocity, pitch, roll, yaw
        r, theta, phi, V, gamma, psi, omega_x, omega_y, omega_z, pitch, roll, yaw = x

        omega = self.planet.omega # angular velocity of the planet

        m = self.mass # vehicle mass
        g = self.planet.get_gravity(r) # planet gravity

        # thrust vector, not used
        epsilon = 0.0 # first thrust vector angle
        zeta = 0.0    # second thrust vector angle
        T = 0.0       # thrust

        # moments of inertia
        I_xx, I_yy, I_zz, I_xy, I_yz, I_zx = self.I

        # geometric ralationships
        beta = asin( cos(gamma)*(cos(roll)*sin(yaw - psi) + sin(pitch)*sin(roll)*cos(yaw - psi)) - sin(gamma)*cos(pitch)*sin(roll) ) # sideslip angle
        alpha = asin( (cos(gamma)*(sin(pitch)*cos(roll)*cos(yaw - psi) - sin(roll)*sin(yaw - psi)) - sin(gamma)*cos(pitch)*cos(roll)) / cos(beta) ) # angle of attack
        sigma = asin( (cos(alpha)*sin(beta)*sin(pitch) - sin(alpha)*sin(beta)*cos(roll)*cos(pitch) + cos(beta)*sin(roll)*cos(pitch)) / cos(gamma) ) # bank angle

        altitude = r - self.planet.radius
        # atmospheric properties
        p, rho, a, mu = self.planet.atmosphere(t, theta, phi, altitude) # freestream pressure, density, speed of sound, viscosity

        # mach number
        Ma = V / a # freestream mach number

        # reynold numbers
        Re = (rho * V * self.L) / mu

        # aerodynamic coefficients
        C_L, C_D, C_S, C_N, C_A, C_M = self.get_aero_coefficients(Ma, V, p, rho, alpha, beta)

        # aerodynamic forces
        L = 0.5 * rho * C_L * self.A * (V**2) # lift
        D = 0.5 * rho * C_D * self.A * (V**2) # drag
        S = 0.5 * rho * C_S * self.A * (V**2) # side (not used)

        # force moments
        M_x = 0.5 * rho * C_A * self.A * self.L * (V**2)
        M_y = 0.5 * rho * C_M * self.A * self.L * (V**2)
        M_z = 0.5 * rho * C_N * self.A * self.L * (V**2)


        # kinematic equations
        r_dot = V * sin(gamma)
        theta_dot = (V * cos(gamma) * cos(psi))/(r * cos(phi))
        phi_dot = (V * cos(gamma) * sin(psi))/r

        pitch_dot = (omega_y * sin(roll)) + (omega_z * cos(roll))
        roll_dot = omega_x - (omega_y*cos(roll) - omega_z*sin(roll))*tan(pitch)
        yaw_dot = (1.0 / cos(pitch)) * (omega_y*cos(roll) - omega_z*sin(roll))


        # dynamical equations
        V_dot = (T/m)*(cos(zeta)*cos(epsilon)) -(D/m) -g*sin(gamma) \
                + r*(omega**2)*cos(phi)*(cos(phi)*sin(gamma) - sin(phi)*sin(psi)*cos(gamma))


        gamma_dot = (1.0 / V)*( (T/m)*(sin(zeta)*sin(sigma) + cos(zeta)*sin(epsilon)*cos(sigma)) \
                            + (L/m)*cos(sigma) -g*cos(gamma) \
                            + ((V**2)/(r))*cos(gamma) + 2.0*V*omega*cos(phi)*cos(psi) \
                            + r*(omega**2)*cos(phi)*(cos(phi)*cos(gamma) + sin(phi)*sin(psi)*sin(gamma)) )


        psi_dot = (1.0 / V)*( (1.0 / (m*cos(gamma)))*(T*(cos(zeta)*sin(epsilon)*sin(sigma) - sin(zeta)*cos(sigma)) + L*sin(sigma)) \
                            - ((V**2)/r)*cos(gamma)*cos(psi)*tan(phi) + 2.0*V*omega*(sin(psi)*cos(phi)*tan(gamma) - sin(phi)) \
                            - ((r*(omega**2))/cos(gamma))*sin(phi)*cos(phi)*cos(psi) )

        # rotational dynamics
        omega_x_dot = (1.0 / (I_xx*I_yy - (I_xy**2))) * (I_yy*M_x + I_xy*(M_y - (I_xx + I_yy - I_zz)*omega_x*omega_z) + ((I_xy**2) + (I_yy**2) - I_yy*I_zz)*omega_y*omega_z)
        omega_y_dot = (1.0 / (I_xx*I_yy - (I_xy**2))) * (I_xx*M_y + I_xy*(M_x - (I_xx + I_yy - I_zz)*omega_y*omega_z) + ((I_xy**2) + (I_yy**2) - I_xx*I_zz)*omega_x*omega_z)
        omega_z_dot = (M_z / I_zz) + ((I_xx - I_yy) / I_zz)*omega_x*omega_y + (I_xy / I_zz)*(omega_x**2 - omega_y**2)


        return [r_dot, theta_dot, phi_dot, V_dot, gamma_dot, psi_dot, omega_x_dot, omega_y_dot, omega_z_dot, pitch_dot, roll_dot, yaw_dot]
