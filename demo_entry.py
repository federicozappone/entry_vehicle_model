import numpy as np
import trimesh

from vehicle import Vehicle
from planet import Mars
from plot_utils import do_plot, do_3d_plot

from flow import Mach_vector, pres_coeff_max, pres_coeff_mod_newton, pres_from_Cp, \
    surface_force, surface_moment, aero_coeff


class Sphere(Vehicle):

    def __init__(self, mass, diameter, planet):
        Vehicle.__init__(self, mass, diameter, diameter, planet)

        # using spherical moments of inertia
        Im = 0.4 * self.mass * (self.L/2.0)**2
        self.I = [Im, Im, Im, 0.0, 0.0, 0.0]

    # constant aerodynamic coefficients
    def get_aero_coefficients(self, Ma, V, p_inf, rho_inf, alpha, beta):
        return 0.0, 0.46, 0.0, 0.0, 0.0, 0.0

    # external aero coefficients
    def get_input_aero_coefficients(self, Ma, V, p_inf, rho_inf, alpha, beta):
        return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0


class Orion(Vehicle):

    def __init__(self, planet):
        Vehicle.__init__(self, 9300, 5.03, 5.03, planet)

        self.mesh = trimesh.load("models/orion.stl")
        rotation = trimesh.transformations.rotation_matrix(0.5 * np.pi, [0, 0, 1])
        self.mesh.apply_transform(rotation)

        self.normals = np.array(self.mesh.facets_normal)
        self.areas = np.array(self.mesh.facets_area)
        self.origins = np.array(self.mesh.facets_origin)
        self.CG = np.array(self.mesh.centroid)
        self.A_ref = np.sum(self.areas)

        self.coeffs = None

        # using spherical moments of inertia
        Im = 0.4 * self.mass * (self.L/2.0)**2
        self.I = [Im, Im, Im, 0.0, 0.0, 0.0]

    def get_aero_coefficients(self, Ma, V_inf, p_inf, rho_inf, alpha, beta):
        if self.coeffs is None:
            M_vector = Mach_vector(M_inf=Ma, alpha=0.0, theta=0.0)
            Cp_max = pres_coeff_max(M=Ma, gamma_var=1.33)
            Cp, delta = pres_coeff_mod_newton(self.normals, M_vector, Cp_max)

            p = pres_from_Cp(Cp, p_inf, rho_inf, V_inf)
            F = surface_force(p, self.normals, self.areas) * -1.0
            M, L = surface_moment(F, self.origins, self.CG)

            self.coeffs = aero_coeff(F, M, self.A_ref, self.L, rho_inf, V_inf, 0.0, 0.0)

            # force moments coefficients not used for now
            self.coeffs[3] = 0.0
            self.coeffs[4] = 0.0
            self.coeffs[5] = 0.0

        return self.coeffs

    # external aero coefficients
    def get_input_aero_coefficients(self, Ma, V_inf, p_inf, rho_inf, alpha, beta):
        return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0


def demo():

    planet = Mars()
    vehicle = Orion(planet)

    # initial conditions (taken from MSL entry parameters)
    entry_interface = 125e3
    longitude = np.radians(137.42)
    latitude = np.radians(-4.49)
    entry_velocity = 5.8e3
    flight_path_angle = np.radians(-15.5)
    heading_angle = np.radians(0.0)
    omega_x = 0.0
    omega_y = 0.0
    omega_z = 0.0
    pitch = 0.0
    roll = 0.0
    yaw = 0.0

    # array of initial conditions
    x0 = [
            entry_interface + planet.radius, longitude, latitude,
            entry_velocity, flight_path_angle, heading_angle,
            omega_x, omega_y, omega_z, pitch, roll, yaw
         ]


    vehicle.set_initial_values(x0)

    # timestep
    dt = 0.1

    x = []
    t = []
    current_altitude = x0[0] - planet.radius

    # simulate until distance from ground is 10km
    while vehicle.r.successful() and current_altitude > 10e3:
        step = vehicle.step(dt)
        x.append(step)
        t.append(vehicle.r.t)

        current_altitude = step[0] - planet.radius
        print("current altitude:", current_altitude)

    x = np.array(x)
    x[:, 0] -= planet.radius

    do_plot(
        "velocity (km/s)", x[:, 3] / 1e3,
        "altitude (km)", x[:, 0] / 1e3,
        "trajectory", "velocity/altitude"
    )


if __name__ == "__main__":
    demo()
