import numpy as np

from vehicle import Vehicle
from planet import Mars
from plot_utils import do_plot, do_3d_plot


class GEV(Vehicle):

    def __init__(self, mass, L, c, planet):
        Vehicle.__init__(self, mass, L, c, planet)

    # constant aerodynamic coefficients for the demo
    def get_aero_coefficients(self, Ma, Re, alpha, beta):
        A = np.pi * (2.25**2)
        BC = 146
        Cd = 3257.0 / (A * BC)
        return 0.0, Cd, 0.0, 0.0, 0.0

    # external aero coefficients
    def get_input_aero_coefficients(self, Ma, Re, alpha, beta):
        return 0.0, 0.0, 0.0, 0.0, 0.0


def ballistic_entry_demo():

    planet = Mars()
    vehicle = GEV(3257.0, 2.25, 2.25, planet)

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

    x = np.array(x)
    x[:, 0] -= planet.radius

    do_plot(
        "velocity (km/s)", x[:, 3] / 1e3,
        "altitude (km)", x[:, 0] / 1e3,
        "trajectory", "velocity/altitude"
    )

    do_plot(
        "velocity (km/s)", x[:, 3] / 1e3,
        "time (s)", t,
        "trajectory", "longitude/latitude"
    )


    do_3d_plot("trajectory", "lon", np.degrees(x[:, 1]), "lat", np.degrees(x[:, 2]), "h", x[:, 0])


if __name__ == "__main__":
    ballistic_entry_demo()
