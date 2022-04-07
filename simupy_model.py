import simupy_flight
import numpy as np
import math

from simupy.block_diagram import BlockDiagram
from simupy.block_diagram import DEFAULT_INTEGRATOR_OPTIONS, SimulationResult

from mars_atmosphere import Mars_Atmosphere
from plot_utils import do_plot, do_3d_plot


int_opts = DEFAULT_INTEGRATOR_OPTIONS.copy()
int_opts["nsteps"] = 100000


Ixx = 1.0 # kg-m2
Iyy = 1.0 # kg-m2
Izz = 1.0 # kg-m2
Ixy = 0.0 # kg-m2
Iyz = 0.0 # kg-m2
Izx = 0.0 # kg-m2

m = 3257.0 # kg

x = 0.0
y = 0.0
z = 0.0

nose_radius = 2.25 # m

S_A = np.pi * (nose_radius**2) # reference area
b_l = nose_radius
c_l = nose_radius
a_l = b_l

flight_path_angle = np.radians(15.5)

lat_ic = -5.49 * np.pi / 180
long_ic = 127.42 * np.pi / 180
h_ic = 125.0e3
V_N_ic = 0.0
V_E_ic = 5.6e3 * math.cos(flight_path_angle)
V_D_ic = 5.6e3 * math.sin(flight_path_angle)
psi_ic = 15.0 * np.pi / 180
theta_ic = 0.0 * np.pi / 180
phi_ic = 0.0 * np.pi / 180
omega_X_ic = 0.0 * np.pi / 180
omega_Y_ic = 0.0 * np.pi / 180
omega_Z_ic = 0.0 * np.pi / 180

mars_rotation_rate = 1.06e-7
mars_equitorial_radius = 3389.0e3
mars_gravitational_constant = 4.282837e13
mars_planetodetic = simupy_flight.Planetodetic(
                        a=mars_equitorial_radius,
                        omega_p=mars_rotation_rate,
                        f=0.0)

mars_atmosphere = Mars_Atmosphere()

planet = simupy_flight.Planet(
    gravity=simupy_flight.get_spherical_gravity(mars_gravitational_constant),
    winds=simupy_flight.get_constant_winds(),
    atmosphere=mars_atmosphere,
    planetodetics=mars_planetodetic,
)

vehicle = simupy_flight.Vehicle(
    base_aero_coeffs=simupy_flight.get_constant_aero(CD_b=0.47),
    m=m,
    I_xx=Ixx,
    I_yy=Iyy,
    I_zz=Izz,
    I_xy=Ixy,
    I_yz=Iyz,
    I_xz=Izx,
    x_com=x,
    y_com=y,
    z_com=z,
    x_mrc=x,
    y_mrc=y,
    z_mrc=z,
    S_A=S_A,
    a_l=a_l,
    b_l=b_l,
    c_l=c_l,
    d_l=b_l,
)

BD = BlockDiagram(planet, vehicle)
BD.connect(planet, vehicle, inputs=np.arange(planet.dim_output))
BD.connect(vehicle, planet, inputs=np.arange(vehicle.dim_output))

planet.initial_condition = planet.ic_from_planetodetic(
    long_ic, lat_ic, h_ic, V_N_ic, V_E_ic, V_D_ic, psi_ic, theta_ic, phi_ic
)
planet.initial_condition[-3:] = omega_X_ic, omega_Y_ic, omega_Z_ic

# simulate
res = BD.simulate(60*2, integrator_options=int_opts)

do_plot(
    "velocity (km/s)", res.y[:, simupy_flight.Planet.V_T_idx] / 1e3,
    "altitude (km)", res.y[:, simupy_flight.Planet.h_D_idx] / 1e3,
    "trajectory", "velocity/altitude"
)

do_3d_plot(
    "trajectory", 
    "lon (degrees)", np.degrees(res.y[:, simupy_flight.Planet.lamda_D_idx]),
    "lat (degrees)", np.degrees(res.y[:, simupy_flight.Planet.phi_D_idx]),
    "h (km)", res.y[:, simupy_flight.Planet.h_D_idx] / 1e3
)

