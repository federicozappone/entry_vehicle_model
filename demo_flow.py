import numpy as np
import trimesh
import time
import matplotlib.pyplot as plt

from flow import Mach_vector, pres_coeff_max, pres_coeff_mod_newton, pres_from_Cp, \
    surface_force, surface_moment, aero_coeff, Mach

from mars_atmosphere import Mars_Atmosphere


mesh = trimesh.load("models/orion.stl")

alpha = - np.radians(15) # 15 degrees angle of attack

# freestream on the x axis
rotation = trimesh.transformations.rotation_matrix(0.5 * np.pi + alpha, [0, 0, 1])
mesh.apply_transform(rotation)

atmosphere = Mars_Atmosphere()

# parameters
V_inf = 5800
altitude = 60e3

p_inf, rho_inf, a, mu = atmosphere(0.0, 0.0, 0.0, altitude)
Ma = Mach(a, V_inf)

normals = np.array(mesh.facets_normal)
areas = np.array(mesh.facets_area)
origins = np.array(mesh.facets_origin)
CG = np.array(mesh.centroid)

M_vector = Mach_vector(M_inf=Ma, alpha=np.radians(0.0), theta=np.radians(0.0))
Cp_max = pres_coeff_max(M=Ma, gamma_var=1.33)
Cp, delta = pres_coeff_mod_newton(normals, M_vector, Cp_max)

p = pres_from_Cp(Cp, p_inf, rho_inf, V_inf)

F = surface_force(p, normals, areas) * -1.0
M, L = surface_moment(F, origins, CG)

A_ref = np.sum(areas)
L_ref = 5.03


print("A_ref:", A_ref)
print("L_ref:", L_ref)


coeff = aero_coeff(F, M, A_ref, L_ref, rho_inf, V_inf, 0.0, 0.0)


coeff_labels = ["C_L", "C_D", "C_S", "C_N", "C_A", "C_M"]

for C, label in zip(coeff, coeff_labels):
    print(label, C)


pressure_colormap = plt.get_cmap("jet")(p / p.max())
pressure_colormap *= 255.0
pressure_colormap.astype(np.uint8)

print("p_max", p.max())

for facet, color in zip(mesh.facets, pressure_colormap):
    mesh.visual.face_colors[facet] = color

print("displaying pressure distrubution on mesh viewer")
mesh.show()
