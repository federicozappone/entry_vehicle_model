import numpy as np


def Mach(a, V):
    """
    Calculates flow Mach number
    """

    Ma = V / a

    return Ma


def Reynolds(rho, U, L, mu):
    """
    Calculates flow Reynolds number
    """

    Re = (rho * U * L) / mu

    return Re


def speed_of_sound(gamma_var, R, T):
    """
    Calculates local speed of sound.
    Input variables:
        gamma_var   :   Ratio of specific heats
        R           :   Specific gas constant for given gas
        T           :   Gas temperature
    """


    a = np.sqrt(gamma_var * R * T)

    return a


def normal_shock_ratios(Ma_1, gamma_var):
    """
    Returns normal shock ratios for static and stagnation pressure and
    temperature, and density.  Also returns the Mach number following the
    shock (http://www.grc.nasa.gov/WWW/k-12/airplane/normal.html).
    Note that input variables are for flow UPSTREAM of shock, while returns
    are for the flow DOWNSTREAM of the shock.  Returned ratios are of the form:
    'Downstream condition' / 'Upstream condition' i.e.
    'Condition beyond shock' / ' Condition in front of shock'
    Input variables:
        Ma_1        :   Mach number upstream of shock
        gamma_var   :   Ratio of specific heats
        
    Returns:
        [0] : Static pressure ratio (p2/p1)
        [1] : Static temperature ratio (T2/T1)
        [2] : Total pressure ratio (p02/p01)
        [3] : Total temperature ratio (always 1.0)
        [4] : Density ratio (rho2/rho1)
        [5] : Post-shock Mach number (Ma2)
        [6] : Stagnation-static pressure ratio (p02 / p1)
    """

    p_ratio = ((2 * gamma_var * (Ma_1**2)) - (gamma_var - 1)) / (gamma_var + 1)

    T_ratio = (((2 * gamma_var * (Ma_1**2)) - (gamma_var - 1)) * \
        (((gamma_var - 1) * (Ma_1**2)) + 2)) / (((gamma_var + 1)**2) * (Ma_1**2))

    rho_ratio = ((gamma_var + 1) * (Ma_1**2)) / (((gamma_var - 1) * \
        (Ma_1**2)) + 2)

    p_0_ratio = ((((gamma_var + 1) * (Ma_1**2)) / (((gamma_var - 1) * \
        (Ma_1**2)) + 2))**(gamma_var / (gamma_var - 1))) * \
        (((gamma_var + 1) / ((2 * gamma_var * (Ma_1**2)) - \
        (gamma_var - 1)))**(1 / (gamma_var - 1)))

    T_0_ratio = 1.0

    Ma_2 = np.sqrt((((gamma_var - 1) * (Ma_1**2)) + 2) / ((2 * gamma_var * \
        (Ma_1**2)) - (gamma_var - 1)))
    
    p02_p1_ratio = ((((gamma_var + 1) * (Ma_1**2)) / 2)**(gamma_var / (gamma_var - 1))) / \
        ((((2 * gamma_var * (Ma_1**2)) / (gamma_var + 1)) - \
        ((gamma_var - 1) / (gamma_var + 1)))**(1 / (gamma_var - 1)))

    return [p_ratio, T_ratio, p_0_ratio, T_0_ratio, rho_ratio, Ma_2, p02_p1_ratio]


def shock_standoff(R_n, rho_inf, rho_s):
    """
    Approximates supersonic shock stand-off distance for the stagnation
    point of an obstacle with equivalent radius R_n
    Input variables:
        R_n     :   Obstacle radius
        rho_inf :   Freestream density
        rho_s   :   Density at point immediately behind shock
    """

    delta = R_n * (rho_inf / rho_s)

    return delta


def p_stag(**kwargs):
    """
    Calculates stagnation pressure based upon either of three input
    variable sets.  Optionally returns the ratio between stagnation
    and freestream pressure if no static term is supplied.
    First method:
        p_stag(rho = fluid density,
            V = fluid velocity,
            p = static pressure)
    Second method:
        p_stag(Ma = fluid Mach number,
            gamma_var = ratio of specific heats,
            p = static pressure)
    Return ratio:
        p_stag(Ma = fluid Mach number,
            gamma_var = ratio of specific heats)
    """

    if ('rho' in kwargs) and ('V' in kwargs) and ('p' in kwargs):
        p_0 = kwargs['p'] + (0.5 * kwargs['rho'] * (kwargs['V']**2))
    elif ('p' in kwargs) and ('gamma_var' in kwargs) and ('Ma' in kwargs):
        p_0 = kwargs['p'] * ((1 + (((kwargs['gamma_var'] - 1) / 2) * \
            (kwargs['Ma']**2)))**(kwargs['gamma_var'] / (kwargs['gamma_var'] - 1)))
    elif ('gamma_var' in kwargs) and ('Ma' in kwargs):
        p_0 = ((1 + (((kwargs['gamma_var'] - 1) / 2) * \
            (kwargs['Ma']**2)))**(kwargs['gamma_var'] / (kwargs['gamma_var'] - 1)))
    else:
        raise KeyError('Incorrect variable assignment')

    return p_0


def p_dyn(**kwargs):
    """
    Calculates dynamic pressure based upon either of two input
    variable sets.
    First method (incompressible flow only):
        p_dyn(rho = fluid density,
            V = fluid velocity)
    Second method (compressible flow only):
        p_dyn(Ma = fluid Mach number,
            gamma_var = ratio of specific heats,
            p = static pressure)
    """

    if ('rho' in kwargs) and ('V' in kwargs):
        q = 0.5 * kwargs['rho'] * (kwargs['V']**2)
    elif ('Ma' in kwargs) and ('p' in kwargs) and \
        (('gamma_var' in kwargs) or ('gamma' in kwargs)):
            q = 0.5 * kwargs['gamma'] * kwargs['p'] * (kwargs['Ma']**2)

    return q


def pres_coeff_max(M, gamma_var):
    """
    Calculates the maximum pressure coefficient on a surface following a normal
    shock at the stagnation point
    Input variables:
        M           :   Freestream Mach number (scalar, float)
        gamma_var   :   Specific heat ratio (scalar, float)
    """

    p01_inf_ratio = p_stag(Ma=M, gamma_var=gamma_var)
    _, _, p02_p01_ratio, _, _, _, _ = normal_shock_ratios(M, gamma_var)
    p02_inf_ratio = p02_p01_ratio * p01_inf_ratio

    Cp_max = (2 / (gamma_var * (M**2))) * (p02_inf_ratio - 1)

    return Cp_max


def pres_from_Cp(Cp, p_inf, rho_inf, V_inf):
    """
    Calculates pressure from freestream conditions based upon local pressure
    coefficient, Cp.
    Input variables:
        Cp          :   Local pressure coefficient
        p_inf       :   Freestream pressure
        rho_inf     :   Freestream density
        V_inf       :   Freestream velocity
    """

    p = p_inf + (Cp * 0.5 * rho_inf * (V_inf**2))

    return p


def Mach_vector(M_inf, alpha, theta):
    """
    Returns vector of components of Mach number based upon pitch and yaw
    angles of freestream flow direction.
    Input variables:
        M_inf   :   Freestream Mach number
        alpha   :   Freestream pitch angle
        theta   :   Freestream yaw angle
    """

    y = np.sin(alpha)
    z = np.sin(theta) * np.cos(alpha)
    x = np.cos(theta) * np.cos(alpha)

    M = M_inf * np.array([np.float(x), np.float(y), np.float(z)])

    return M


def pres_coeff_mod_newton(N, M_vector, Cp_max):
    """
    Calculates pressure coefficient along a surface in supersonic/hypersonic flow
    using the Modified Newtonian method.
    Input variables:
        N       :   Array of face normal vectors for the surface being assessed
        M       :   Freestream Mach number (array or list of components, 3 floats)
        Cp_max  :   Maximum pressure coefficient, evaluated at a stagnation
                    point behind a normal shock wave (scalar, float)
    """

    # Initalise variables
    len_N = len(N)
    delta = np.zeros([len_N])
    mag_N = np.zeros([len_N])
    mag_M = np.linalg.norm(M_vector)

    # Angle between two 3D vectors
    for index, value in enumerate(N):
        mag_N[index] = np.linalg.norm(N[index, :])
        delta[index] = np.arccos(np.dot(N[index], M_vector) / (mag_N[index] * mag_M))

    # Calculate Cp for all elements
    Cp = Cp_max * ((np.cos(delta))**2)

    # Apply element shielding assumption
    # i.e. if delta is greater than 90°, set Cp to 0
    for index, value in enumerate(Cp):
        if abs(delta[index]) > (np.pi / 2):
            Cp[index] = 0

    return Cp, delta


def surface_force(p, N, A):
    """
    Calculates the XYZ components of the force acting normal to a three
    dimensional surface element given the local pressure, element area, and
    element normal vector.
    Input variables:
        p   :   Pressure acting upon surface
        N   :   Surface normal vector
        A   :   Surface area
    """

    F_mag = p * A
    F_x = -F_mag * N[:, 0]
    F_y = -F_mag * N[:, 1]
    F_z = -F_mag * N[:, 2]
    F = np.array([F_x, F_y, F_z]).T

    return F


def surface_shear(t, S, A):
    """
    Calculates the XYZ components of the force acting tangential to a three
    dimensional surface element given the shear force, element area, and
    local shear vector.
    Input variables:
        t   :   Force acting tangential to surface
        S   :   Surface shear vector
        A   :   Surface area
    """

    T_x = -t * S[:, 0]
    T_y = -t * S[:, 1]
    T_z = -t * S[:, 2]
    T = np.array([T_x, T_y, T_z]).T

    return T


def surface_moment(F, C, CG):
    """
    Calculates moment(s) exerted about the centre of gravity of a structure by
    forces acting on individual surface panels.  Moment arm vectors are
    calculated about the centre of gravity based upon the surface element
    centres.
    Input variables:
        F   :   Force XYZ components (array or list, m X 3 floats)
        C   :   Coordinates of surface element centres (array or list, m X 3 floats)
        CG  :   Coordinates of centre of gravity (array or list, 3 floats)
    """

    # Calculate moment arms (XYZ vector components)
    L = C - CG

    # Split moment arms into axis directions
    L_x = L[:, 0]
    L_y = L[:, 1]
    L_z = L[:, 2]

    # Split forces into directional components
    F_x = F[:, 0]
    F_y = F[:, 1]
    F_z = F[:, 2]

    # Calculate moments
    M_x = (F_y * L_z) + (F_z * L_y)
    M_y = (F_x * L_z) + (F_z * L_x)
    M_z = (F_x * L_y) + (F_y * L_x)

    # Combine moments for export
    M = np.array([M_x, M_y, M_z]).T

    return [M, L]


def aero_coeff(F, M, A_ref, L_ref, rho_inf, V_inf, alpha, beta):
    """
    Calculates aerodynamic coefficients from the sum of forces and moments
    acting on a body and about its centre of gravity.
    Input variables:
        F       :   Forces acting on elements (array or list, m X 3 floats)
        M       :   Moments acting on elements about CoG (array or list, m X 3 floats)
        A_ref   :   Reference area
        L_ref   :   Reference length
        p_inf   :   Freestream pressure
        V_inf   :   Freestream velocity
    Returns (as a list):
        C_L     :   Lift coefficient (vertical force)
        C_D     :   Drag coefficient (longitudinal force)
        C_S     :   Side force coefficient (lateral force)
        C_N     :   Normal coefficient (yawing moment)
        C_A     :   Axial coefficient (rolling moment)
        C_M     :   Moment coefficient (pitching moment)
    """

    # Dynamic pressure
    q = p_dyn(rho=rho_inf, V=V_inf)

    # Sum of forces
    if np.shape(F) == (3,):
        F_x = F[0]
        F_y = F[1]
        F_z = F[2]
    else:
        F_x = np.sum(F[:, 0])
        F_y = np.sum(F[:, 1])
        F_z = np.sum(F[:, 2])

    # Sum of moments (in a given direction, NOT about an axis)
    # NB/ Moment components in [x, y, z] direction act about the [z&y, x&z, x&y] axes.
    if np.shape(M) == (3,):
        M_x = M[0]
        M_y = M[1]
        M_z = M[2]
    else:
        M_x = np.sum(M[:, 0])
        M_y = np.sum(M[:, 1])
        M_z = np.sum(M[:, 2])

    M_pitch = M_y
    M_roll = M_x
    M_yaw = M_z
    F_drag = F_x
    F_lift = F_z
    F_lat = F_y

    C_L = F_lift / (q * A_ref)
    C_D = F_drag / (q * A_ref)
    C_S = F_lat / (q * A_ref)
    C_N = M_yaw / (q * A_ref * L_ref)
    C_A = M_roll / (q * A_ref * L_ref)
    C_M = M_pitch / (q * A_ref * L_ref)

    return [C_L, C_D, C_S, C_N, C_A, C_M]
