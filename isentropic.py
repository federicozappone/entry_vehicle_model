from math import sqrt


def isentropic_relations(Ma, p, rho, T, gamma):
    Ma1 = sqrt( ((gamma - 1) * Ma**2 + 2.0) / (2.0 * gamma * Ma**2 - (gamma - 1.0)) )
    p1 = ( (2.0* gamma * Ma**2 - (gamma - 1.0)) / (gamma + 1.0) ) * p
    rho1 = ( ((gamma + 1.0) * Ma**2) / ((gamma - 1.0) * Ma**2 + 2.0) ) * rho
    T1 = ( (2.0 * gamma * Ma**2 - (gamma - 1.0)) * ((gamma - 1.0) * Ma**2 + 2.0) / ((gamma + 1.0)**2 * Ma**2) ) * T

    return Ma1, p1, rho1, T1