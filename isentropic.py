from math import sqrt


def isentropic_relations(Ma, rho, gamma):
    Ma1 = sqrt( ((gamma - 1) * Ma**2 + 2.0) / (2.0 * gamma * Ma**2 - (gamma - 1.0)) )
    rho1 = ( ((gamma + 1.0) * Ma**2) / ((gamma - 1.0) * Ma**2 + 2.0) ) * rho

    return Ma1, rho1