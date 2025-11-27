import numpy as np


def rho_w(S):
    rho_w = 62.368 + 0.438603 * S + 1.60074 * 10 ** (-3) * S**2  # lb/ft3
    rho_w = rho_w * 16.0185  # kg/m3
    return rho_w


def Rsw(Ppsi, Tf):
    # Coeficientes A
    A0 = 8.15839
    A1 = -6.12265e-2
    A2 = 1.91663e-4
    A3 = -2.1654e-7

    # Coeficientes B
    B0 = 1.01021e-2
    B1 = -7.44241e-5
    B2 = 3.05553e-7
    B3 = -2.94883e-10

    # Coeficientes C
    C0 = -9.02505
    C1 = 0.130237
    C2 = -8.53425e-4
    C3 = 2.34122e-6
    C4 = -2.37049e-9

    A = A0 + A1 * Tf + A2 * Tf**2 + A3 * Tf**3
    B = B0 + B1 * Tf + B2 * Tf**2 + B3 * Tf**3
    C = (C0 + C1 * Tf + C2 * Tf**2 + C3 * Tf**3 + C4 * Tf**4) * 1e-7

    Rsw = A + B * Ppsi + C * Ppsi**2

    return Rsw * 0.178108  # sm3/sm3


def Bww(Ppsi, Tf):
    delta_VwT = (
        -1.0001 * 10 ** (-2) + 1.333191 * 10 ** (-4) * Tf + 5.50654 * 10 ** (-7) * Tf**2
    )
    delta_VwP = (
        -1.95301 * 10 ** (-9) * Ppsi * Tf
        - 1.72834 * 10 ** (-13) * Ppsi**2 * Tf
        - 3.58922 * 10 ** (-7) * Ppsi
        - 2.25341 * 10 ** (-10) * Ppsi**2
    )

    Bw = (1 + delta_VwT) * (1 + delta_VwP)

    return Bw  # bbl/STB


def mu_w(P, T):
    A = 109.527
    B = -1.12166
    mu_w1 = A * T**B
    mu_w = (0.9994 + 4.0295 * 10 ** (-5) * P + 3.1062 * 10 ** (-9) * P**2) * mu_w1
    return mu_w
