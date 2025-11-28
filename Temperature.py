import numpy as np


def calc_Temp(dl, m_mass, T_0):
    g = 9.81
    SG = 0.72

    if dl <= 579.55:
        T_inf = 80 - 0.1277 * dl
        theta_deg = 10
        TEC = 2
    elif dl <= 867.49:
        T_inf = 6
        theta_deg = 15
        TEC = 1
    else:
        T_inf = 6 + 0.01 * dl
        theta_deg = 90
        TEC = 1

    theta_rad = np.deg2rad(theta_deg)

    C_p = ((2e-3 * T_0 - 1.429) * SG + (2.67e-3) * T_0 + 3.049) * 1000

    Termo_A = T_inf - (m_mass * g * np.sin(theta_rad)) / TEC

    Exp_term = np.exp(-TEC * dl / (m_mass * C_p))

    Termo_B = Exp_term * (Termo_A - T_0)

    T = Termo_A - Termo_B

    return T
