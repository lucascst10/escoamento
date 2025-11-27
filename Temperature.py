import numpy as np

def calc_Temp(dl, m_total, T_0):
    g = 9.81 
    SG = 0.72

    if dl <= 579.55:
        T_inf = 80 - 0.1277 * dl
        theta_rad = 10 * np.pi / 180
        TEC = 2
    elif dl <= 867.49:
        T_inf = 6
        theta_rad = 15 * np.pi / 180
        TEC = 1
    else:
        T_inf = 6 + 0.01 * dl
        theta_rad = 90 * np.pi / 180
        TEC = 1

    C_p = ((2*10**(-3)*T_0 - 1.429)*SG + (2.67*10**(-3))*T_0 + 3.049) * 10**(3)

    T = (
        T_inf
        - m_total * g * np.sin(theta_rad) / TEC
        - np.exp((-TEC * dl) / (m_total * C_p))
        * (T_inf - (m_total * g * np.sin(theta_rad) / TEC) - T_0)
    )

    return T
