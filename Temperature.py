import numpy as np


def calc_Temp(dl, m_total, T_0):
    g = 9.81  # m/s²
    SG = 0.72
    # Cálculo da temperatura externa
    if dl < 579.55:
        T_inf = 80 - 0.1277 * dl
        theta_rad = 10 * np.pi / 180  # rad
        TEC = 2  # W/mK
    elif dl < 867.49:
        T_inf = 6
        theta_rad = 15 * np.pi / 180  # rad
        TEC = 1  # W/mK
    else:
        T_inf = 6 + 0.01 * dl
        theta_rad = 90 * np.pi / 180  # rad
        TEC = 1  # W/mK

    C_p = ((2*10**(-3)*T_0 - 1.429)*SG + (2.67*10**(-3))*T_0 + 0.352)*10**(-3) #J/(Kg*°C)

    # T_0 é a temperatura anterior, a inicial no poço é igual à 80°C.
    T = (
        T_inf
        - m_total * g * np.sin(theta_rad) / TEC
        - np.exp((-TEC * dl) / (m_total * C_p))
        * (T_inf - (m_total * g * np.sin(theta_rad) / TEC) - T_0)
    )

    return T
