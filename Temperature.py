import numpy as np


def calc_Temp(dl, C_p, m_gas, m_liquido, T_0):
    g = 9.81  # m/s²

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

    m_total = m_gas + m_liquido

    # T_0 é a temperatura anterior, a inicial no poço é igual à 80°C.
    T = (
        T_inf
        - m_total * g * np.sin(theta_rad) / TEC
        - np.exp((-TEC * dl) / (m_total * C_p))
        * (T_inf - (m_total * g * np.sin(theta_rad) / TEC) - T_0)
    )

    return T
