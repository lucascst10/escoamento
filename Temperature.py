import numpy as np

def calc_Temp(dl, m_mass, T_0, bsw, holdup):
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

    d15 = SG
    d_ot = d15 - 5.93e-4 * (T_0 - 15)
    C_po = ((2e-3 * T_0 - 1.429) * d_ot + (2.67e-3) * T_0 + 3.049) * 1000

    C_p_metano = 2.2537 * 1000
    C_p_etano = 1.7662 * 1000
    C_pg = 0.7 * C_p_metano + 0.3 * C_p_etano

    C_p_water = 4180
    C_pl = (1 - bsw) * C_po + bsw * C_p_water
    C_pm = holdup * C_pl + (1 - holdup) * C_pg

    Termo_A = T_inf - (m_mass * g * np.sin(theta_rad)) / TEC
    Exp_term = np.exp(-TEC * dl / (m_mass * C_pm))
    Termo_B = Exp_term * (Termo_A - T_0)

    T = Termo_A - Termo_B
    return T
