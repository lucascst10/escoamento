import numpy as np


def driftflux(v_sg, v_m, v_sl, rho_g, rho_l, p_atm, P, g, d_h, sigma_l, theta):
    def wold_ghajar():
        C_0 = (v_sg / v_m) * (1 + ((v_sl / v_sg) ** ((rho_g / rho_l) ** 0.1)))
        k7 = (
            (g * d_h * sigma_l * (1 + np.cos(theta)) * (rho_l - rho_g)) / (rho_l**2)
        ) ** 0.25
        v_d = (3.583 * (1 + np.sin(theta)) * (p_atm / P)) * k7
        return C_0, v_d

    C_0, v_d = wold_ghajar()

    v_g = C_0 * v_m + v_d
    h_g = v_sg / (C_0 * v_m + v_d)

    return v_g, h_g

def homogeneo(holdup_l_ns, rho_l, rho_g, mu_l, mu_g):
    rho_m = holdup_l_ns * rho_l + (1-holdup_l_ns)* rho_g
    mu_m = holdup_l_ns * mu_l + (1-holdup_l_ns)*mu_g

    return rho_m, mu_m