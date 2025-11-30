import numpy as np


def driftflux(
    v_sg, v_m, v_sl, rho_g, rho_l, p_atm, P, g, d_h, sigma_l, theta, mu_l, mu_g
):
    def Bendiksen():
        # print(f"v_sg:{v_sg,}, vm:{v_m}, v_sl:{v_sl}, v_sg:{v_sg}, rho_g:{rho_g}, rho_l:{rho_l}")
        Fr = v_m / np.sqrt(g * d_h)
        if Fr < 3.5:
            C_0 = 1.05 + 0.15 * np.sin(theta)
            v_d = np.sqrt(g * d_h) * ((0.35 * np.sin(np.deg2rad(theta))) + 0.54 * np.cos(np.deg2rad(theta)))
        else:
            C_0 = 1.2
            v_d = 0.35 * np.sqrt(g * d_h) * np.sin(theta)

        return C_0, v_d

    C_0, v_d = Bendiksen()

    v_g = C_0 * v_m + v_d
    h_g = v_sg / (C_0 * v_m + v_d)
    # print(f"h_g:{h_g}")
    # print(f"rho_l:{rho_l}, rho_g:{rho_g}, mu_l:{mu_l}, mu_g:{mu_g}")

    rho_m = (1 - h_g) * rho_l + h_g * rho_g
    mu_m = (1 - h_g) * mu_l + h_g * mu_g
    return rho_m, mu_m, h_g


def homogeneo(holdup_l_ns, rho_l, rho_g, mu_l, mu_g):
    rho_m = holdup_l_ns * rho_l + (1 - holdup_l_ns) * rho_g
    mu_m = holdup_l_ns * mu_l + (1 - holdup_l_ns) * mu_g

    return rho_m, mu_m
