import numpy as np
from DriftFlux import driftflux

v_sg, v_m, v_sl, rho_g, rho_l, p_atm, P, g, d_h, sigma_l, theta = 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1

if __name__ == "__main__":
    v_lsc = 6000 / 86400
    theta_1 = 15 #Trecho tubulação - reservatório
    