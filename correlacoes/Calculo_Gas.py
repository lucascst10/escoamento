import numpy as np


def gas_conversoes_precalc(M_air_g, Tf, dg):
    Mg_g = dg * M_air_g  # g/mol
    Tr = Tf + 459.67  # °R
    Tk = (Tf + 459.67) * (5 / 9)  # K
    return Mg_g, Tr, Tk


def psedo_critica(dg):
    if dg < 0.75:
        P_pc = 677 + 15 * dg - 37.5 * dg**2
        T_pc = 168 + 325 * dg - 12.5 * dg**2
    else:
        P_pc = 706 - 51.7 * dg - 11.1 * dg**2
        T_pc = 187 + 330 * dg - 71.5 * dg**2
    return P_pc, T_pc


def psedo_reduzida(Ppsi, Tr, dg):
    P_pc, T_pc = psedo_critica(dg)
    P_pr = Ppsi / P_pc
    T_pr = Tr / T_pc
    return P_pr, T_pr


def Z_Brill(Ppsi, Tr, dg):
    P_pr, T_pr = psedo_reduzida(Ppsi, Tr, dg)
    A = 1.39 * (T_pr - 0.92) ** 0.5 - 0.36 * T_pr - 0.101

    B = (
        (0.62 - 0.23 * T_pr) * P_pr
        + (0.066 / (T_pr - 0.86) - 0.037) * P_pr**2
        + (0.32 / (10 ** (9 * (T_pr - 1)))) * P_pr**6
    )

    C = 0.132 - 0.32 * np.log10(T_pr)

    D = 10 ** (0.3106 - 0.49 * T_pr + 0.1824 * T_pr**2)

    Z = A + (1 - A) / np.exp(B) + C * P_pr**D
    return Z


def fator_formação_gas(Ppsi, Tr, dg, Tsc_r, Psc_psi):
    Z = Z_Brill(Ppsi, Tr, dg)
    return (Psc_psi * Tr * Z) / (Ppsi * Tsc_r)


def massa_especifica_gas(Ppsi, dg, Tr, Mg_g, R_lb):
    Z = Z_Brill(Ppsi, Tr, dg)
    rho_gas = (Ppsi * Mg_g) / (Z * R_lb * Tr)
    return rho_gas


def visco_gas_Lee(Ppsi, Mg_g, dg, R_lb, Tr):
    rho_gas = massa_especifica_gas(Ppsi, dg, Tr, Mg_g, R_lb)
    xv = 3.448 + (986.4 / Tr) + 0.01009 * Mg_g
    yv = 2.4 - 0.2 * xv
    kv = ((9.379 + 0.0160 * Mg_g) * Tr**1.5) / (209.2 + 19.26 * Mg_g + Tr)

    mu_gas = (1 * 10**-4) * kv * np.exp(xv * (rho_gas / 62.4) ** yv)
    return mu_gas


def compressibilidade_gas(P_pc, P_pr, T_pr, Z):
    dZ = (-3.53 / (10 ** (0.9813 * T_pr))) + (
        (2 * 0.274 * P_pr) / (10 ** (0.8157 * T_pr))
    )
    Cg = (1 / (P_pr * P_pc)) - (1 / (Z * P_pc) * dZ)
    return Cg


'''def compressibilidade_gas(Ppsi, Tr, dg, delta_p=5.0):
    """
    Calcula Cg [1/psi] usando derivada numérica e Z_Brill.
    
    Ppsi : pressão [psia]
    Tr   : temperatura [°R]
    dg   : densidade relativa do gás
    delta_p : incremento para derivada [psia]
    """
    # Z no ponto
    Z = Z_Brill(Ppsi, Tr, dg)
    
    # Derivada numérica (diferença central)
    Z_plus = Z_Brill(Ppsi + delta_p, Tr, dg)
    Z_minus = Z_Brill(max(Ppsi - delta_p, 1.0), Tr, dg)
    dZdp = (Z_plus - Z_minus) / (2 * delta_p)
    
    # Compressibilidade isotérmica
    Cg = (1.0 / Ppsi) * (1.0 - (Ppsi / Z) * dZdp)
    
    return Cg'''
