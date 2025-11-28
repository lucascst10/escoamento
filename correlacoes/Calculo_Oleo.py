import numpy as np
import Calculo_Gas as cg


def oleo_conversoes_precalc(Tf, api):
    Tr = Tf + 460  # °R
    do = 141.5 / (api + 131.5)
    return Tr, do


def ponto_bolha_stand(api, Tf, dg, RGO):
    a = 0.00091 * Tf - 0.0125 * api
    return 18.2 * (((RGO / dg) ** 0.83) * 10**a - 1.4)  # 𝑝𝑠𝑖a


def razao_solubilidade_STANDING(dg, api, Tf, Ppsi, RGO):
    if Ppsi > ponto_bolha_stand(api, Tf, dg, RGO):
        Rs = RGO
    else:
        Rs = dg * (((Ppsi / 18.2) + 1.4) * 10 ** ((0.0125 * api) - (0.00091 * Tf))) ** (
            1 / 0.83
        )
    return Rs  # SCF/STB


def compressibilidade_oleo(Rs, dg, api, Tf, Tr, Ppsi, do, Tsc_r, Psc_psi, RGO, Pb):
    if Ppsi >= Pb:
        Bob = 0.9759 + 0.00012 * ((RGO * (dg / do) ** 0.5) + (1.25 * Tf)) ** 1.2
        rho_ob = (62.4 * do + 0.0136 * RGO * dg) / Bob
        Co = 1e-6 * np.exp(
            (rho_ob + 0.004347 * (Ppsi - Pb) - 79.1)
            / (0.0007141 * (Ppsi - Pb) - 12.938)
        )

    else:
        Bg_m3 = cg.fator_formação_gas(Ppsi, Tr, dg, Tsc_r, Psc_psi)
        Bg_bbl = Bg_m3 / 5.615  # bbl/SCF
        #print(f"Bg = {Bg_bbl} bbl/SCF")
        Bo = 0.9759 + 0.00012 * ((Rs * (dg / do) ** 0.5) + (1.25 * Tf)) ** 1.2
        dp = 1e-1
        P_dp = Ppsi + dp
        Rs_var = razao_solubilidade_STANDING(dg, api, Tf, P_dp, Pb)
        dRs_dp = (Rs_var - Rs) / dp
        Bo_var = 0.9759 + 0.00012 * ((Rs_var * (dg / do) ** 0.5) + (1.25 * Tf)) ** 1.2
        dBo_dp = (Bo_var - Bo) / dp
        #if dBo_dp > (Bg_bbl * dRs_dp):
            #print("O resultado será inconsistente")
        Co = -((1 / Bo) * dBo_dp) + ((Bg_bbl / Bo) * dRs_dp)
        Rs = razao_solubilidade_STANDING(dg, api, Tf, Ppsi, Pb)
        Bo = 0.9759 + 0.00012 * ((Rs * (dg / do) ** 0.5) + (1.25 * Tf)) ** 1.2

        Co = (-Rs / (Bo * (0.83 * Ppsi + 21.75))) * (
            0.00014 * (dg / do) ** 0.5 * (Rs * (dg / do) ** 0.5 + 1.25 * (Tf)) ** 0.12
            - Bg_bbl
        )
    return Co  # 1/psia


def fator_formação_STANDING(do, dg, Rs, Tf, Ppsi, Co, api, RGO):
    Pb = ponto_bolha_stand(api, Tf, dg, RGO)
    if Ppsi <= Pb:
        Bo = 0.9759 + 0.00012 * ((Rs * (dg / do) ** 0.5) + (1.25 * Tf)) ** 1.2
    else:
        Bob = 0.9759 + 0.00012 * ((RGO * (dg / do) ** 0.5) + (1.25 * Tf)) ** 1.2
        Bo = Bob * np.exp(-Co * (Ppsi - Pb))
    return Bo  # bbl/STB


def massa_especifica_oleo(do, Rs, dg, Bo, Ppsi, Co, api, Tf, RGO):
    Pb = ponto_bolha_stand(api, Tf, dg, RGO)
    if Ppsi <= Pb:
        rho_o = (62.4 * do + 0.0136 * Rs * dg) / Bo
    else:
        rho_ob = (62.4 * do + 0.0136 * RGO * dg) / Bo
        rho_o = rho_ob * np.exp(Co * (Ppsi - Pb))
    return rho_o  # lb/ft^3


def visco_oleoD_BEAL_STAN(api, Tr):
    A = 10 ** (0.43 + (8.33 / api))
    return 0.32 + (((1.8e7) / (api**4.53)) * (360 / (Tr - 260)) ** A)  # cP


def visco_oleoS_BEAL_STAN(mu_oleoD, Rs):
    a = 10 ** ((-7.4e-4 * Rs) + (2.2e-7 * Rs**2))
    b = (
        (0.68 / (10 ** (8.62e-5 * Rs)))
        + (0.25 / (10 ** (1.1e-3 * Rs)))
        + (0.062 / (10 ** (3.74e-3 * Rs)))
    )
    return a * (mu_oleoD**b)  # cP


def visco_oleoSubS_BEAL_STAN(mu_oleoD, RGO, Ppsi, Pb):
    mu_oS_Pb = visco_oleoS_BEAL_STAN(mu_oleoD, RGO)
    #print(f"mu_oSatuPb = {mu_oS_Pb} cP")
    return mu_oS_Pb + (
        0.001 * (Ppsi - Pb) * (0.024 * mu_oS_Pb**1.6 + 0.038 * mu_oS_Pb**0.56)
    )

