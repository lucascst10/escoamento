import numpy as np
import Calculo_Gas as cg


def oleo_conversoes_precalc(do, Tf):
    Tr = Tf + 460  # °R
    api = (141.5 / do) - 131.5
    return api, Tr


def razao_solubilidade_STANDING(dg, api, Tf, Ppsi, Pb):
    if Ppsi > Pb:
        Rs = dg * (((Pb / 18.2) + 1.4) * 10 ** ((0.0125 * api) - (0.00091 * Tf))) ** (
            1 / 0.83
        )
    else:
        Rs = dg * (((Ppsi / 18.2) + 1.4) * 10 ** ((0.0125 * api) - (0.00091 * Tf))) ** (
            1 / 0.83
        )
    return Rs  # SCF/STB


def compressibilidade_oleo(Rs, dg, api, Tf, Tr, Ppsi, do, Pb):
    if Ppsi >= Pb:  # Usará a correlação de Petrosky e Farshad (1993)
        Co = (
            (1.705e-7)
            * (Rs**0.69357)
            * (dg**0.1885)
            * (api**0.3272)
            * (Tf**0.6729)
            * (Ppsi**-0.5906)
        )
    else:
        """Bg = 0.005035*((Z*Tr)/Ppsi)
         print(f'Bg = {Bg} bbl/SCF')
         Bo = 0.9759 + 0.00012*((Rs*(dg/do)**0.5) + (1.25*Tf))**1.2
         dp = 1e-1
         P_dp = (Ppsi + dp)
         Rs_var = razao_solubilidade_STANDING(dg,api,Tf,P_dp,Pb)
         dRs_dp = (Rs_var - Rs)/dp
         Bo_var = 0.9759 + 0.00012*((Rs_var*(dg/do)**0.5) + (1.25*Tf))**1.2
         dBo_dp = (Bo_var - Bo)/dp
         if dBo_dp > (Bg*dRs_dp):
             print('O resultado será inconsistente')
         Co = -((1/Bo)*dBo_dp) + ((Bg/Bo)*dRs_dp)"""
        Bg = cg.fator_formação_gas(Ppsi, Tr, dg)
        Rs = razao_solubilidade_STANDING(dg, api, Tf, Ppsi, Pb)
        Bo = 0.9759 + 0.00012 * ((Rs * (dg / do) ** 0.5) + (1.25 * Tf)) ** 1.2
        Co = (-Rs / (Bo * (0.83 * Ppsi + 21.75))) * (
            0.00014 * (dg / do) ** 0.5 * (Rs * (dg / do) ** 0.5 + 1.25 * (Tf)) ** 0.12
            - Bg
        )
    return Co  # 1/psia


def fator_formação_STANDING(do, dg, Rs, Tf, Ppsi, Pb, Co):
    if Ppsi <= Pb:
        Bo = 0.9759 + 0.00012 * ((Rs * (dg / do) ** 0.5) + (1.25 * Tf)) ** 1.2
    else:
        Bob = 0.9759 + 0.00012 * ((Rs * (dg / do) ** 0.5) + (1.25 * Tf)) ** 1.2
        Bo = Bob * np.exp(-Co * (Ppsi - Pb))
    return Bo  # bbl/STB


def massa_especifica_oleo(do, Rs, dg, Bo, Ppsi, Pb, Co):
    if Ppsi <= Pb:
        rho_oleo = (62.4 * do + 0.0136 * Rs * dg) / Bo
    else:
        rho_oleo_b = (62.4 * do + 0.0136 * Rs * dg) / Bo
        rho_oleo = rho_oleo_b * np.exp(Co * (Ppsi - Pb))
    return rho_oleo  # lb/ft^3


def visco_oleoD_BEAL_STAN(api, Tr):
    A = 10 ** (0.43 + (8.33 / api))
    mu_oleoDeth = 0.32 + (((1.8e7) / (api**4.53)) * (360 / (Tr - 260)) ** A)
    return mu_oleoDeth  # cP


def visco_oleoSaturado_BERGMAN(mu_oleoD_Berg, Rs):
    a = np.exp(4.768 - (0.8359 * np.log(Rs + 300)))
    b = 0.555 + (133.5 / (Rs + 300))
    mu_oleoSatu = a * (mu_oleoD_Berg**b)
    return mu_oleoSatu  # cP


# Para comparar as correlações de Co
def RAZAO_PB(dg, api, Tf, Pb):
    RsPB = dg * (((Pb / 18.2) + 1.4) * 10 ** ((0.0125 * api) - (0.00091 * Tf))) ** (
        1 / 0.83
    )
    return RsPB  # SCF/STB


def BO_PB(do, dg, RsPB, Tf):
    BoPB = 0.9759 + 0.00012 * ((RsPB * (dg / do) ** 0.5) + (1.25 * Tf)) ** 1.2
    return BoPB


def RHO_PB(do, dg, RsPB, BoPB):
    rho_Pb = (62.4 * do + 0.0136 * RsPB * dg) / BoPB
    return rho_Pb


def Co_PB(rho_Pb, Pb, Ppsi):
    CoPB = (10e-6) * np.exp(
        (rho_Pb + (0.004347 * (Ppsi - Pb)) - 79.1)
        / ((0.0007141 * (Ppsi - Pb)) - 12.938)
    )
    return CoPB
