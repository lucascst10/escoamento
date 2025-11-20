import numpy as np
import Calculo_Gas as cg


def oleo_conversoes_precalc(do, Tf):
    Tr = Tf + 460  # °R
    api = (141.5 / do) - 131.5
    return api, Tr

#Começar com Pb=0 antes do while
def ponto_bolha_stand(api, Tf, dg, Ppsi, Pb, Rs, rgl):
    a = (7.916e-4)*api**1.541 - (4.561e-5)*Tf**1.3911
    if Ppsi > Pb:
        return (112.727 * (rgl ** 0.577421)) / ((dg ** 0.8439) * (10 ** a)) - 1391.051
    else:
        return (112.727 * (Rs ** 0.577421)) / ((dg ** 0.8439) * (10 ** a)) - 1391.051

def razao_solubilidade_STANDING(dg, api, Tf, Ppsi, Pb):
    if Ppsi > Pb:
        Rs = dg * (((Pb / 18.2) + 1.4) * 10 ** ((0.0125 * api) - (0.00091 * Tf))) ** (1 / 0.83)

    else:
        Rs = dg * (((Ppsi / 18.2) + 1.4) * 10 ** ((0.0125 * api) - (0.00091 * Tf))) ** (1 / 0.83)
    return Rs  # SCF/STB


def compressibilidade_oleo(Rs, dg, api, Tf, Tr, Ppsi, do, Pb, Tsc_r, Psc_psi):
    
    if Ppsi >= Pb:  # Usará a correlação de Petrosky e Farshad (1993)
        Co = ((1.705e-7) * (Rs**0.69357) * (dg**0.1885) * (api**0.3272) * (Tf**0.6729) * (Ppsi**-0.5906))
    else:
        Bg_m3 = cg.fator_formação_gas(Ppsi, Tr, dg, Tsc_r, Psc_psi)
        Bg_bbl = Bg_m3 / 5.615  # bbl/SCF
        print(f'Bg = {Bg_bbl} bbl/SCF')
        Bo = 0.9759 + 0.00012*((Rs*(dg/do)**0.5) + (1.25*Tf))**1.2
        dp = 1e-1
        P_dp = (Ppsi + dp)
        Rs_var = razao_solubilidade_STANDING(dg,api,Tf,P_dp,Pb)
        dRs_dp = (Rs_var - Rs)/dp
        Bo_var = 0.9759 + 0.00012*((Rs_var*(dg/do)**0.5) + (1.25*Tf))**1.2
        dBo_dp = (Bo_var - Bo)/dp
        if dBo_dp > (Bg_bbl*dRs_dp):
             print('O resultado será inconsistente')
        Co = -((1/Bo)*dBo_dp) + ((Bg_bbl/Bo)*dRs_dp)
        Rs = razao_solubilidade_STANDING(dg, api, Tf, Ppsi, Pb)
        Bo = 0.9759 + 0.00012 * ((Rs * (dg / do) ** 0.5) + (1.25 * Tf)) ** 1.2
        Co = (-Rs / (Bo * (0.83 * Ppsi + 21.75))) * (
            0.00014 * (dg / do) ** 0.5 * (Rs * (dg / do) ** 0.5 + 1.25 * (Tf)) ** 0.12
            - Bg_bbl
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
    return 0.32 + (((1.8e7) / (api**4.53)) * (360 / (Tr - 260)) ** A) # cP

def visco_oleoS_BEAL_STAN(mu_oleoD, Rs):
    a = 10 ** ((-7.4e-4 * Rs) + (2.2e-7 * Rs**2))
    b = (0.68 / (10 ** (8.62e-5 * Rs))) + \
        (0.25 / (10 ** (1.1e-3 * Rs))) + \
        (0.062 / (10 ** (3.74e-3 * Rs)))
    return a * (mu_oleoD ** b) # cP
