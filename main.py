# main.py

import numpy as np
from DriftFlux import driftflux 

from correlacoes import Calculo_agua as ca
from correlacoes import Calculo_Oleo as co
from correlacoes import Calculo_Gas as cg


def primeiro_trecho(
    comprimento, elementos
):  # Função para calcular tanto a variação da temperatura e da perda de carga para o trecho reservatório - cabeça de poço
    dl = comprimento / elementos

    def calcular_temperatura():
        temperatura = np.zeros(dl)
        for i in range(len(temperatura)):
            temperatura[i] = 6 + 0.4933 * (comprimento - 1650)
            comprimento += dl
        return temperatura

    def calcular_perda_de_carga():
        return


if __name__ == "__main__":
    v_lsc = 6000 / 86400
    theta_1 = 15  # Trecho reservatório
    theta_2 = 10  # Trecho manifold
    bsw = 0.2
    rgl = 300  # Razão Gás-líquido
    api = 23  # Grau API
    p = 500 * 100000  # Pressão no reservatório [Pa]
    sigma_og = 0.00841  # [N/m]
    sigma_wg = 0.004  # [N/m]
    tec_poco = 2  # TEC no poço [W/m.K]
    tec_marinho = 1  # TEC no ambiente marinho [W/m.K]
    S = 0
    comprimento_primeiro_trecho = 150 / np.sin(theta_1)  # metros

    # Condição Standard
    Tsc_f = 60  # °F
    Tsc_r = Tsc_f + 459.67  # °R
    Psc_psi = 14.7  # psia
    # CONSTANTES
    R_lb = 10.73  # psi ⋅ ft^3 ⋅ lbmol^−1⋅ °R^−1
    M_air_g = 28.966  # g/mol

    # PARAMETROS DE ENTRADA
    dg = 0.72
    api = 23

    # Local de observação
    Tc = 80  # °C
    Tf = (Tc * 1.8) + 32  # °F

    Pbar = 500  # bar
    Ppsi = Pbar * 14.5037738  # psia'''

    BSW = 0.2

    RGL = 300  # [sm3/sm3]
    RGO = RGL * 5.61458 / (1 - BSW)  # [scf/stb]

    S = 0
    # Resultados do Gás

    Mg_g, Tr, Tk = cg.gas_conversoes_precalc(M_air_g, Tf, dg)

    P_pc, T_pc = cg.psedo_critica(dg)  # psia e °R
    P_pr, T_pr = cg.psedo_reduzida(Ppsi, Tr, dg)
    Z = cg.Z_Brill(Ppsi, Tr, dg)
    Bg_m3 = cg.fator_formação_gas(Ppsi, Tr, dg, Tsc_r, Psc_psi)
    Bg_bbl = Bg_m3 / 5.615  # bbl/SCF

    rho_gas = cg.massa_especifica_gas(Ppsi, dg, Tr, Mg_g, R_lb)
    mu_gas = cg.visco_gas_Lee(Ppsi, Mg_g, dg, R_lb, Tr)

    # Resultados Do Óleo
    Tr, do = co.oleo_conversoes_precalc(Tf, api)
    Pb = co.ponto_bolha_stand(api, Tf, dg, RGO)
    Rs = co.razao_solubilidade_STANDING(dg, api, Tf, Ppsi, RGO)
    Co = co.compressibilidade_oleo(
        Rs, dg, api, Tf, Tr, Ppsi, do, Tsc_r, Psc_psi, RGO, Pb
    )
    Bo = co.fator_formação_STANDING(do, dg, Rs, Tf, Ppsi, Co, api, RGO)
    rho_o = co.massa_especifica_oleo(do, Rs, dg, Bo, Ppsi, Co, api, Tf, RGO)
    mu_oleoD = co.visco_oleoD_BEAL_STAN(api, Tr)
    mu_oleoS = co.visco_oleoS_BEAL_STAN(mu_oleoD, Rs)

    mu_oleoSubS = co.visco_oleoSubS_BEAL_STAN(mu_oleoD, RGO, Ppsi, Pb)
    # Resultados Da Água
    rho_w = ca.rho_w(S)
    Rsw = ca.Rsw(Ppsi, Tf)

    Bw = ca.Bww(Ppsi, Tf)
