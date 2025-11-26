# main.py

import numpy as np
import pandas as pd
from modelos import driftflux, homogeneo

from correlacoes import Calculo_agua as ca
from correlacoes import Calculo_Oleo as co
from correlacoes import Calculo_Gas as cg

def formulacao(rho, v_m, mu_m, ap, d_h, epsilon, theta, vazao_massica_m, titulo, temperatura, M_g, P):
    R = 8.314
    reynolds = rho * v_m * d_h / mu_m
    fator_atrito = 0.0055 * (1 + ((2e4*epsilon /d_h) + (1e6/reynolds))**(1/3))
    dp_dl_atrito = -fator_atrito * rho * v_m**2 / (2*d_h)
    dp_dl_gravidade = -rho * 9.71 * np.sin(theta)
    E = -vazao_massica_m**2 * titulo * R * temperatura/(ap**2 * M_g * P**2) 

    dp_dl_total = (dp_dl_atrito + dp_dl_gravidade) / (1 + E)
    dp_dl_aceleracao = E * dp_dl_total
    return dp_dl_total


def perda_de_carga(v_lsc, bsw, rho_o, rho_w, rho_g, mu_g, mu_o, mu_w, rgl, bo, bw, bg, rs, rsw, ap, d_h, epsilon, theta_1, theta_2, M_g, P, comprimento, elementos):
    rho_l = bsw * rho_w + (1-bsw) * rho_o
    mu_l = bsw * mu_w + (1-bsw) * mu_o
    

    def vazoes():
        v_wsc = v_lsc * bsw
        v_osc = v_lsc (1- bsw)
        v_gsc = rgl * v_lsc

        vazao_l = (v_osc * bo) + (v_wsc * bw)
        vazao_g = bg* (v_gsc - (v_osc*rs) - (v_wsc*rsw))

        vazao_massica_l = rho_l * vazao_l
        vazao_massica_g = rho_l * vazao_g
        v_m = vazao_l + vazao_g
        vazao_massica_m = vazao_massica_g + vazao_massica_l

        return vazao_l, vazao_g, vazao_massica_l, vazao_massica_g, v_m, vazao_massica_m

    def velocidades():
        v_l, v_g = vazoes
        v_sl = v_l / ap
        v_sg = v_g / ap

        holdup_l_ns = v_sl / (v_l + v_g)

        return v_sl, v_sg, holdup_l_ns


    vazao_l, vazao_g, vazao_massica_l, vazao_massica_g, v_m, vazao_massica_m = vazoes()
    titulo = vazao_massica_g / vazao_massica_m
    v_sl, v_sg, holdup_l_ns = velocidades()

    def primeiro_trecho():  # Função para calcular tanto a variação da temperatura e da perda de carga para o trecho reservatório - cabeça de poço
        dl = comprimento / elementos

        def calcular_temperatura():
            temperatura = np.zeros(dl)
            for i in range(len(temperatura)):
                temperatura[i] = 6 + 0.4933 * (comprimento - 1650)
                comprimento += dl
            return temperatura

        temperatura = calcular_temperatura()

        def calcular_perda_de_carga(pressao):
            rho_m, mu_m = homogeneo(holdup_l_ns, rho_l, rho_g, mu_l, mu_g)
            temperatura = calcular_temperatura()
            dp_dl_total = formulacao(rho_m, v_m, mu_m, ap, d_h, epsilon, theta_1 , vazao_massica_m, titulo, temperatura, M_g, pressao)
            return dp_dl_total
        
        L=0
        PRESSAO = []
        for i in range(len(temperatura)):
            dp_dl = calcular_perda_de_carga(P)    
            p = p + dp_dl*L
            PRESSAO.append(p)
            L += dl
        return PRESSAO
    pressaaaaao = primeiro_trecho()
    return rho_l


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

    rho_g = cg.massa_especifica_gas(Ppsi, dg, Tr, Mg_g, R_lb)
    mu_g = cg.visco_gas_Lee(Ppsi, Mg_g, dg, R_lb, Tr)

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
    

    teste = perda_de_carga(v_lsc, bsw, rho_o, rho_w, rho_g, mu_g, mu_o, mu_w, rgl, bo, bw, bg, rs, rsw, ap, d_h, epsilon, theta_1, theta_2, M_g, P, comprimento, elementos)
