# main.py

import numpy as np
import pandas as pd


from correlacoes import Calculo_agua as ca
from correlacoes import Calculo_Oleo as co
from correlacoes import Calculo_Gas as cg
from modelos import *
from Temperature import *


def formulacao(
    rho,
    v_m,
    mu_m,
    ap,
    d_h,
    epsilon,
    theta,
    vazao_massica_m,
    titulo,
    temperatura,
    M_g,
    P,
):
    R = 8.314
    reynolds = rho * v_m * d_h / mu_m
    fator_atrito = 0.0055 * (1 + ((2e4 * epsilon / d_h) + (1e6 / reynolds)) ** (1 / 3))
    dp_dl_atrito = -fator_atrito * rho * v_m**2 / (2 * d_h)
    dp_dl_gravidade = -rho * 9.71 * np.sin(theta)
    E = -(vazao_massica_m**2) * titulo * R * temperatura / (ap**2 * M_g * P**2)

    dp_dl_total = (dp_dl_atrito + dp_dl_gravidade) / (1 + E)
    dp_dl_aceleracao = E * dp_dl_total
    return dp_dl_total


def calcula_PVT(P, T):
    S = 0
    api = 23
    Mg_g, Tr, Tk = cg.gas_conversoes_precalc(M_air_g, Tf, dg)

    P_pc, T_pc = cg.psedo_critica(dg)  # psia e °R
    P_pr, T_pr = cg.psedo_reduzida(Ppsi, Tr, dg)
    Z = cg.Z_Brill(Ppsi, Tr, dg)
    Bg_m3 = cg.fator_formação_gas(Ppsi, Tr, dg, Tsc_r, Psc_psi) #m3/sm3
    Bg_bbl = Bg_m3 / 5.615  # bbl/SCF

    rho_g = cg.massa_especifica_gas(Ppsi, dg, Tr, Mg_g, R_lb) #lb/ft^3
    mu_g = cg.visco_gas_Lee(Ppsi, Mg_g, dg, R_lb, Tr) #cP

    # Resultados Do Óleo
    Tr, do = co.oleo_conversoes_precalc(Tf, api) 
    Pb = co.ponto_bolha_stand(api, Tf, dg, RGO) # psia
    Rs = co.razao_solubilidade_STANDING(dg, api, Tf, Ppsi, RGO) # SCF/STB
    Co = co.compressibilidade_oleo(
        Rs, dg, api, Tf, Tr, Ppsi, do, Tsc_r, Psc_psi, RGO, Pb
    ) #  1/psia
    Bo = co.fator_formação_STANDING(do, dg, Rs, Tf, Ppsi, Co, api, RGO) # bbl/STB
    rho_o = co.massa_especifica_oleo(do, Rs, dg, Bo, Ppsi, Co, api, Tf, RGO) #lb/ft^3"
    mu_oleoD = co.visco_oleoD_BEAL_STAN(api, Tr) #cP
    mu_oleoS = co.visco_oleoS_BEAL_STAN(mu_oleoD, Rs) # cP

    mu_oleoSubS = co.visco_oleoSubS_BEAL_STAN(mu_oleoD, RGO, Ppsi, Pb) #cP
    mu_w = ca.mu_w(P, T) #cP
    
    # Resultados Da Água
    rho_w = ca.rho_w(S) # lb/SCF
    Rsw = ca.Rsw(Ppsi, Tf) # sm3/sm3
    Bw = ca.Bww(Ppsi, Tf) # bbl/STB
    mu_w = ca.mu_w(P, T) #cP


    #Conversões
    rho_w_SI = rho_w * 16.01846 #kg/m3
    rho_g_SI = rho_g * 16.01846 #kg/m3
    rho_o_SI = rho_o * 16.01846 #kg/m3

    mu_oleoD_SI = mu_oleoD * 0.001 #Pa.s
    mu_oleoS_SI = mu_oleoS * 0.001 #Pa.s
    mu_oleoSubS_SI = mu_oleoSubS * 0.001 #pa.s
    mu_g_SI = mu_g * 0.001 #Pa.s
    
    Rs_SI = Rs * 0.178108 #sm3/sme
    Rsw_SI = Rsw #sm3/sm3
    Bw_SI = Bw * 0.158987 #m3/sm3
    Bo_SI = Bo * 0.158987 #m3/sm3
    Bg_SI = Bg_m3 #m3/sm3
    Pb_SI = Pb * 6894.75729 #Pa
    Co_SI = Co/6894.75729 #1/Pa
    mu_w_SI = mu_w * 0.001 # Pa.s

    return (
    rho_w_SI,
    rho_g_SI,
    rho_o_SI,
    mu_oleoS_SI,
    mu_oleoD_SI,
    mu_oleoSubS_SI,
    mu_g_SI,
    Rsw_SI,
    Bw_SI,
    Bo_SI,
    Bg_SI,
    Rs_SI,
    Pb_SI,
    Co_SI,
    Z,
    mu_w_SI,
    api,
    do
    )


def perda_de_carga(
    v_lsc, bsw, rgl, ap, d_h, epsilon, theta, M_g, pressao, temperatura, dl, L
):
    def vazoes(v_lsc, bo, bw, bg, rs, rsw, rho_l):
        v_wsc = v_lsc * bsw
        v_osc = v_lsc*(1 - bsw)
        v_gsc = rgl * v_lsc

        vazao_l = (v_osc * bo) + (v_wsc * bw)
        vazao_g = bg * (v_gsc - (v_osc * rs) - (v_wsc * rsw))

        vazao_massica_l = rho_l * vazao_l
        vazao_massica_g = rho_l * vazao_g
        v_m = vazao_l + vazao_g
        vazao_massica_m = vazao_massica_g + vazao_massica_l

        return vazao_l, vazao_g, vazao_massica_l, vazao_massica_g, v_m, vazao_massica_m

    def velocidades(vazao_l, vazao_g, ap):
        #vazao_l, vazao_g = vazoes
        v_sl = vazao_l / ap
        v_sg = vazao_g / ap

        holdup_l_ns = v_sl / (vazao_l + vazao_g)

        return v_sl, v_sg, holdup_l_ns

    

    def calcular_perda_de_carga():
        (
            rho_w,
            rho_g,
            rho_o,
            mu_oleoS,
            mu_oleoD,
            mu_oleoSubS,
            mu_g,
            Rsw,
            Bw,
            Bo,
            Bg_m3,
            Rs,
            Pb,
            Co,
            Z,
            mu_w,
            api,
            do
        ) = calcula_PVT(pressao, temperatura)
        rho_l = bsw * rho_w + (1 - bsw) * rho_o

        vazao_l, vazao_g, vazao_massica_l, vazao_massica_g, v_m, vazao_massica_m = (
            vazoes(v_lsc, Bo, Bw, Bg_m3, Rs, Rsw, rho_l)
        )
        v_sl, v_sg, holdup_l_ns = velocidades(vazao_l, vazao_g, ap)

        if pressao > Pb:
            mu_l = bsw * mu_w(1 - bsw) * mu_oleoSubS
            rho_m, mu_m = homogeneo(holdup_l_ns, rho_l, rho_g, mu_l, mu_g)
        else:
            mu_l = bsw * mu_w*(1 - bsw) * mu_oleoS
            rho_m, mu_m = driftflux(
                v_sg,
                v_m,
                v_sl,
                rho_g,
                rho_l,
                p_atm,
                pressao,
                g,
                d_h,
                sigma_lg,
                theta,
                mu_l,
                mu_g,
            )

        vazao_l, vazao_g, vazao_massica_l, vazao_massica_g, v_m, vazao_massica_m = (
            vazoes(v_lsc, Bo, Bw, Bg_m3, Rs, Rsw, rho_l)
        )
        titulo = vazao_massica_g / vazao_massica_m

        dp_dl_total = formulacao(
            rho_m,
            v_m,
            mu_m,
            ap,
            d_h,
            epsilon,
            theta_1,
            vazao_massica_m,
            titulo,
            temperatura,
            M_g,
            pressao,
        )
        return dp_dl_total, do, vazao_massica_m

    dp_dl, do, vazao_massica_m = calcular_perda_de_carga()

    temp = calc_Temp(L, vazao_massica_m, Tc)

    return dp_dl, temp


if __name__ == "__main__":
    v_lsc = 6000 / 86400
    theta_1 = 15  # Trecho reservatório
    theta_2 = 10  # Trecho manifold
    theta_3 = 90
    bsw = 0.2
    rgl = 300  # Razão Gás-líquido
    api = 23  # Grau API
    p = 500 * 100000  # Pressão no reservatório [Pa]
    sigma_og = 0.00841  # [N/m]
    sigma_wg = 0.004  # [N/m]
    sigma_lg = bsw * sigma_wg + (1 - bsw) * sigma_og
    tec_poco = 2  # TEC no poço [W/m.K]
    tec_marinho = 1  # TEC no ambiente marinho [W/m.K]
    p_atm = 100000
    g = 9.81
    S = 0

    # Condição Standard
    Tsc_f = 60  # °F
    Tsc_r = Tsc_f + 459.67  # °R
    Psc_psi = 14.7  # psia
    # CONSTANTES
    R_lb = 10.73  # psi ⋅ ft^3 ⋅ lbmol^−1⋅ °R^−1
    M_air_g = 28.966  # g/mol

    # PARAMETROS DE ENTRADA
    api = 23
    do = 141.5/(api + 131.5)
    dg = 0.72

    # Local de observação
    Tc = 80  # °C
    Tf = (Tc * 1.8) + 32  # °F

    Pbar = 500  # bar
    Ppsi = Pbar * 14.5037738  # psia'''

    BSW = 0.2

    RGL = 300  # [sm3/sm3]
    RGO = RGL * 5.61458 / (1 - BSW)  # [scf/stb]

    comprimento_primeiro_trecho = 150 / np.sin(theta_1)  # metros
    comprimento_segundo_trecho = 150 / np.sin(theta_2)
    comprimento_terceiro_trecho = 1500
    comprimento_total = (
        comprimento_terceiro_trecho
        + comprimento_primeiro_trecho
        + comprimento_segundo_trecho
    )

    elementos = 1000
    dl = comprimento_total / elementos
    pressao = np.zeros(elementos)
    temperatura = np.zeros(elementos)
    pressao[0] = p
    temperatura[0] = Tc
    d_h = 0.0254 * 7
    ap = np.pi * d_h**2 / 4
    epsilon = 0.0075 * 0.0254
    Mg_ar = 0.02896
    M_g = dg * Mg_ar
    for L in range(len(pressao)-1):
        if L < comprimento_primeiro_trecho:
            dp_dl, t_now = perda_de_carga(
                v_lsc,
                bsw,
                rgl,
                ap,
                d_h,
                epsilon,
                theta_1,
                M_g,
                pressao[L],
                temperatura[L],
                dl,
                L
            )
        elif L < comprimento_segundo_trecho:
            dp_dl, t_now = perda_de_carga(
                v_lsc,
                bsw,
                rgl,
                ap,
                d_h,
                epsilon,
                theta_2,
                M_g,
                pressao[L],
                temperatura[L],
                dl,
                L
            )
        else:
            dp_dl, t_now = perda_de_carga(
                v_lsc,
                bsw,
                rgl,
                ap,
                d_h,
                epsilon,
                theta_3,
                M_g,
                pressao[L],
                temperatura[L],
                dl,
                L
            )

        pressao[L+1] = pressao[L] + dp_dl * dl
        temperatura[L+1] = t_now
        
