import numpy as np
import matplotlib.pyplot as plt
import correlacoes.Calculo_agua as ca
import correlacoes.Calculo_Oleo as co
import correlacoes.Calculo_Gas as cg
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
    g = 9.81
    reynolds = rho * v_m * d_h / mu_m
    fator_atrito = 0.0055 * (1 + ((2e4 * epsilon / d_h) + (1e6 / reynolds)) ** (1 / 3))
    dp_dl_atrito = -fator_atrito * rho * v_m**2 / (2 * d_h)
    dp_dl_gravidade = -rho * g * np.sin(theta)
    E = -(vazao_massica_m**2) * titulo * R * temperatura / (ap**2 * M_g * P**2)

    dp_dl_total = (dp_dl_atrito + dp_dl_gravidade) / (1 + E)
    dp_dl_aceleracao = E * dp_dl_total
    return dp_dl_total


def calcula_PVT(P, T):
    Ppsi = P * 0.000145037738
    Tf = (T * 1.8) + 32

    global RGO, dg, M_air_g, api, Tsc_r, Psc_psi, R_lb, S

    Mg_g, Tr, Tk = cg.gas_conversoes_precalc(M_air_g, Tf, dg)

    P_pc, T_pc = cg.psedo_critica(dg)
    P_pr, T_pr = cg.psedo_reduzida(Ppsi, Tr, dg)
    Z = cg.Z_Brill(Ppsi, Tr, dg)
    Bg_m3 = cg.fator_formação_gas(Ppsi, Tr, dg, Tsc_r, Psc_psi)
    Bg_bbl = Bg_m3 / 5.615

    rho_g = cg.massa_especifica_gas(Ppsi, dg, Tr, Mg_g, R_lb)
    mu_g = cg.visco_gas_Lee(Ppsi, Mg_g, dg, R_lb, Tr)

    Tr_o, do = co.oleo_conversoes_precalc(Tf, api)
    Pb = co.ponto_bolha_stand(api, Tf, dg, RGO)
    Rs = co.razao_solubilidade_STANDING(dg, api, Tf, Ppsi, RGO)
    Co = co.compressibilidade_oleo(
        Rs, dg, api, Tf, Tr_o, Ppsi, do, Tsc_r, Psc_psi, RGO, Pb
    )
    Bo = co.fator_formação_STANDING(do, dg, Rs, Tf, Ppsi, Co, api, RGO)
    rho_o = co.massa_especifica_oleo(do, Rs, dg, Bo, Ppsi, Co, api, Tf, RGO)
    mu_oleoD = co.visco_oleoD_BEAL_STAN(api, Tr_o)
    mu_oleoS = co.visco_oleoS_BEAL_STAN(mu_oleoD, Rs)
    mu_oleoSubS = co.visco_oleoSubS_BEAL_STAN(mu_oleoD, RGO, Ppsi, Pb)

    rho_w = ca.rho_w(S)
    Rsw = ca.Rsw(Ppsi, Tf)
    Bw = ca.Bww(Ppsi, Tf)
    mu_w = ca.mu_w(P, T)

    rho_w_SI = rho_w * 16.01846
    rho_g_SI = rho_g * 16.01846
    rho_o_SI = rho_o * 16.01846

    mu_oleoD_SI = mu_oleoD * 0.001
    mu_oleoS_SI = mu_oleoS * 0.001
    mu_oleoSubS_SI = mu_oleoSubS * 0.001
    mu_g_SI = mu_g * 0.001

    Rs_SI = Rs * 0.178108
    Rsw_SI = Rsw
    Bw_SI = Bw
    Bo_SI = Bo
    Bg_SI = Bg_m3
    Pb_SI = Pb * 6894.75729
    Co_SI = Co / 6894.75729
    mu_w_SI = mu_w * 0.001

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
        do,
    )


def perda_de_carga(
    v_lsc, bsw, rgl, ap, d_h, epsilon, theta, M_g, pressao, temperatura, dl, L
):
    def vazoes(v_lsc, bo, bw, bg, rs, rsw, rho_l):
        v_wsc = v_lsc * bsw
        v_osc = v_lsc * (1 - bsw)
        v_gsc = rgl * v_lsc

        vazao_l = (v_osc * bo) + (v_wsc * bw)
        vazao_g = bg * (v_gsc - (v_osc * rs) - (v_wsc * rsw))

        vazao_massica_l = rho_l * vazao_l
        vazao_massica_g = rho_l * vazao_g
        v_m = (vazao_l + vazao_g) / ap
        vazao_massica_m = vazao_massica_g + vazao_massica_l

        return vazao_l, vazao_g, vazao_massica_l, vazao_massica_g, v_m, vazao_massica_m

    def velocidades(vazao_l, vazao_g, ap):
        v_sl = vazao_l / ap
        v_sg = vazao_g / ap
        holdup_l_ns = v_sl / (v_sl + v_sg)
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
            Pb_SI,
            Co,
            Z,
            mu_w,
            api,
            do,
        ) = calcula_PVT(pressao, temperatura)

        rho_l = bsw * rho_w + (1 - bsw) * rho_o

        vazao_l, vazao_g, vazao_massica_l, vazao_massica_g, v_m, vazao_massica_m = (
            vazoes(v_lsc, Bo, Bw, Bg_m3, Rs, Rsw, rho_l)
        )
        v_sl, v_sg, holdup_l_ns = velocidades(vazao_l, vazao_g, ap)

        if pressao > Pb_SI:
            mu_l = bsw * mu_w + (1 - bsw) * mu_oleoSubS
            rho_m, mu_m = homogeneo(holdup_l_ns, rho_l, rho_g, mu_l, mu_g)
        else:
            mu_l = bsw * mu_w + (1 - bsw) * mu_oleoS
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
            theta,
            vazao_massica_m,
            titulo,
            temperatura,
            M_g,
            pressao,
        )

        return dp_dl_total, do, vazao_massica_m, holdup_l_ns, Bo, Bg_m3

    dp_dl, do, vazao_massica_m, holdup_l_ns, Bo, Bg_m3 = calcular_perda_de_carga()
    temp = calc_Temp(L, vazao_massica_m, Tc)
    return dp_dl, temp, holdup_l_ns, Bo, Bg_m3


def plot_pressure_contour(
    comprimento_total,
    elementos,
    pressao_array,
    theta_1_deg,
    theta_2_deg,
    theta_3_deg,
    L1,
    L2,
):
    dl = comprimento_total / elementos
    H_reservoir = 1800

    X_pos = np.zeros(elementos)
    H_pos = np.zeros(elementos)

    L_ANM = L1
    L_Manifold = L1 + L2

    X_ANM = L1 * np.cos(np.deg2rad(theta_1_deg))
    X_Manifold = X_ANM + L2 * np.cos(np.deg2rad(theta_2_deg))

    H_ANM = 1650
    H_Manifold = 1500

    for i in range(elementos):
        posicao_L = i * dl

        if posicao_L <= L_ANM:
            theta_deg = theta_1_deg
            H_start = H_reservoir
            X_start = 0
            Delta_L = posicao_L

            H_pos[i] = H_start - Delta_L * np.sin(np.deg2rad(theta_deg))
            X_pos[i] = X_start + Delta_L * np.cos(np.deg2rad(theta_deg))

        elif posicao_L <= L_Manifold:
            theta_deg = theta_2_deg
            H_start = H_ANM
            X_start = X_ANM
            Delta_L = posicao_L - L_ANM

            H_pos[i] = H_start - Delta_L * np.sin(np.deg2rad(theta_deg))
            X_pos[i] = X_start + Delta_L * np.cos(np.deg2rad(theta_deg))

        else:
            theta_deg = theta_3_deg
            H_start = H_Manifold
            X_start = X_Manifold
            Delta_L = posicao_L - L_Manifold

            H_pos[i] = H_start - Delta_L
            X_pos[i] = X_start

    plt.figure(figsize=(10, 8))
    sc = plt.scatter(
        X_pos, H_pos, c=pressao_array / 1e5, cmap="inferno", s=20, marker="o"
    )
    plt.plot(X_pos, H_pos, "k--", alpha=0.5, linewidth=1)

    plt.colorbar(sc, label="Pressão [bar]")
    plt.xlabel("Distância Horizontal (X) [m]")
    plt.ylabel("Profundidade Vertical (H) [m]")
    plt.ylim(H_reservoir + 50, -50)
    plt.title("Distribuição de Pressão: Caminho da Tubulação (X vs H)")
    plt.grid(True, linestyle=":", alpha=0.6)
    plt.show()


if __name__ == "__main__":
    v_lsc = 6000 / 86400
    theta_1 = 15
    theta_2 = 10
    theta_3 = 90
    bsw = 0.2
    api = 23
    p = 300 * 100000
    sigma_og = 0.00841
    sigma_wg = 0.004
    sigma_lg = bsw * sigma_wg + (1 - bsw) * sigma_og
    tec_poco = 2
    tec_marinho = 1
    p_atm = 100000
    g = 9.81
    S = 0

    Tsc_f = 60
    Tsc_r = Tsc_f + 459.67
    Psc_psi = 14.7
    R_lb = 10.73
    M_air_g = 28.966

    api = 23
    do = 141.5 / (api + 131.5)
    dg = 0.72

    Tc = 80
    Tf = (Tc * 1.8) + 32

    BSW = 0.2
    RGL = 200
    RGO = RGL * 5.61458 / (1 - BSW)

    comprimento_primeiro_trecho = 150 / np.sin(np.deg2rad(theta_1))
    comprimento_segundo_trecho = 150 / np.sin(np.deg2rad(theta_2))
    comprimento_terceiro_trecho = 1500

    L1 = comprimento_primeiro_trecho
    L2 = comprimento_segundo_trecho
    compr1 = L1
    compr2 = L1 + L2
    comprimento_total = compr1 + L2 + comprimento_terceiro_trecho

    elementos = 1000
    dl = comprimento_total / elementos
    pressao = np.zeros(elementos)
    temperatura = np.zeros(elementos)
    holdups = np.zeros(elementos)
    Bo_array = np.zeros(elementos)
    Bg_array = np.zeros(elementos)
    pressao[0] = p
    temperatura[0] = Tc
    d_h = 0.0254 * 8
    ap = np.pi * d_h**2 / 4
    epsilon = 0.0075 * 0.0254
    Mg_ar = 0.02896
    M_g = dg * Mg_ar

    for i in range(elementos - 1):
        posicao = i * dl
        if posicao < compr1:
            theta = theta_1
        elif posicao < compr2:
            theta = theta_2
        else:
            theta = theta_3

        dp_dl, t_now, holdup, Bo_local, Bg_local = perda_de_carga(
            v_lsc,
            bsw,
            RGL,
            ap,
            d_h,
            epsilon,
            np.deg2rad(theta),
            M_g,
            pressao[i],
            temperatura[i],
            dl,
            posicao,
        )

        if pressao[i] < 0:
            break

        pressao[i + 1] = pressao[i] + dp_dl * dl
        temperatura[i + 1] = t_now
        holdups[i] = holdup
        Bo_array[i] = Bo_local
        Bg_array[i] = Bg_local

    L_vector = np.linspace(0, comprimento_total, elementos)

    plot_pressure_contour(
        comprimento_total, elementos, pressao, theta_1, theta_2, theta_3, L1, L2
    )

    plt.figure(figsize=(10, 5))
    plt.plot(L_vector, pressao / 1e5)
    plt.xlabel("Comprimento L [m]")
    plt.ylabel("Pressão [bar]")
    plt.title("Perfil de pressão ao longo do duto")
    plt.grid(True)
    plt.show()

    plt.figure(figsize=(10, 5))
    plt.plot(L_vector, temperatura)
    plt.xlabel("Comprimento L [m]")
    plt.ylabel("Temperatura [°C]")
    plt.title("Perfil de temperatura ao longo do duto")
    plt.grid(True)
    plt.show()

    plt.figure(figsize=(10, 5))
    plt.plot(L_vector[:-1], holdups[:-1])
    plt.xlabel("Comprimento L [m]")
    plt.ylabel("Hold-up de líquido [-]")
    plt.title("Perfil de hold-up de líquido ao longo do duto")
    plt.grid(True)
    plt.show()

    plt.figure(figsize=(10, 5))
    plt.plot((pressao[:-1] / 1e5), Bg_array[:-1], label="Bg")
    plt.plot((pressao[:-1] / 1e5), Bo_array[:-1], label="Bo")
    plt.xlabel("Pressão [bar]")
    plt.ylabel("Fatores de formação [m³/sm³]")
    plt.title("Perfil de Bg e Bo ao longo de L")
    plt.legend()
    plt.grid(True)
    plt.show()

    plt.figure(figsize=(10, 5))
    plt.plot((L_vector[:-1]), Bg_array[:-1], label="Bg")
    plt.plot((L_vector[:-1]), Bo_array[:-1], label="Bo")
    plt.xlabel("Comprimento L [m]")
    plt.ylabel("Fatores de formação [m³/sm³]")
    plt.title("Perfil de Bg e Bo ao longo de L")
    plt.legend()
    plt.grid(True)
    plt.show()
