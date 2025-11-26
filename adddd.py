from PVT_código import *
from Modelo_Homogênio import *
from Modelo_Drift_Flux import *
import math
import numpy as np
import matplotlib.pyplot as plt

# Dados geométricos do problema.
#######################################
Vazao_total = 4500  # [sm3/d], objetivo a ser alcançado.
H_final = 1100  # [m]
H_bomba = 750  # [m]
H_manif = 850  # [m]
theta = 37
d_max = 2.5 * 0.0254  # [pol] pra m O GANHO DE PRESSÃO DA BOMBA NÃO FICOU ITERATIVO POR PROBLEMAS DE CONVERGÊNCIA(Linha 283)
########################################
# Caracterização do fluído quando.
########################################
BSW = 0.2  # “Basic Sediment and Water”.
RGL = 150  # [sm3/sm3], Razão gás-líquido.
RGO = RGL / (1 - BSW)
API = 25  # como <30 o modelo Black-Oil pode ser aplicado.
dg = 0.75  # >= 0.75, temos então gás úmido.
sigma_og = 0.00841  # [N/m], tensão interfacial de óleo e gás.
sigma_wg = 0.004  # [N/m], tensão interfacial de água e gás.
TEC_poço = 2  # [W/mK]
TEC_mar = 1  # [W/mK]
S = 0  # % salinidade
rho_ar = 0.076474  # lb/ft3
cp_w = 4180  # J/kg°C
cp_g_met = 2220  # J/kg°C
cp_g_eta = 1750  # J/kg°C
e = 0.152 / 1000  # Ferro galvanizado exemplo
########################################
# Dados standard obtidos apartir das Vazões volumétricas
V_L_sc = (4500 / (60 * 60 * 24))  # m3/s
V_o_sc = V_L_sc * (1 - BSW)
V_w_sc = V_L_sc * BSW
V_g_sc = RGL * V_L_sc


def V_l(RGO, do, P, T_F, API, dg, S):
    Rs = Standing_Rs(RGO, P, dg, T_F, API)
    Bo = Bo_func(RGO, Rs, do, P, T_F, API, dg, S)
    Bww = Bw(P, T_F)
    V_L = V_o_sc * Bo + V_w_sc * Bww

    return V_L


def V_g(P, T_F, V_g_sc, V_o_sc, V_w_sc):
    T_R = T_F + 459.67
    P_pr, T_pr = P_pr_T_pr(dg, T_R, P)
    Z = Brill_Begs_Z(T_pr, P_pr)
    Bgg = Bg(Z, T_F, P)
    Rsww = Rsw(P, T_F)
    Rs = Standing_Rs(RGO, P, dg, T_F, API)
    termo1 = V_o_sc * Rs
    termo2 = V_w_sc * Rsww
    V_G = (V_g_sc - termo1 - termo2) * Bgg
    return V_G


#########################################
# Hold Up no slip

def sigmal(Vl, VG): # Holdup no slip
    HL = Vl / (VG + Vl)
    return HL


# sigma_l = sigmal(Vl, VG)
#########################################
# Teste perfil de temperatura
do_15_15 = do_API(API)
T_teste = np.arange(80, 17, -0.1)  


def do_T(do_15_15, T_teste):
    do = do_15_15 - (5.93 * 10 ** (-4)) * (T_teste - 15)  # Temperatura em °C
    return do


#########################################
# Variação do do em função da Temperatura em °C
do_fT = do_T(do_15_15, T_teste)
rho_o = rho_w(S) * do_fT  # Densidade do óleo em lb/ft3
rho_g = rho_ar * dg  # Densidade do gás em lb/ft3

# Chute inicial para a interpolação
R = 10.7316  # psi·ft^3 / lb-mol·R (constante dos gases no sistema inglês)
T = 519.67  # Temperatura em Rankine (60°F)
P = 14.7  # Pressão em psi (1atm)
Z = 0.8  # Fator de compressibilidade (teste PVT)
MM_met = 16.043  # lb-mol (métano)
MM_eta = 30.070  # lb-mol (etano)

# Cálculo da massa molar do gás
MM_g = (rho_g * Z * R * T) / P

# Interpolação para fração de metano (X) e etano (Y = 1 - X)
X = (MM_g - MM_eta) / (MM_met - MM_eta)
Y = 1 - X

HL = 1


# Cálculo das capacidades térmicas
def cp_o(do_15_15, T_teste):
    do = do_T(do_15_15, T_teste)
    cp_o = (2 * 10 ** (-3) * T_teste - 1.429) * do + (2.67 * 10 ** (-3)) * T_teste + 3.049  # kJ/kg°C
    return cp_o * 1000 # J/kg°C


# Capacidade térmica da mistura
cp_o_values = [cp_o(do_15_15, T) for T in T_teste]
cp_g = X * cp_g_met + Y * cp_g_eta


def cp_m(cp_w, cp_g, cp_o_values, HL, BSW, do_15_15, T_teste):
    cp_mlist = []
    for i in cp_o_values:
        cp_m = HL * ((1 - BSW) * i + (BSW * cp_w)) + (1 - HL) * cp_g
        cp_mlist.append(cp_m)
    return cp_mlist


def vazaoMass(Vl, Vg, rho_w, rho_o, rho_g):
    rhow = rho_w * 16.0185 # lb/ft³ pra kg/m³
    rhoo = rho_o * 16.0185
    rhog = rho_g * 16.0185
    rhol = rhow + rhoo
    massico = Vl * rhol + Vg * rhog

    return massico


cp_m_value = cp_m(cp_w, cp_g, cp_o_values, 1, BSW, do_15_15, T_teste)  # J/kg°C
# chute inicial para o primeiro ponto da tubulação
###################################################################################
# Condições de contorno e temperaturas
T_plat = 17  # [°C]
T_reser = 80  # [°C]
P_plat = 14.5038  # [psi]
P_reser = 400 * 14.5038  # [psi]
T_inf_poço = 80  # [°C]
T_inf_mar = 4  # [°C]

# Ângulos e dimensões geométricas
theta_1 = 90  # Graus ascendente
theta_2 = 323  # Graus descendente
L1 = 349  # [m], Vertical
L3 = 849  # [m], Vertical
L2 = 165

# Teste código main
# Cálculo da variação de temperatura ao longo dos trechos
i = 0
i_new = 0
i_old = 0
Tinf1_old = [T_inf_poço]

HL_list = []
P_perca = []
P_perca.append(P_reser)
T = [T_reser]
do_15_15 = 0.8

# PVT vetores do óleo
do = []
Pb = []
Rs = []
Co = []
Bo = []
rho_o = []
mi_o = []
# PVT vetores do gás
P_pr = []
T_pr = []
Z = []
vm = []

# Frações de perda de carga
delta_p_acc1 = []
delta_p_f = []
delta_p_g = []
delta_p_t = []
Trechos = [L1, L2, L3]


a = (T_inf_mar - T_inf_poço) / L3
b = (T_plat - T_inf_mar)/L2
step = 1

while i != round(Trechos[0] + 1):
    P_atual = P_perca[-1]
    T_atual = T[-1]
    T_atual_F = [T[-1] * 9 / 5 + 32]
    # PVT
    do1 = do_T(do_15_15, T_atual_F[0])
    do.append(do1)
    Rs1 = Standing_Rs(RGO, P_atual, dg, T_atual_F[0], API)
    Rs.append(Rs1)
    Pb1 = Standing_Pb(RGO, dg, T_atual_F[0], API)
    Pb.append(Pb1)
    Co1 = Co_func(P_atual, RGO, dg, do1, T_atual_F[0], API, S, y_CO2=0, y_H2S=0, y_N2=0)
    Co.append(Co1)
    Bo1 = Bo_func(RGO, Rs1, do1, P_atual, T_atual_F[0], API, dg, S)
    Bo.append(Bo1)
    rho_o1 = rho_ofunc(RGO, do1, dg, T_atual_F[0], API, P_atual, S)
    rho_o.append(rho_o1)
    mio = mi_o_func(RGO, dg, API, T_atual_F[0], P_atual)
    mi_o.append(mio)

    # Fase Gás
    Ppr1, Tpr1 = P_pr_T_pr(dg, T_atual_F[0], P_atual, y_CO2=0, y_H2S=0, y_N2=0)
    P_pr.append(Ppr1)
    T_pr.append(Tpr1)
    Z1 = Brill_Begs_Z(Tpr1, Ppr1)
    Z.append(Z1)
    Bg1 = Bg(Z1, T_atual_F[0], P_atual)
    rho_g1 = rho_g * 16.0185
    mig = lee_et_al_mu_g(T_atual_F[0], MM_g, rho_g1)

    # Fase água
    rhow = rho_w(0)
    Rsw1 = Rsw(P_atual, T_atual_F[0])
    Bw1 = Bw(P_atual, T_atual_F[0])
    miw = mi_w(0, P_atual, T_atual_F[0])

    # Calcular as vazões
    Vl = V_o_sc * Bo1 + V_w_sc * Bw1
    Vg = (V_g_sc - V_o_sc * Rs1 - V_w_sc * Rsw1) * Bg1
    if Vg<0:
        Vg = 0
    Hl = sigmal(Vl, Vg)
    HL_list.append(Hl)
    mil = mio + miw

    Cpo = [cp_o(do_15_15, T_atual_F[0])]
    Cpm = cp_m(cp_w, cp_g, Cpo, HL, BSW, do_15_15, T_atual_F[0])  # q aconteceu aq
    m_m = (vazaoMass(Vl, Vg, rhow, rho_o1, rho_g1))
    rhol = rho_o1 + rhow
    ml = Vl * (rho_o1 + rhow)
    mg = Vg * rho_g1
    rho_m = Hl * rhol + (1 - Hl) * rho_g1
    Ap = math.pi * (d_max / 2)**2
    vm1 = (ml + mg) / (rho_m * Ap)
    vm.append(vm1)

    dl = i_new - i_old
    Tinf1 = T_inf_poço + (a) * i  # Variação da Temperatura externa em função do dL
    Tf_1 = m_m * 9.81 * np.sin(math.radians(theta_1)) / TEC_poço
    Tf_2 = np.exp((((-TEC_poço) / (m_m * Cpm[0])) * dl))  # arrumar essa bomba
    Tf_3 = ((Tinf1) - (Tf_1) - (T_atual))
    Tf = Tinf1 - Tf_1 - (Tf_2 * Tf_3)
    T.append(Tf)
    if Pb1<P:
        dp_dl_f = calcular_gradiente_atrito_H(e, d_max, Hl, mil, mig, rho_m, vm)
        dp_dl_g = calcular_gradiente_gravidade(rho_m, Theta_1)
        dp_dl_a = 0
        dp_dl = dp_dl_f + dp_dl_g
    else:
        dp_dl, delta_p_acc, delta_p_fric, delta_p_grav = Modelo_Homogêneo(rho_m, vm1, dg, d_max, e, theta_1, ml, mg, Hl, rhol, rho_g1, mil, mig, Tf, Z1)
    delta_p_t.append(dp_dl)
    delta_p_acc1.append(delta_p_acc)
    delta_p_f.append(delta_p_fric)
    delta_p_g.append(delta_p_grav)
    dp_dl_psi = dp_dl / 6894.76
    DP = dp_dl_psi * dl
    P_perca.append(P_atual - DP)
    i_old = i
    if i == 350:
        print(i)
    i += 1
    i_new = i
    

i = 0

while i != round(Trechos[1] + 1):
    if i == 0:
        P_atual = P_perca[-1] + 195*14.5 # 195 bar
    else:
        P_atual = P_perca[-1]
    T_atual = T[-1]
    T_atual_F = [T[-1] * 9 / 5 + 32]
    # PVT
    do1 = do_T(do_15_15, T_atual_F[0])
    do.append(do1)
    Rs1 = Standing_Rs(RGO, P_atual, dg, T_atual_F[0], API)
    Rs.append(Rs1)
    Pb1 = Standing_Pb(RGO, dg, T_atual_F[0], API)
    Pb.append(Pb1)
    Co1 = Co_func(P_atual, RGO, dg, do1, T_atual_F[0], API, S, y_CO2=0, y_H2S=0, y_N2=0)
    Co.append(Co1)
    Bo1 = Bo_func(RGO, Rs1, do1, P_atual, T_atual_F[0], API, dg, S)
    Bo.append(Bo1)
    rho_o1 = rho_ofunc(RGO, do1, dg, T_atual_F[0], API, P_atual, S)
    rho_o.append(rho_o1)
    mio = mi_o_func(RGO, dg, API, T_atual_F[0], P_atual)
    mi_o.append(mio)

    # Fase Gás
    Ppr1, Tpr1 = P_pr_T_pr(dg, T_atual_F[0], P_atual, y_CO2=0, y_H2S=0, y_N2=0)
    P_pr.append(Ppr1)
    T_pr.append(Tpr1)
    Z1 = Brill_Begs_Z(Tpr1, Ppr1)
    Z.append(Z1)
    Bg1 = Bg(Z1, T_atual_F[0], P_atual)
    rho_g1 = rho_g * 16.0185
    mig = lee_et_al_mu_g(T_atual_F[0], MM_g, rho_g1)

    # Fase água
    rhow = rho_w(0)
    Rsw1 = Rsw(P_atual, T_atual_F[0])
    Bw1 = Bw(P_atual, T_atual_F[0])
    miw = mi_w(0, P_atual, T_atual_F[0])

    # Calcular as vazões
    Vl = V_o_sc * Bo1 + V_w_sc * Bw1
    Vg = (V_g_sc - V_o_sc * Rs1 - V_w_sc * Rsw1) * Bg1
    if Vg<0:
        Vg = 0
    Hl = sigmal(Vl, Vg)
    HL_list.append(Hl)
    mil = mio + miw

    Cpo = [cp_o(do_15_15, T_atual_F[0])]
    Cpm = cp_m(cp_w, cp_g, Cpo, HL, BSW, do_15_15, T_atual_F[0])  # 
    m_m = (vazaoMass(Vl, Vg, rhow, rho_o1, rho_g1))
    rhol = rho_o1 + rhow
    ml = Vl * (rho_o1 + rhow)
    mg = Vg * rho_g1
    rho_m = Hl * rhol + (1 - Hl) * rho_g1
    Ap = math.pi * (d_max / 2)**2
    vm1 = (ml + mg) / (rho_m * Ap)
    vm.append(vm1)

    dl = i_new - i_old
    Tinf1 = T_inf_poço + (b) * i  # Variação da Temperatura externa em função do dL
    Tf_1 = m_m * 9.81 * np.sin(math.radians(theta_2)) / TEC_mar
    Tf_2 = np.exp((((-TEC_mar) / (m_m * Cpm[0])) * dl))  # arrumar essa bomba
    Tf_3 = ((Tinf1) - (Tf_1) - (T_atual))
    Tf = Tinf1 - Tf_1 - (Tf_2 * Tf_3)
    T.append(Tf)
    if Pb1<P:
        dp_dl, delta_p_acc, delta_p_fric, delta_p_grav = Modelo_Homogêneo(rho_m, vm1, dg, d_max, e, theta_1, ml, mg, Hl, rhol, rho_g1, mil, mig, Tf, Z1)
    else:
        vm2, Hl2 = calcular_velocidade_mistura(vm1, d_max, theta_2, Hl)
        rho_m = Hl * rhol + (1 - Hl) * rho_g1
        dp_dl, delta_p_acc, delta_p_fric, delta_p_grav = Modelo_Homogêneo(rho_m, vm2, dg, d_max, e, theta_2, ml, mg, Hl2, rhol, rho_g1, mil, mig, Tf, Z1)
    delta_p_t.append(dp_dl)
    delta_p_acc1.append(delta_p_acc)
    delta_p_f.append(delta_p_fric)
    delta_p_g.append(delta_p_grav)
    dp_dl_psi = dp_dl / 6894.76
    DP = dp_dl_psi * dl
    P_perca.append(P_atual - DP)
    if i == 166:
        print(i)
    i_old = i
    i += 1
    i_new = i
    
i = 0

while i != round(Trechos[2] + 1):
    P_atual = P_perca[-1]
    T_atual = T[-1]
    T_atual_F = [T[-1] * 9 / 5 + 32]
    # PVT
    do1 = do_T(do_15_15, T_atual_F[0])
    do.append(do1)
    Rs1 = Standing_Rs(RGO, P_atual, dg, T_atual_F[0], API)
    Rs.append(Rs1)
    Pb1 = Standing_Pb(RGO, dg, T_atual_F[0], API)
    Pb.append(Pb1)
    Co1 = Co_func(P_atual, RGO, dg, do1, T_atual_F[0], API, S, y_CO2=0, y_H2S=0, y_N2=0)
    Co.append(Co1)
    Bo1 = Bo_func(RGO, Rs1, do1, P_atual, T_atual_F[0], API, dg, S)
    Bo.append(Bo1)
    rho_o1 = rho_ofunc(RGO, do1, dg, T_atual_F[0], API, P_atual, S)
    rho_o.append(rho_o1)
    mio = mi_o_func(RGO, dg, API, T_atual_F[0], P_atual)
    mi_o.append(mio)

    # Fase Gás
    Ppr1, Tpr1 = P_pr_T_pr(dg, T_atual_F[0], P_atual, y_CO2=0, y_H2S=0, y_N2=0)
    P_pr.append(Ppr1)
    T_pr.append(Tpr1)
    Z1 = Brill_Begs_Z(Tpr1, Ppr1)
    Z.append(Z1)
    Bg1 = Bg(Z1, T_atual_F[0], P_atual)
    rho_g1 = rho_g * 16.0185
    mig = lee_et_al_mu_g(T_atual_F[0], MM_g, rho_g1)

    # Fase água
    rhow = rho_w(0)
    Rsw1 = Rsw(P_atual, T_atual_F[0])
    Bw1 = Bw(P_atual, T_atual_F[0])
    miw = mi_w(0, P_atual, T_atual_F[0])

    # Calcular as vazões
    Vl = V_o_sc * Bo1 + V_w_sc * Bw1
    Vg = (V_g_sc - V_o_sc * Rs1 - V_w_sc * Rsw1) * Bg1
    if Vg<0:
        Vg = 0
    Hl = sigmal(Vl, Vg)
    HL_list.append(Hl)
    mil = mio + miw

    Cpo = [cp_o(do_15_15, T_atual_F[0])]
    Cpm = cp_m(cp_w, cp_g, Cpo, HL, BSW, do_15_15, T_atual_F[0])  
    m_m = (vazaoMass(Vl, Vg, rhow, rho_o1, rho_g1))
    rhol = rho_o1 + rhow
    ml = Vl * (rho_o1 + rhow)
    mg = Vg * rho_g1
    rho_m = Hl * rhol + (1 - Hl) * rho_g1
    Ap = math.pi * (d_max / 2)**2
    vm1 = (ml + mg) / (rho_m * Ap)
    vm.append(vm1)

    dl = i_new - i_old
    Tinf1 = T_inf_poço + (b) * i  # Variação da Temperatura externa em função do dL
    Tf_1 = m_m * 9.81 * np.sin(math.radians(theta_1)) / TEC_mar
    Tf_2 = np.exp((((-TEC_mar) / (m_m * Cpm[0])) * dl))  # arrumar essa bomba
    Tf_3 = ((Tinf1) - (Tf_1) - (T_atual))
    Tf = Tinf1 - Tf_1 - (Tf_2 * Tf_3)
    T.append(Tf)
    if Pb1<P:
        dp_dl, delta_p_acc, delta_p_fric, delta_p_grav = Modelo_Homogêneo(rho_m, vm1, dg, d_max, e, theta_1, ml, mg, Hl, rhol, rho_g1, mil, mig, Tf, Z1)
    else:
        vm2, Hl2 = calcular_velocidade_mistura(vm1, d_max, theta_1, Hl)
        rho_m = Hl * rhol + (1 - Hl) * rho_g1
        dp_dl, delta_p_acc, delta_p_fric, delta_p_grav = Modelo_Homogêneo(rho_m, vm2, dg, d_max, e, theta_1, ml, mg, Hl2, rhol, rho_g1, mil, mig, Tf, Z1)
    delta_p_t.append(dp_dl)
    delta_p_acc1.append(delta_p_acc)
    delta_p_f.append(delta_p_fric)
    delta_p_g.append(delta_p_grav)
    dp_dl_psi = dp_dl / 6894.76
    DP = dp_dl_psi * dl
    P_perca.append(P_atual - DP)
    if i == 850:
        print(i)
    i_old = i
    i += 1
    i_new = i

aceleração = sum(delta_p_acc1)/1366
gavitacional = sum(delta_p_g)/1366
fricção = sum(delta_p_f)/1366
total = sum(delta_p_t)/1366



x = np.linspace(0, L3 + L2 + L1 , len(T))
plt.plot(x, T)  # Garantindo tamanhos compatíveis
plt.xlabel('Comprimento [m]')
plt.ylabel('Temperatura [°C]')
plt.grid(True)
plt.show()

y = np.linspace(0, L3 + L2 + L1, len(P_perca))
plt.plot(y, P_perca)  # Garantindo tamanhos compatíveis
plt.xlabel('Comprimento [m]')
plt.ylabel('Pressão [psi]')
plt.grid(True)
plt.show()

z = np.linspace(0, L3 + L2 + L1, len(HL_list))
plt.plot(z, HL_list)  # Garantindo tamanhos compatíveis
plt.xlabel('Comprimento [m]')
plt.ylabel('Hl')
plt.grid(True)
plt.show()

# Plot PVT

a = np.linspace(0, L3 + L2 + L1 , len(vm))
plt.plot(a, vm)  # Garantindo tamanhos compatíveis
plt.xlabel('Comprimento [m]')
plt.ylabel('Velocidade da mistura [m/s]')
plt.grid(True)
plt.show() 


b = np.linspace(0, L3 + L2 + L1 , len(Bo))
plt.plot(b, Bo)  # Garantindo tamanhos compatíveis
plt.xlabel('Comprimento [m]')
plt.ylabel('FVF do óleo [bbl/STB]')
plt.grid(True)
plt.show() 


c = np.linspace(0, L3 + L2 + L1 , len(Z))
plt.plot(c, Z)  # Garantindo tamanhos compatíveis
plt.xlabel('Comprimento [m]')
plt.ylabel('Fator de compressibilidade [ad]')
plt.grid(True)
plt.show() 

d = np.linspace(0, L3 + L2 + L1 , len(rho_o))
plt.plot(d, rho_o)  # Garantindo tamanhos compatíveis
plt.xlabel('Comprimento [m]')
plt.ylabel('Densidade do óleo [lb/ft³]')
plt.grid(True)
plt.show() 

e = np.linspace(0, L3 + L2 + L1 , len(delta_p_acc1))
plt.plot(e, delta_p_acc1)  # Garantindo tamanhos compatíveis
plt.xlabel('Comprimento [m]')
plt.ylabel('Perda de carga (aceleração) [Pa/m]')
plt.grid(True)
plt.show() 

f = np.linspace(0, L3 + L2 + L1 , len(delta_p_f))
plt.plot(f, delta_p_f)  # Garantindo tamanhos compatíveis
plt.xlabel('Comprimento [m]')
plt.ylabel('Perda de carga (fricção) [Pa/m]')
plt.grid(True)
plt.show() 

g = np.linspace(0, L3 + L2 + L1 , len(delta_p_g))
plt.plot(g, delta_p_g)  # Garantindo tamanhos compatíveis
plt.xlabel('Comprimento [m]')
plt.ylabel('Perda de carga (gravitacional) [Pa/m]')
plt.grid(True)
plt.show() 

h = np.linspace(0, L3 + L2 + L1 , len(delta_p_t))
plt.plot(h, delta_p_t)  # Garantindo tamanhos compatíveis
plt.xlabel('Comprimento [m]')
plt.ylabel('Perda de carga (Total) [Pa/m]')
plt.grid(True)
plt.show() 