import Calculo_Gas as cg
import Calculo_Oleo as co
import Calculo_agua as ca

if __name__ == "__main__":
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
    print(f"Ppsi = {Ppsi} psia")

    BSW = 0.2

    RGL = 300  # [sm3/sm3]
    RGO = RGL * 5.61458 / (1 - BSW)  # [scf/stb]
    print(f"RGO = {RGO / 5.61458} Sm3/Sm3")

    S = 0

    # Resultados do Gás

    print("=" * 30)
    print("       RESULTADOS DO GÁS       ")
    print("=" * 30)

    Mg_g, Tr, Tk = cg.gas_conversoes_precalc(M_air_g, Tf, dg)

    P_pc, T_pc = cg.psedo_critica(dg)  # psia e °R
    print(f"P_pc = {P_pc} psia")
    print(f"T_pc = {T_pc} °R")

    P_pr, T_pr = cg.psedo_reduzida(Ppsi, Tr, dg)
    print("P_pr = ", P_pr)
    print("T_pr = ", T_pr)

    Z = cg.Z_Brill(Ppsi, Tr, dg)
    print("Z = ", Z)

    Bg_m3 = cg.fator_formação_gas(Ppsi, Tr, dg, Tsc_r, Psc_psi)
    Bg_bbl = Bg_m3 / 5.615  # bbl/SCF
    print(f"Bg = {Bg_m3} m3/Sm3 = CF/SCF")
    print(f"Bg = {Bg_bbl} bbl/SCF")

    print("_" * 30)
    print("Principais")
    print("_" * 30)

    rho_gas = cg.massa_especifica_gas(Ppsi, dg, Tr, Mg_g, R_lb)
    print(f"rho_gas = {rho_gas} lb/ft^3")

    mu_gas = cg.visco_gas_Lee(Ppsi, Mg_g, dg, R_lb, Tr)
    print(f"mu_gas = {mu_gas} cP")

    # Resultados Do Óleo

    print("\n" + "=" * 30)
    print("       RESULTADOS DO ÓLEO       ")
    print("=" * 30)

    Tr, do = co.oleo_conversoes_precalc(Tf, api)
    print(f"°API = {api}")

    Pb = co.ponto_bolha_stand(api, Tf, dg, RGO)
    print(f"Pb = {Pb} psia")

    Rs = co.razao_solubilidade_STANDING(dg, api, Tf, Ppsi, RGO)
    print(f"Rs = {Rs} SCF/STB")

    Co = co.compressibilidade_oleo(
        Rs, dg, api, Tf, Tr, Ppsi, do, Tsc_r, Psc_psi, RGO, Pb
    )
    print(f"Co = {Co} 1/psia")

    Bo = co.fator_formação_STANDING(do, dg, Rs, Tf, Ppsi, Co, api, RGO)
    print(f"Bo = {Bo} bbl/STB")

    print("_" * 30)
    print("Principais")
    print("_" * 30)

    rho_o = co.massa_especifica_oleo(do, Rs, dg, Bo, Ppsi, Co, api, Tf, RGO)
    print(f"rho_oleo = {rho_o} lb/ft^3")

    mu_oleoD = co.visco_oleoD_BEAL_STAN(api, Tr)
    print(f"mu_oleoDeth = {mu_oleoD} cP")

    mu_oleoS = co.visco_oleoS_BEAL_STAN(mu_oleoD, Rs)
    print(f"mu_oleoSatu = {mu_oleoS} cP")

    mu_oleoSubS = co.visco_oleoSubS_BEAL_STAN(mu_oleoD, RGO, Ppsi, Pb)
    print(f"mu_oleoSubSatu = {mu_oleoSubS} cP")

    print(f"Co = {Co} 1/psia")

    # Resultados Do Óleo

    print("\n" + "=" * 30)
    print("       RESULTADOS DO ÁGUA       ")
    print("=" * 30)

    rho_w = ca.rho_w(S)
    print(f"rho_agua = {rho_w} lb/SCF")

    Rsw = ca.Rsw(Ppsi, Tf)
    print(f"Rsw = {Rsw} sm3/sm3")

    Bw = ca.Bww(Ppsi, Tf)
    print(f"Bw = {Bw} bbl/STB")
