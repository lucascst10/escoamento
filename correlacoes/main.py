import Calculo_Gas as cg
import Calculo_Oleo as co

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
    do = 0.86

    # Local de observação
    Tc = 80  # °C
    Tf = (Tc * 1.8) + 32  # °F

    Pbar = 500  # bar
    Ppsi = Pbar * 14.5037738  # psia'''

    S = 2

    Pb = 5000  # psia

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
    print(f"Bg = {Bg_m3} m3/Sm3 = cft/SCF")
    print(f"Bg = {Bg_bbl} bbl/SCF")

    print("_" * 30)
    print("Principais")
    print("_" * 30)

    rho_gas = cg.massa_especifica_gas(Ppsi, dg, Tr, Mg_g, R_lb)
    print(f"rho_gas = {rho_gas} lb/ft^3")

    mu_gas = cg.visco_gas_Lee(Ppsi, Mg_g, dg, R_lb, Tr)
    print(f"mu_gas = {mu_gas} cP")

    Cg = cg.compressibilidade_gas(P_pc, P_pr, T_pr, Z)
    print(f"Cg = {Cg} 1/psia")

    """Cg = cg.compressibilidade_gas(Ppsi, Tr, dg, delta_p=1.0)
    print(f'Cg = {Cg} 1/psia')"""

    # Resultados Do Óleo

    print("\n" + "=" * 30)
    print("       RESULTADOS DO ÓLEO       ")
    print("=" * 30)

    api, Tr = co.oleo_conversoes_precalc(do, Tf)
    print(f"°API = {api}")

    Rs = co.razao_solubilidade_STANDING(dg, api, Tf, Ppsi, Pb)
    print(f"Rs = {Rs} SCF/STB")
    RsPB = co.razao_solubilidade_STANDING(dg, api, Tf, Pb, Pb)
    print(f"RsPB = {RsPB} SCF/STB")

    Co = co.compressibilidade_oleo(Rs, dg, api, Tf, Tr, Ppsi, do, Pb)

    Bo = co.fator_formação_STANDING(do, dg, Rs, Tf, Ppsi, Pb, Co)
    print(f"Bo = {Bo} bbl/STB")

    print("_" * 30)
    print("Principais")
    print("_" * 30)

    rho_oleo = co.massa_especifica_oleo(do, Rs, dg, Bo, Ppsi, Pb, Co)
    print(f"rho_oleo = {rho_oleo} lb/ft^3")

    mu_oleoDeth = co.visco_oleoD_BEAL_STAN(api, Tr)
    print(f"mu_oleoDeth = {mu_oleoDeth} cP")

    mu_oleoSb = co.visco_oleoSaturado_BERGMAN(mu_oleoDeth, Rs)
    print(f"mu_oleoSatu = {mu_oleoSb} cP")

    print(f"Co = {Co} 1/psia")

    # Para comparar as correlações de Co
    if Ppsi >= Pb:
        RsPB = co.RAZAO_PB(dg, api, Tf, Pb)
        BoPB = co.BO_PB(do, dg, RsPB, Tf)
        rho_Pb = co.RHO_PB(do, dg, RsPB, BoPB)
        CoPB = co.Co_PB(rho_Pb, Pb, Ppsi)
        print(f"CoPB = {CoPB} 1/psia")
