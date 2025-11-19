import numpy as np

"""def razao_solubilidade_gas_agua(Tf, Ppsi):
    A0 = 8.15839
    A1 = -6.12265e-2
    A2 = 1.91663e-4
    A3 = -2.1654e-7

    B0 = 1.01021e-2
    B1 = -7.44241e-5
    B2 = 3.05553e-7
    B3 = -2.94883e-10

    C0 = -9.02505
    C1 = 0.130237
    C2 = -8.53425e-4
    C3 = 2.34122e-6
    C4 = -2.37049e-9

    # Cálculo dos parâmetros A, B e C com T em °F
    A = A0 + A1*Tf + A2*Tf**2 + A3*Tf**3
    B = B0 + B1*Tf + B2*Tf**2 + B3*Tf**3
    C = (C0 + C1*Tf + C2*Tf**2 + C3*Tf**3 + C4*Tf**4) * 1e-7

    # Cálculo de Rsw
    Rsgw = A + B*Ppsi + C*Ppsi**2
    return Rsgw
print("Razão de Solubilidade Gás-Água (Rsgw):", razao_solubilidade_gas_agua(122, 3526), "SCF/STB")"""


def razao_gas_agua(Tf, Ppsi):
    A = 2.12 + (3.45e-3) * Tf - (3.59e-5) * Tf**2
    B = 0.0107 - (5.26e-5) * Tf + (1.48e-7) * Tf**2
    C = 8.75e-7 + (3.9e-9) * Tf - (1.02e-11) * Tf**2
    Rsw = A + B * Ppsi + C * Ppsi**2
    return Rsw


print("Razão de Solubilidade Gás-Água (Rsgw):", razao_gas_agua(122, 3526), "SCF/STB")
