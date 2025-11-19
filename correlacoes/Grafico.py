import numpy as np
import matplotlib.pyplot as plt
import Calculo_Oleo as co

# Dados fixos
dg = 0.84
do = 0.86
Tf = 122
Pb = 5000
Pmin = 100  # pressão mínima do gráfico
Pmax = Pb - 1  # até pressão de bolha
n_pontos = 100  # resolução do gráfico

# Pré-cálculos
api, Tr = co.oleo_conversoes_precalc(do, Tf)

# Lista para armazenar os valores de compressibilidade do óleo
Co_lista = []
pressao_lista = np.linspace(Pmin, Pmax, n_pontos)

for P in pressao_lista:
    Rs = co.razao_solubilidade_STANDING(dg, api, Tf, P, Pb)
    Co_valor = co.compressibilidade_oleo(Rs, dg, api, Tf, Tr, P, do, Pb)
    Co_lista.append(Co_valor)

# Plot
plt.figure(figsize=(8, 5))
plt.plot(pressao_lista, Co_lista, label="Compressibilidade do óleo (Co)")
plt.axvline(x=Pb, color="red", linestyle="--", label="Pressão de bolha (Pb)")
plt.xlabel("Pressão (psia)")
plt.ylabel("Co (1/psia)")
plt.title("Compressibilidade do óleo vs Pressão")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
