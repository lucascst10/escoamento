import numpy as np
from DriftFlux import driftflux

def primeiro_trecho(comprimento, elementos): #Função para calcular tanto a variação da temperatura e da perda de carga para o trecho reservatório - cabeça de poço
    dl = comprimento / elementos
    def calcular_temperatura():  
        temperatura = np.zeros(dl)

if __name__ == "__main__":
    v_lsc = 6000 / 86400
    theta_1 = 15 #Trecho reservatório
    theta_2 = 10 #Trecho manifold
    bsw = 0.2 
    rgl = 300 #Razão Gás-líquido
    api = 23 #Grau API
    p = 500 * 100000 #Pressão no reservatório [Pa]
    sigma_og = 0.00841 #[N/m]
    sigma_wg = 0.004 #[N/m]
    tec_poco = 2 #TEC no poço [W/m.K]
    tec_marinho = 1 #TEC no ambiente marinho [W/m.K]
    comprimento_primeiro_trecho = 150 / np.sin(theta_1) #metros
