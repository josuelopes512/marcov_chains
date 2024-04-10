import numpy as np
import matplotlib.pyplot as plt

# Parâmetros
mu = 0.1  # Retorno esperado
sigma = 0.2  # Volatilidade
S0 = 100  # Preço inicial
n_years = 5  # Número de anos
n_simulations = 1000  # Número de simulações
dt = 1/252  # Intervalo de tempo (considerando dias úteis)

# Função para gerar simulações de caminhos aleatórios
def monte_carlo_simulation(mu, sigma, S0, n_years, n_simulations, dt):
    # Inicialização da matriz para armazenar os caminhos simulados
    simulations = np.zeros((n_simulations, int(n_years/dt)))
    for i in range(n_simulations):
        # Geração de caminhos aleatórios
        noise = np.random.normal(0, 1, int(n_years/dt))
        path = S0 * np.exp(np.cumsum((mu - 0.5 * sigma ** 2) * dt + sigma * np.sqrt(dt) * noise))
        simulations[i, :] = path
    return simulations

# Execução da simulação de Monte Carlo
simulated_paths = monte_carlo_simulation(mu, sigma, S0, n_years, n_simulations, dt)

# Plotagem dos resultados
plt.figure(figsize=(10, 6))
plt.plot(simulated_paths.T, color='gray', alpha=0.1)
plt.plot(np.mean(simulated_paths, axis=0), color='blue', linewidth=2, label='Média das Simulações')
plt.title('Simulação de Monte Carlo para Previsão de Séries Temporais Financeiras')
plt.xlabel('Dias Úteis')
plt.ylabel('Preço')
plt.legend()
plt.grid(True)
plt.show()
