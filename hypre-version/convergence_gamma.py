# Charger les données de convergence
convergence_data = pd.read_csv("convergence_gamma.csv")

# Tracé pour différentes valeurs de gamma
plt.figure(figsize=(10, 6))

for gamma_value in convergence_data['gamma'].unique():
    subset = convergence_data[convergence_data['gamma'] == gamma_value]
    plt.plot(subset['iteration'], subset['error'], label=f'γ = {gamma_value}')

plt.xlabel("Nombre d'itérations")
plt.ylabel("Erreur résiduelle")
plt.yscale('log')
plt.title("Influence de γ sur la convergence")
plt.legend()
plt.grid()
plt.savefig("gamma_convergence.png", dpi=300)
plt.show()

