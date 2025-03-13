import pandas as pd
import matplotlib.pyplot as plt

# Lecture du fichier CSV contenant les solutions
data = pd.read_csv('solutions.csv')

plt.figure(figsize=(10,6))
plt.plot(data['x'], data['u_explicit'], marker='o', label='Méthode explicite')
plt.plot(data['x'], data['u_implicit'], marker='s', label='Méthode implicite')
plt.plot(data['x'], data['u_newton'], marker='^', label='Méthode de Newton')
plt.xlabel('x')
plt.ylabel('u')
plt.title('Solution de l\'équation de diffusion non linéaire')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('solution_plot.png', dpi=300)  # Sauvegarde en haute résolution
plt.show()

