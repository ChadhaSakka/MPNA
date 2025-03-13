# Makefile pour projet_diffusion_nonlineaire

CC = mpicc
CFLAGS = -I./include -I/usr/include/hypre -O2 -Wall -std=c99
LDFLAGS = -L/usr/lib/x86_64-linux-gnu -lHYPRE -lm

# Liste des fichiers source
SRC = src/main.c src/diffusion.c src/solveur.c src/utils.c
OBJ = $(SRC:.c=.o)
TARGET = diffusion_solver

# Règle par défaut : compilation de l'exécutable principal
all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Cible pour compiler les tests unitaires
test: $(TARGET)
	@echo "Compilation des tests..."
	$(CC) $(CFLAGS) -o test/test_diffusion test/test_diffusion.c src/diffusion.c src/utils.c $(LDFLAGS)
	$(CC) $(CFLAGS) -o test/test_solveur test/test_solveur.c src/solveur.c src/diffusion.c src/utils.c $(LDFLAGS)
	@echo "Tests compilés. Vous pouvez exécuter 'test/test_diffusion' et 'test/test_solveur'."

# Nettoyage des fichiers objets et exécutables
clean:
	rm -f $(OBJ) $(TARGET) test/test_diffusion test/test_solveur

