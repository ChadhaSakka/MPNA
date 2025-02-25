# Compilateur MPI
CC = mpicc

# Dossier des fichiers d'en-tête
INCLUDE_DIR = include

# Dossier des fichiers source
SRC_DIR = src

# Fichiers source
SRCS = $(SRC_DIR)/sparse_matvec.c $(SRC_DIR)/main.c

# Fichier exécutable
TARGET = sparse_matvec

# Options de compilation
CFLAGS = -I$(INCLUDE_DIR) -Wall
LDFLAGS = -lm

# Fichier de données par défaut
DATA_FILE = data/bcsstk03.mtx

# Règle par défaut : compiler le programme
all: $(TARGET)

# Compilation du programme
$(TARGET): $(SRCS)
	$(CC) $(CFLAGS) $(SRCS) -o $(TARGET) $(LDFLAGS)

# Exécuter en mode sériel (1 processus)
run-serial: $(TARGET)
	mpirun -np 1 ./$(TARGET) $(DATA_FILE)

# Exécuter en mode distribué (4 processus)
run-distributed: $(TARGET)
	mpirun -np 4 ./$(TARGET) $(DATA_FILE)

# Nettoyer les fichiers générés
clean:
	rm -f $(TARGET) *.o

# Phony targets (pas de fichiers associés)
.PHONY: all run-serial run-distributed clean
