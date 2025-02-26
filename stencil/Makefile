# Compiler
CC = gcc

# Compiler flags
CFLAGS = -Wall -Wextra -O3 -march=native -std=c11 -fopenmp
LDFLAGS = -lm  # Link math library

# Executable names
CG_EXEC = cg_solver
GMRES_EXEC = gmres_solver

# Source files
CG_SRC = cg.c
GMRES_SRC = gmres.c

# Object files
CG_OBJ = $(CG_SRC:.c=.o)
GMRES_OBJ = $(GMRES_SRC:.c=.o)

# Default target (compile all)
all: $(CG_EXEC) $(GMRES_EXEC)

# Compile Conjugate Gradient Solver
$(CG_EXEC): $(CG_OBJ)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

# Compile GMRES Solver
$(GMRES_EXEC): $(GMRES_OBJ)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

# Build object files
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Debug build
debug: CFLAGS += -g
debug: clean all

# Clean the compiled files
clean:
	rm -f $(CG_EXEC) $(GMRES_EXEC) $(CG_OBJ) $(GMRES_OBJ)

# Run the CG solver
run_cg: $(CG_EXEC)
	./$(CG_EXEC)

# Run the GMRES solver
run_gmres: $(GMRES_EXEC)
	./$(GMRES_EXEC)

# Phony targets (not real files)
.PHONY: all clean debug run_cg run_gmres

