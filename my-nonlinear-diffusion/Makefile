##############################
# Makefile for Nonlinear Diffusion
##############################

# Compiler and flags
CC       = gcc
CFLAGS   = -Wall -Wextra -O3 -std=c11 -Iinclude
LDFLAGS  = -lm   # Link math library

# Directories
SRC_DIR  = src
INC_DIR  = include
BUILD_DIR= build
TEST_DIR = test

# Executable name
TARGET   = diffusion_solver

# Source files
SRCS = $(SRC_DIR)/main.c \
       $(SRC_DIR)/diffusion_linear_implicit.c \
       $(SRC_DIR)/diffusion_newton.c \
       $(SRC_DIR)/utils.c

# Object files (place them in build/ directory)
OBJS = $(SRCS:$(SRC_DIR)/%.c=$(BUILD_DIR)/%.o)

# Default target
all: $(TARGET)

# How to build final executable
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

# Rule to build each .o file
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Create the build directory if it doesn't exist
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Testing (example)
test:
	$(CC) $(CFLAGS) -o test_linear_implicit $(TEST_DIR)/test_linear_implicit.c \
	    $(SRC_DIR)/diffusion_linear_implicit.c \
	    $(SRC_DIR)/utils.c \
	    $(LDFLAGS)
	./test_linear_implicit

# Clean up
clean:
	rm -rf $(BUILD_DIR)/*.o $(TARGET) test_linear_implicit

.PHONY: all clean test

