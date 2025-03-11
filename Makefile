CC = gcc
MPICC = mpicc
CFLAGS = -Wall #-O2 
SRC_DIR = src
BIN_DIR = bin

# Archivos fuente
SECUENCIAL_SRC = $(SRC_DIR)/FourierS.c
SECUENCIAL_SRC2 = $(SRC_DIR)/FourierSOPT.c
MPI_SRC = $(SRC_DIR)/FourierMPI.c

# Ejecutables
SECUENCIAL_EXE = $(BIN_DIR)/FourierS
SECUENCIAL_EXE2 = $(BIN_DIR)/FourierSOPT
MPI_EXE = $(BIN_DIR)/FourierMPI

all: $(SECUENCIAL_EXE) $(MPI_EXE) $(SECUENCIAL_EXE2)

$(SECUENCIAL_EXE2): $(SECUENCIAL_SRC2)
	$(CC) $(CFLAGS) $< -o $@ -lm

$(SECUENCIAL_EXE): $(SECUENCIAL_SRC)
	$(CC) $(CFLAGS) $< -o $@ -lm

$(MPI_EXE): $(MPI_SRC)
	$(MPICC) $(CFLAGS) $< -o $@ -lm

clean:
	rm -rf $(BIN_DIR)/*
