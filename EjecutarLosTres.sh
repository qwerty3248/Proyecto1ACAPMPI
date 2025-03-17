#!/bin/bash


./bin/FourierS && ./bin/FourierSOPT && mpirun -np $1 ./bin/FourierMPI