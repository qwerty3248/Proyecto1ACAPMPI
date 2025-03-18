#!/bin/bash

echo "Secuencial CFT"
cat txt/ContinuoDFT2.txt 
echo "Secuencial desenrollado CFT"
cat txt/ContinuoDFTOPT2.txt 
echo "CFT en paralelo (MPI)"
cat txt/ContinuoDFTMPI2.txt