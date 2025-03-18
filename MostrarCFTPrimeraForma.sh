#!/bin/bash
echo "Secuencial CFT"
cat txt/ContinuoDFT.txt 
echo "Secuencial desenrollado CFT"
cat txt/ContinuoDFTOPT.txt 
echo "CFT en paralelo (MPI)"
cat txt/ContinuoDFTMPI.txt