#!/bin/bash

echo "EL formato de los datos salidas es:"
echo "Número de elementos en el vector de muestras - tiempo que ha tardado en ejecutarlo (ms) "

echo ""
echo "DFT"

./MostrarDFT.sh 
echo ""
echo "CFT Primera Forma (Suma de rectangulos)"

./MostrarCFTPrimeraForma.sh
echo ""
echo "CFT Segunda Forma (Simpson)"

./MostrarCFT2.sh
echo ""
echo "CFT Tercera Forma (Trapezio)"

./MostrarCFT3.sh