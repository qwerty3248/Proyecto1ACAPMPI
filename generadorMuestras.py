import numpy as np

cota_abajo = -1
cota_arriba = 1


def generar_muestra(num_elementos, nombre_archivo):
    # Generar una señal de ejemplo (en este caso una señal aleatoria)
    muestras = np.random.uniform(cota_abajo,cota_arriba,num_elementos)  # Generamos una señal aleatoria de tamaño 'num_elementos'
    
    # Guardar en un archivo
    with open(nombre_archivo, 'w') as file:
        file.write(f"{num_elementos}\n")  # Escribir el número de elementos
        for muestra in muestras:
            file.write(f"{muestra}\n")  # Escribir cada muestra en una nueva línea

# Usar la función
num_elementos = 500000  # Cambia este número según lo que necesites
nombre_archivo = "txt/MuestraGenerada500000.txt"
generar_muestra(num_elementos, nombre_archivo)
