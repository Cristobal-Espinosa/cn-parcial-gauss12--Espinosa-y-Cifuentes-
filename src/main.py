import os
from gauss import gauss

carpeta = os.path.dirname(os.path.abspath(__file__))
raiz = os.path.dirname(carpeta)
ruta_a = os.path.join(raiz, "datos", "matriz_a.csv")
ruta_b = os.path.join(raiz, "datos", "vector_b.csv")
ruta_resultados = os.path.join(raiz, "resultados")

A = []
for linea in open(ruta_a):
    linea = linea.strip()
    if linea != "":
        fila = []
        for valor in linea.split(","):
            fila.append(float(valor))
        A.append(fila)

b = []
for linea in open(ruta_b):
    linea = linea.strip()
    if linea != "":
        b.append(float(linea))

print("Matriz leida:", len(A), "x", len(A[0]))
print("Vector b leido:", len(b), "valores")

sin = gauss(A, b, False)
con = gauss(A, b, True)

print()
print("Corrida sin pivoteo")
print("  multiplicador maximo:", sin["mult_max"])
print("  factor de crecimiento:", sin["crecimiento"])
print("  flops:", sin["contadores"]["mult_div"] + sin["contadores"]["sum_res"])

print()
print("Corrida con pivoteo parcial escalado")
print("  multiplicador maximo:", con["mult_max"])
print("  factor de crecimiento:", con["crecimiento"])
print("  flops:", con["contadores"]["mult_div"] + con["contadores"]["sum_res"])
print("  intercambios de fila:", con["contadores"]["intercambios"])

archivo = open(os.path.join(ruta_resultados, "traza_eliminacion.txt"), "w")
for nombre, corrida in (("SIN PIVOTEO", sin), ("PIVOTEO PARCIAL ESCALADO", con)):
    archivo.write("=== " + nombre + " ===\n")
    paso = 1
    for matriz in corrida["traza"]:
        archivo.write("\nMatriz despues del paso " + str(paso) + "\n")
        for fila in matriz:
            texto = ""
            for valor in fila:
                texto = texto + "%14.6e" % valor
            archivo.write(texto + "\n")
        paso = paso + 1
    archivo.write("\n")
archivo.close()

archivo = open(os.path.join(ruta_resultados, "solucion.csv"), "w")
archivo.write("incognita,valor\n")
i = 1
for valor in con["x"]:
    archivo.write("x" + str(i) + "," + repr(valor) + "\n")
    i = i + 1
archivo.close()

print()
print("Listo: se escribieron traza_eliminacion.txt y solucion.csv en resultados/")
